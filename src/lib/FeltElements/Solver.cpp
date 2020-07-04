#include "Solver.hpp"
#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
#endif

#include <spdlog/spdlog.h>

#include "Attributes.hpp"
#include "Derivatives.hpp"
#include "Format.hpp"

namespace FeltElements::Solver
{
namespace
{
auto constexpr const index_of = [](auto const & haystack, auto && needle) {
	auto const & it =
		std::find(haystack.cbegin(), haystack.cend(), std::forward<decltype(needle)>(needle));
	return static_cast<FeltElements::Tensor::Index>(std::distance(haystack.cbegin(), it));
};

constexpr Scalar penalty = 10e5;
}  // namespace

void update_elements_stiffness_and_forces(Mesh const & mesh, FeltElements::Attributes & attributes)
{
	for (auto cellh : boost::make_iterator_range(mesh.cells()))
	{
		auto const & cell_vtxhs = attributes.vtxh[cellh];
		auto [K, T, v] = Derivatives::KTv(
			attributes.x.for_element(cell_vtxhs),
			attributes.dN_by_dX[cellh],
			attributes.V[cellh],
			(*attributes.material).lambda,
			(*attributes.material).mu);
		attributes.v[cellh] = v;
		attributes.T[cellh] = T;
		attributes.K[cellh] = K;
	}
}

namespace LDLT
{
std::size_t solve(Mesh & mesh, Attributes & attrib, std::size_t max_steps)
{
	constexpr Scalar epsilon = 0.00005;

	EigenFixedDOFs const & mat_fixed_dof = ([&attrib,
											 rows = static_cast<Eigen::Index>(mesh.n_vertices()),
											 cols = static_cast<Eigen::Index>(Node::dim)]() {
		using EigenMapVectorFixedDOFs = Eigen::Map<Eigen::VectorXd>;
		EigenMapTensorVertices const & map{attrib.fixed_dof[Vtxh{0}].data(), rows, cols};
		// Copy to remove per-row alignment.
		EigenFixedDOFs mat = EigenMapVectorFixedDOFs{VerticesMatrix{map}.data(), rows * cols};
		return mat;
	})();

	Scalar const rho = attrib.material->rho;
	Eigen::Vector3d const & F = 1.0 / 4.0 * rho * EigenConstTensorMap<1, 3>{(*attrib.f).data()};

	Eigen::VectorXd mat_R{3 * mesh.n_vertices()};
	Eigen::MatrixXd mat_K{3 * mesh.n_vertices(), 3 * mesh.n_vertices()};
	Eigen::VectorXd mat_u{3 * mesh.n_vertices()};

	std::size_t step;
	for (step = 0; step < max_steps; step++)
	{
		SPDLOG_DEBUG("LDLT iteration {}", step);
		Solver::update_elements_stiffness_and_forces(mesh, attrib);

		mat_R.setZero();
		mat_K.setZero();

		for (auto vtxh : boost::make_iterator_range(mesh.vertices()))
		{
			auto const vtx_idx = vtxh.idx();
			for (auto cellh : boost::make_iterator_range(mesh.vertex_cells(vtxh)))
			{
				auto const & cell_vtxhs = attrib.vtxh[cellh];
				auto const & cell_vtx_idx = index_of(cell_vtxhs, vtxh);
				auto const & cell_T = attrib.T[cellh];

				auto const & Ta = EigenConstTensorMap<4, 3>{cell_T.data()}.block<1, 3>(
					static_cast<Eigen::Index>(cell_vtx_idx), 0);

				auto const & Fa = attrib.v[cellh] * F.transpose();

				mat_R.block<3, 1>(3 * vtx_idx, 0) += Ta - Fa;
			}
		}

		auto const update_submatrix =
			[&mat_K, &attrib](auto const vtxh_src, auto const vtxh_dst, auto const cellh_) {
				auto const & cell_K = attrib.K[cellh_];
				auto const & cell_vtxhs = attrib.vtxh[cellh_];
				auto const cell_a = index_of(cell_vtxhs, vtxh_src);
				auto const cell_b = index_of(cell_vtxhs, vtxh_dst);
				auto const a = vtxh_src.idx();
				auto const b = vtxh_dst.idx();

				using Tensor::Func::all;
				Tensor::Matrix<3> const Kab = cell_K(cell_a, all, cell_b, all);

				mat_K.block<3, 3>(3 * a, 3 * b) += EigenConstTensorMap<3, 3>{Kab.data()};
			};

		for (auto vtxh : boost::make_iterator_range(mesh.vertices()))
		{
			for (auto cellh : boost::make_iterator_range(mesh.vertex_cells(vtxh)))
				update_submatrix(vtxh, vtxh, cellh);
		}
		for (auto heh : boost::make_iterator_range(mesh.halfedges()))
		{
			auto const & halfedge = mesh.halfedge(heh);
			auto const & vtxh_src = halfedge.from_vertex();
			auto const & vtxh_dst = halfedge.to_vertex();
			for (auto cellh : boost::make_iterator_range(mesh.halfedge_cells(heh)))
				update_submatrix(vtxh_src, vtxh_dst, cellh);
		}

		//		mat_R.array() *= (Eigen::VectorXd::Ones(mat_R.size()) - mat_fixed_dof).array();
		SPDLOG_DEBUG("R\n{}", mat_R);
		mat_K += penalty * mat_fixed_dof.asDiagonal();
		SPDLOG_DEBUG("K (constrained)\n{}", mat_K);

		//		assert(std::abs(mat_K.determinant()) > 0.00001);

		mat_u = mat_K.ldlt().solve(-mat_R);
		mat_u.array() *= (Eigen::VectorXd::Ones(mat_u.size()) - mat_fixed_dof).array();
		SPDLOG_DEBUG("u (constrained)\n{}", mat_u);

		for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
		{
			auto const vtx_idx = (*itvtxh).idx();
			attrib.x[*itvtxh] += Tensor::Map<3>{mat_u.block<3, 1>(3 * vtx_idx, 0).data()};
		}

		if (mat_u.lpNorm<Eigen::Infinity>() < epsilon)
			break;
	}

	return step;
}
}  // namespace LDLT

namespace Gauss
{
std::size_t solve(Mesh & mesh, Attributes & attrib, std::size_t max_steps)
{
	constexpr Scalar epsilon = 0.00001;
	using Tensor::Func::all;

	Scalar const rho = attrib.material->rho;
	Node::Force const & F = 1.0 / 4.0 * rho * (*attrib.f);

	std::vector<Node::Pos> u{};
	u.resize(mesh.n_vertices());

	std::size_t step;
	for (step = 0; step < max_steps; step++)
	{
		SPDLOG_DEBUG("Gauss iteration {}", step);

		Solver::update_elements_stiffness_and_forces(mesh, attrib);

		for (auto & u_a : u) u_a.zeros();

		Scalar max_norm = 0;

		for (auto vtxh_src : boost::make_iterator_range(mesh.vertices()))
		{
			auto const a = static_cast<Tensor::Index>(vtxh_src.idx());
			auto const & fixed_dof = attrib.fixed_dof[vtxh_src];

			Node::Force Ra = 0;
			Tensor::Matrix<3> Kaa = 0;
			Node::Force Ka_u = 0;
			for (auto cellh : boost::make_iterator_range(mesh.vertex_cells(vtxh_src)))
			{
				auto const & cell_vtxhs = attrib.vtxh[cellh];
				auto const & cell_a = index_of(cell_vtxhs, vtxh_src);
				auto const & cell_T = attrib.T[cellh];
				auto const & cell_K = attrib.K[cellh];

				auto const & Fa = attrib.v[cellh] * F;

				Ra += cell_T(cell_a, all) - Fa;
				Kaa += cell_K(cell_a, all, cell_a, all);
			}
			diag(Kaa) += penalty * fixed_dof;

			for (auto heh : boost::make_iterator_range(mesh.outgoing_halfedges(vtxh_src)))
			{
				Tensor::Multi<Node::dim, Node::dim> Kab = 0;
				auto const & halfedge = mesh.halfedge(heh);
				auto const & vtxh_dst = halfedge.to_vertex();
				auto const b = static_cast<Tensor::Index>(vtxh_dst.idx());

				for (auto cellh : boost::make_iterator_range(mesh.halfedge_cells(heh)))
				{
					auto const & cell_K = attrib.K[cellh];
					auto const & cell_vtxhs = attrib.vtxh[cellh];
					auto const cell_a = index_of(cell_vtxhs, vtxh_src);
					auto const cell_b = index_of(cell_vtxhs, vtxh_dst);
					Kab += cell_K(cell_a, all, cell_b, all);
				}
				Ka_u += Kab % u[b];
			}

			using Tensor::Func::inv;
			u[a] = inv(Kaa) % (-Ra - Ka_u) * (1.0 - fixed_dof);
			attrib.x[vtxh_src] += u[a];
			using Tensor::Func::abs;
			using Tensor::Func::max;
			max_norm = std::max(max(abs(u[a])), max_norm);
			SPDLOG_DEBUG("\nR[{}] = {}", a, Ra);
			SPDLOG_DEBUG("\nu[{}] = {}", a, u[a]);
		}

		if (max_norm < epsilon)
			break;
	}

	return step;
}
}  // namespace Gauss
}  // namespace FeltElements::Solver
