// For spdlog
#include "internal/Format.hpp"
#include "Solver.hpp"
#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
#endif

#include <spdlog/spdlog.h>
#include <boost/range/irange.hpp>

#include "Attributes.hpp"
#include "Derivatives.hpp"

namespace FeltElements::Solver
{
void update_elements_stiffness_and_residual(
	Mesh const & mesh, FeltElements::Attributes & attributes)
{
	for (auto cellh : boost::make_iterator_range(mesh.cells()))
	{
		auto const & cell_vtxhs = attributes.vtxhs[cellh];
		auto const & boundary_faces_vtxh_idxs = attributes.boundary_faces_vtxh_idxs[cellh];

		auto const & [K, R] = Derivatives::KR(
			attributes.x.for_element(cell_vtxhs),
			boundary_faces_vtxh_idxs,
			attributes.x.for_elements(cell_vtxhs, boundary_faces_vtxh_idxs),
			attributes.dN_by_dX[cellh],
			*attributes.material,
			*attributes.forces);
		attributes.K[cellh] = K;
		attributes.R[cellh] = R;
	}
}

namespace LDLT
{
std::size_t solve(Mesh & mesh, Attributes & attrib, std::size_t max_steps)
{
	constexpr Scalar epsilon = 0.00005;
	constexpr Scalar penalty = std::numeric_limits<Scalar>::max() / 2;

	EigenFixedDOFs const & mat_fixed_dof = ([&attrib,
											 rows = static_cast<Eigen::Index>(mesh.n_vertices()),
											 cols = static_cast<Eigen::Index>(Node::dim)]() {
		EigenMapTensorVertices const & map{attrib.fixed_dof[Vtxh{0}].data(), rows, cols};
		// Copy to remove per-row alignment.
		EigenFixedDOFs vec = Eigen::Map<EigenFixedDOFs>{VerticesMatrix{map}.data(), rows * cols};
		return vec;
	})();

	Eigen::VectorXd mat_R{3 * mesh.n_vertices()};
	Eigen::MatrixXd mat_K{3 * mesh.n_vertices(), 3 * mesh.n_vertices()};
	Eigen::VectorXd mat_u{3 * mesh.n_vertices()};

	std::size_t step;
	for (step = 0; step < max_steps; step++)
	{
		SPDLOG_DEBUG("LDLT iteration {}", step);
		Solver::update_elements_stiffness_and_residual(mesh, attrib);

		mat_R.setZero();
		mat_K.setZero();

		for (auto vtxh : boost::make_iterator_range(mesh.vertices()))
		{
			auto const vtx_idx = vtxh.idx();
			for (auto cellh : boost::make_iterator_range(mesh.vertex_cells(vtxh)))
			{
				auto const & cell_vtxhs = attrib.vtxhs[cellh];
				auto const & cell_vtx_idx = index_of(cell_vtxhs, vtxh);
				auto const & cell_R = attrib.R[cellh];

				// Forces for node `a`.
				auto const & Ra = EigenConstTensorMap<4, 3>{cell_R.data()}.block<1, 3>(
					static_cast<Eigen::Index>(cell_vtx_idx), 0);

				// Update global residual with difference between internal and external forces.
				mat_R.block<3, 1>(3 * vtx_idx, 0) += Ra;
			}
		}

		// Compute sum to check equilibrium.
//		Eigen::Vector3d sum; sum.setZero();
//		Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor>> mat_R_Nx3{mat_R.data(), mat_R.rows() / 3, 3};
//		for (Eigen::Index const node_idx : boost::irange(mat_R_Nx3.rows()))
//			sum += mat_R_Nx3.row(node_idx);

		auto const update_submatrix =
			[&mat_K, &attrib](auto const vtxh_src, auto const vtxh_dst, auto const cellh_) {
				auto const & cell_K = attrib.K[cellh_];
				auto const & cell_vtxhs = attrib.vtxhs[cellh_];
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

		assert(std::abs(mat_K.determinant()) > 0.00001);

		mat_u = mat_K.ldlt().solve(-mat_R);
		mat_u.array() *= (Eigen::VectorXd::Ones(mat_u.size()) - mat_fixed_dof).array();
		SPDLOG_DEBUG("u (constrained)\n{}", mat_u);

		for (auto const & vtxh : boost::make_iterator_range(mesh.vertices()))
		{
			attrib.x[vtxh] += Tensor::Map<3>{mat_u.block<3, 1>(3 * vtxh.idx(), 0).data()};
			//			SPDLOG_DEBUG("({}, {}, {})\n", attrib.x[vtxh](0), attrib.x[vtxh](1),
			// attrib.x[vtxh](2));
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

	std::vector<Node::Pos> u{};
	u.resize(mesh.n_vertices());

	std::size_t step;
	for (step = 0; step < max_steps; step++)
	{
		SPDLOG_DEBUG("Gauss iteration {}", step);

		Solver::update_elements_stiffness_and_residual(mesh, attrib);

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
				auto const & cell_vtxhs = attrib.vtxhs[cellh];
				auto const & cell_a = index_of(cell_vtxhs, vtxh_src);
				auto const & cell_R = attrib.R[cellh];
				auto const & cell_K = attrib.K[cellh];

				Ra += cell_R(cell_a, all);
				Kaa += cell_K(cell_a, all, cell_a, all);
			}
			//			diag(Kaa) += penalty * fixed_dof;

			for (auto heh : boost::make_iterator_range(mesh.outgoing_halfedges(vtxh_src)))
			{
				Tensor::Multi<Node::dim, Node::dim> Kab = 0;
				auto const & halfedge = mesh.halfedge(heh);
				auto const & vtxh_dst = halfedge.to_vertex();
				auto const b = static_cast<Tensor::Index>(vtxh_dst.idx());

				for (auto cellh : boost::make_iterator_range(mesh.halfedge_cells(heh)))
				{
					auto const & cell_K = attrib.K[cellh];
					auto const & cell_vtxhs = attrib.vtxhs[cellh];
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
			SPDLOG_DEBUG("R[{}] = {}", a, Ra);
			SPDLOG_DEBUG("u[{}] = {}", a, u[a]);
		}

		if (max_norm < epsilon)
			break;
	}

	return step;
}
}  // namespace Gauss
}  // namespace FeltElements::Solver
