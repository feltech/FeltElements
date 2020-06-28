#include "Solver.hpp"
#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
#endif

#include "Attributes.hpp"
#include "Derivatives.hpp"
#include "Format.hpp"
#include <spdlog/spdlog.h>

namespace FeltElements::Solver
{
namespace
{
auto constexpr const index_of = [](auto const& haystack, auto&& needle) {
	auto const& it =
		std::find(haystack.cbegin(), haystack.cend(), std::forward<decltype(needle)>(needle));
	return static_cast<FeltElements::Tensor::Index>(std::distance(haystack.cbegin(), it));
};

constexpr Scalar epsilon = 0.00001;
}  // namespace

void update_elements_stiffness_and_internal_forces(
	Mesh const& mesh, FeltElements::Attributes& attributes, Scalar const lambda, Scalar const mu)
{
	for (auto itcellh = mesh.cells_begin(); itcellh != mesh.cells_end(); itcellh++)
	{
		auto const& cellh = *itcellh;
		auto const& cell_vtxhs = attributes.vtxh[cellh];
		auto [K, T] = Derivatives::KT(
			attributes.x.for_element(cell_vtxhs), attributes.dN_by_dX[cellh], lambda, mu);
		attributes.T[*itcellh] = T;
		attributes.K[*itcellh] = K;
	}
}
namespace LDLT
{
std::size_t solve(Mesh& mesh, Attributes& attrib, std::size_t max_steps, Scalar lambda, Scalar mu)
{
	EigenFixedDOFs const& mat_fixed_dof = ([&attrib,
											rows = static_cast<Eigen::Index>(mesh.n_vertices()),
											cols = static_cast<Eigen::Index>(Node::dim)]() {
		using EigenMapVectorFixedDOFs = Eigen::Map<Eigen::VectorXd>;
		EigenMapTensorVertices const& map{attrib.fixed_dof[Vtxh{0}].data(), rows, cols};
		// Copy to remove per-row alignment.
		EigenFixedDOFs mat = EigenMapVectorFixedDOFs{VerticesMatrix{map}.data(), rows * cols};
		return mat;
	})();

	std::size_t step;
	for (step = 0; step < max_steps; step++)
	{
		SPDLOG_DEBUG("LDLT iteration {}", step);
		Solver::update_elements_stiffness_and_internal_forces(mesh, attrib, lambda, mu);

		Eigen::VectorXd mat_T{3 * mesh.n_vertices()};
		mat_T.setZero();
		for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
		{
			auto const vtx_idx = (*itvtxh).idx();
			for (auto itcellh = mesh.vc_iter(*itvtxh); itcellh.valid(); itcellh++)
			{
				auto const& cell_vtxhs = attrib.vtxh[*itcellh];
				auto const& cell_vtx_idx = index_of(cell_vtxhs, *itvtxh);
				auto const& cell_T = attrib.T[*itcellh];

				mat_T.block<3, 1>(3 * vtx_idx, 0) +=
					EigenConstTensorMap<4, 3>{cell_T.data()}.block<1, 3>(
						static_cast<Eigen::Index>(cell_vtx_idx), 0);
			}
		}
		mat_T.array() *= (Eigen::VectorXd::Ones(mat_fixed_dof.size()) - mat_fixed_dof).array();

		Eigen::MatrixXd mat_K{3 * mesh.n_vertices(), 3 * mesh.n_vertices()};
		mat_K.setZero();

		auto const update_submatrix =
			[&mat_K, &attrib](auto const vtxh_src, auto const vtxh_dst, auto const cellh) {
				auto const& cell_K = attrib.K[cellh];
				auto const& cell_vtxhs = attrib.vtxh[cellh];
				auto const cell_a = index_of(cell_vtxhs, vtxh_src);
				auto const cell_b = index_of(cell_vtxhs, vtxh_dst);
				auto const a = vtxh_src.idx();
				auto const b = vtxh_dst.idx();

				using Tensor::Func::all;
				Tensor::Matrix<3> const Kab = cell_K(cell_a, all, cell_b, all);

				mat_K.block<3, 3>(3 * a, 3 * b) += EigenConstTensorMap<3, 3>{Kab.data()};
			};

		for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
		{
			for (auto itcellh = mesh.vc_iter(*itvtxh); itcellh.valid(); itcellh++)
			{ update_submatrix(*itvtxh, *itvtxh, *itcellh); }
		}
		for (auto itheh = mesh.halfedges_begin(); itheh != mesh.halfedges_end(); itheh++)
		{
			auto const& halfedge = mesh.halfedge(*itheh);
			auto const& vtxh_src = halfedge.from_vertex();
			auto const& vtxh_dst = halfedge.to_vertex();
			for (auto itcellh = mesh.hec_iter(*itheh); itcellh.valid(); itcellh++)
			{ update_submatrix(vtxh_src, vtxh_dst, *itcellh); }
		}

		SPDLOG_DEBUG("T\n{}", mat_T);
		mat_K += 10e5 * mat_fixed_dof.asDiagonal();
		SPDLOG_DEBUG("K (constrained)\n{}", mat_K);

		assert(std::abs(mat_K.determinant()) > 0.00001);

		Eigen::VectorXd mat_u{3 * mesh.n_vertices()};
		mat_u.setZero();
		mat_u = mat_K.ldlt().solve(-mat_T);
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
}  // namespace FeltElements::Solver
