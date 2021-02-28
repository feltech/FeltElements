#include "Solver.hpp"

#include <spdlog/spdlog.h>
#include <Eigen/IterativeLinearSolvers>
#include <boost/range/irange.hpp>
#ifdef NDEBUG
#include <boost/preprocessor/stringize.hpp>
#endif

#include "Attributes.hpp"
#include "Derivatives.hpp"
#include "MeshFacade.hpp"
#include "internal/Format.hpp"	// For logging

#ifdef NDEBUG
#define FE_PARALLEL_PREAMBLE omp parallel for default(none)
#define FE_PARALLEL(args) _Pragma(BOOST_PP_STRINGIZE(FE_PARALLEL_PREAMBLE args))
#else
#define FE_PARALLEL(args)
#endif

namespace FeltElements::Solver
{
void Base::update_elements_stiffness_and_residual()
{
	auto const num_cells = static_cast<int>(m_mesh.n_cells());

	FE_PARALLEL(shared(num_cells, m_attrs))
	for (int cell_idx = 0; cell_idx < num_cells; cell_idx++)
	{
		Cellh cellh{cell_idx};
		auto const & cell_vtxhs = m_attrs.vtxhs[cellh];
		auto const & boundary_faces_vtxh_idxs = m_attrs.boundary_faces_vtxh_idxs[cellh];

		auto const & [K, R] = Derivatives::KR(
			m_attrs.x.for_element(cell_vtxhs),
			boundary_faces_vtxh_idxs,
			m_attrs.x.for_elements(cell_vtxhs, boundary_faces_vtxh_idxs),
			m_attrs.dN_by_dX[cellh],
			*m_attrs.material,
			*m_attrs.forces);
		m_attrs.K[cellh] = K;
		m_attrs.R[cellh] = R;
	}
}

Scalar Base::find_epsilon() const
{
	Scalar total_V = 0;
	Scalar min_V = std::numeric_limits<Scalar>::max();
	Scalar max_V = 0;
	Element::NodePositions min_X;

	for (auto const & X : MeshIters{m_mesh, m_attrs}.Xs())
	{
		auto const V = Derivatives::V(X);
		total_V += V;
		if (V > max_V)
			max_V = V;
		if (V < min_V)
		{
			min_X = X;
			min_V = V;
		}
	}
	Scalar const avg_V = total_V / static_cast<Scalar>(m_mesh.n_cells());
	// Edge length assuming a regular tetrahedron, divided by a constant.
	Scalar const epsilon = std::pow(6 * std::sqrt(2.0) * min_V, 1.0 / 3.0) * 1e-5;
	spdlog::info(
		"total V = {}; mean V = {}; max V = {}, min V = {} @\n{}\nepsilon = {}",
		total_V,
		avg_V,
		max_V,
		min_V,
		min_X,
		epsilon);
	return epsilon;
}

namespace
{
void log_xs(Mesh const & mesh, Attributes const & attrs)
{
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_DEBUG
	std::string xs_str;
	for (auto const & x : MeshIters{mesh, attrs}.xs()) xs_str += fmt::format("{}\n", x);

	spdlog::debug(xs_str);
#else
	(void)mesh;
	(void)attrs;
#endif
}
}  // namespace

void Matrix::solve()
{
	Scalar const epsilon = find_epsilon();
	//	constexpr Scalar penalty = std::numeric_limits<Scalar>::max() / 100;

	EigenFixedDOFs const & mat_fixed_dof = ([&attrs = m_attrs,
											 rows = static_cast<Eigen::Index>(m_mesh.n_vertices()),
											 cols = static_cast<Eigen::Index>(Node::dim)]() {
		EigenMapTensorVertices const & map{attrs.fixed_dof[Vtxh{0}].data(), rows, cols};
		// Copy to remove per-row alignment.
		EigenFixedDOFs vec = Eigen::Map<EigenFixedDOFs>{VerticesMatrix{map}.data(), rows * cols};
		return vec;
	})();

	Eigen::VectorXd mat_R{3 * m_mesh.n_vertices()};
	Eigen::MatrixXd mat_K{3 * m_mesh.n_vertices(), 3 * m_mesh.n_vertices()};
	Eigen::VectorXd mat_u{3 * m_mesh.n_vertices()};

	Node::Force const force_increment = m_attrs.forces->F_by_m / m_params.num_force_increments;
	m_attrs.forces->F_by_m = 0;

	Scalar max_norm;

	for (std::size_t increment_num = 0; increment_num < m_params.num_force_increments;
		 increment_num++)
	{
		m_attrs.forces->F_by_m += force_increment;
		stats.force_increment_counter++;

		for (std::size_t step = 0; step < m_params.num_steps; step++)
		{
			stats.step_counter++;
			SPDLOG_DEBUG("Matrix solver iteration {}:{}", increment_num, step);

			update_elements_stiffness_and_residual();

			mat_R.setZero();
			mat_K.setZero();

			for (auto vtxh : boost::make_iterator_range(m_mesh.vertices()))
			{
				auto const vtx_idx = vtxh.idx();
				for (auto cellh : boost::make_iterator_range(m_mesh.vertex_cells(vtxh)))
				{
					auto const & cell_vtxhs = m_attrs.vtxhs[cellh];
					auto const & cell_vtx_idx = index_of<Eigen::Index>(cell_vtxhs, vtxh);
					auto const & cell_R = m_attrs.R[cellh];

					// Forces for node `a`.
					auto const & Ra =
						EigenConstTensorMap<4, 3>{cell_R.data()}.block<1, 3>(cell_vtx_idx, 0);

					// Update global residual with difference between internal and external forces.
					mat_R.block<3, 1>(3 * vtx_idx, 0) += Ra;
				}
			}

			// Compute sum to check equilibrium.
			//		Eigen::Vector3d sum; sum.setZero();
			//		Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor>>
			// mat_R_Nx3{mat_R.data(), mat_R.rows() / 3, 3}; 		for (Eigen::Index const node_idx
			// : boost::irange(mat_R_Nx3.rows())) 			sum += mat_R_Nx3.row(node_idx);

			auto const update_submatrix =
				[&mat_K, &attrs = m_attrs](
					auto const vtxh_src, auto const vtxh_dst, auto const cellh_) {
					auto const & cell_K = attrs.K[cellh_];
					auto const & cell_vtxhs = attrs.vtxhs[cellh_];
					auto const cell_a = index_of(cell_vtxhs, vtxh_src);
					auto const cell_b = index_of(cell_vtxhs, vtxh_dst);
					auto const a = vtxh_src.idx();
					auto const b = vtxh_dst.idx();

					using Tensor::Func::all;
					Tensor::Matrix<3> const Kab = cell_K(cell_a, all, cell_b, all);

					mat_K.block<3, 3>(3 * a, 3 * b) += EigenConstTensorMap<3, 3>{Kab.data()};
				};

			for (auto vtxh : boost::make_iterator_range(m_mesh.vertices()))
			{
				for (auto cellh : boost::make_iterator_range(m_mesh.vertex_cells(vtxh)))
					update_submatrix(vtxh, vtxh, cellh);
			}
			for (auto heh : boost::make_iterator_range(m_mesh.halfedges()))
			{
				auto const & halfedge = m_mesh.halfedge(heh);
				auto const & vtxh_src = halfedge.from_vertex();
				auto const & vtxh_dst = halfedge.to_vertex();
				for (auto cellh : boost::make_iterator_range(m_mesh.halfedge_cells(heh)))
					update_submatrix(vtxh_src, vtxh_dst, cellh);
			}
			// Inspect matrix to calculate penalty of relative size.
			Scalar const penalty = mat_K.lpNorm<Eigen::Infinity>() * 10000;
			// Zero-out penalised degrees of freedom.
			mat_K.diagonal() = mat_K.diagonal().cwiseProduct(
				Eigen::VectorXd::Ones(mat_fixed_dof.size()) - mat_fixed_dof);
			// Set penalised degrees of freedom to penalty value.
			mat_K += penalty * mat_fixed_dof.asDiagonal();

			//			SPDLOG_DEBUG("K (constrained)\n{}", mat_K);
			//			auto const detK = mat_K.determinant();
			//			if (std::abs(detK) < 0.00001)
			//				throw std::invalid_argument("Stiffness matrix |K| ~ 0");
			//		if (!std::isfinite(detK))
			//			throw std::invalid_argument{fmt::format(
			//				"Stiffness matrix |K| = {} with max = {}", detK,

			// mat_K.lpNorm<Eigen::Infinity>())};
			//		mat_u = mat_K.ldlt().solve(-mat_R);
			//				Eigen::ConjugateGradient< decltype(mat_K), Eigen::Lower|Eigen::Upper>
			// cg; cg.compute(mat_K); 				mat_u = cg.solve(-mat_R);
			mat_u = mat_K.partialPivLu().solve(-mat_R);

			//			if (!mat_u.array().isFinite().all())
			//				throw std::logic_error(
			//					fmt::format("Resulting displacement is not finite u = \n{}.",
			// mat_u));

			mat_u.array() *= (Eigen::VectorXd::Ones(mat_u.size()) - mat_fixed_dof).array();
			//			SPDLOG_DEBUG("u (constrained)\n{}", mat_u);

			for (auto const & vtxh : boost::make_iterator_range(m_mesh.vertices()))
				m_attrs.x[vtxh] += Tensor::Map<3>{mat_u.block<3, 1>(3 * vtxh.idx(), 0).data()};

			max_norm = mat_u.lpNorm<Eigen::Infinity>();
			stats.max_norm = max_norm;
			log_xs(m_mesh, m_attrs);
			if (max_norm < epsilon)
				break;
		}
	}
}

void Gauss::solve()
{
	Scalar const epsilon = find_epsilon() / 1000;
	using Tensor::Func::all;

	std::vector<Node::Pos> u(m_mesh.n_vertices(), {0, 0, 0});
	Scalar max_norm = 0;

	Node::Force const force_increment = m_attrs.forces->F_by_m / m_params.num_force_increments;
	m_attrs.forces->F_by_m = 0;

	for (std::size_t increment_num = 0; increment_num < m_params.num_force_increments;
		 increment_num++)
	{
		m_attrs.forces->F_by_m += force_increment;
		stats.force_increment_counter++;

		for (std::size_t step = 0; step < m_params.num_steps; step++)
		{
			stats.step_counter++;
			pauser.wait_while_paused();

			SPDLOG_DEBUG("Gauss iteration {}:{}", increment_num, step);

			update_elements_stiffness_and_residual();

			max_norm = 0;
			for (auto & u_a : u) u_a *= 0.99;  //.zeros();

			FE_PARALLEL(shared(m_mesh, m_attrs, u, force_increment) reduction(max : max_norm))
			for (Tensor::Index a = 0; a < m_mesh.n_vertices(); a++)
			{
				Vtxh vtxh_src{static_cast<int>(a)};
				auto const & fixed_dof = m_attrs.fixed_dof[vtxh_src];
				if (Tensor::Func::all_of(fixed_dof != 0))
					continue;

				Node::Force Ra = 0;
				Tensor::Matrix<3> Kaa = 0;
				Node::Force Ka_u = 0;
				for (auto cellh : boost::make_iterator_range(m_mesh.vertex_cells(vtxh_src)))
				{
					auto const & cell_vtxhs = m_attrs.vtxhs[cellh];
					auto const & cell_a = index_of(cell_vtxhs, vtxh_src);
					auto const & cell_R = m_attrs.R[cellh];
					auto const & cell_K = m_attrs.K[cellh];

					Ra += cell_R(cell_a, all);
					Kaa += cell_K(cell_a, all, cell_a, all);
				}
				//			diag(Kaa) += penalty * fixed_dof;

				for (auto heh : boost::make_iterator_range(m_mesh.outgoing_halfedges(vtxh_src)))
				{
					Tensor::Multi<Node::dim, Node::dim> Kab = 0;
					auto const & halfedge = m_mesh.halfedge(heh);
					auto const & vtxh_dst = halfedge.to_vertex();
					auto const b = static_cast<Tensor::Index>(vtxh_dst.idx());

					for (auto cellh : boost::make_iterator_range(m_mesh.halfedge_cells(heh)))
					{
						auto const & cell_K = m_attrs.K[cellh];
						auto const & cell_vtxhs = m_attrs.vtxhs[cellh];
						auto const cell_a = index_of(cell_vtxhs, vtxh_src);
						auto const cell_b = index_of(cell_vtxhs, vtxh_dst);
						Kab += cell_K(cell_a, all, cell_b, all);
					}
					Ka_u += Kab % u[b];
				}

				using Tensor::Func::inv;
				u[a] = inv(Kaa) % (-Ra - Ka_u) * (1.0 - fixed_dof);
				m_attrs.x[vtxh_src] += u[a];
				using Tensor::Func::abs;
				using Tensor::Func::max;
				// Note: double-reduction in case OpenMP disabled.
				max_norm = std::max(max_norm, max(abs(u[a])));
				//				SPDLOG_DEBUG("R[{}] = {}", a, Ra);
				//				SPDLOG_DEBUG("u[{}] = {}", a, u[a]);
			}

			SPDLOG_DEBUG("Max norm: {}", max_norm);
			stats.max_norm = max_norm;

			log_xs(m_mesh, m_attrs);

			if (max_norm < epsilon)
				break;
		}
	}
}
}  // namespace FeltElements::Solver
