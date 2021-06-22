#include "Solver.hpp"

#include <iostream>

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
void Base::update_elements_stiffness_and_residual(Scalar const lambda)
{
	auto const num_cells = static_cast<int>(m_mesh.n_cells());

	FE_PARALLEL(shared(num_cells, m_attrs, lambda))
	for (int cell_idx = 0; cell_idx < num_cells; cell_idx++)
	{
		Cellh cellh{cell_idx};
		auto const & cell_vtxhs = m_attrs.vtxhs[cellh];
		auto const & boundary_faces_vtxh_idxs = m_attrs.boundary_faces_vtxh_idxs[cellh];

		auto const & [K, F, R] = Derivatives::KR(
			m_attrs.x.for_element(cell_vtxhs),
			boundary_faces_vtxh_idxs,
			m_attrs.x.for_elements(cell_vtxhs, boundary_faces_vtxh_idxs),
			m_attrs.dN_by_dX[cellh],
			*m_attrs.material,
			*m_attrs.forces,
			lambda);
		m_attrs.F[cellh] = F;
		m_attrs.R[cellh] = R;
		m_attrs.K[cellh] = K;
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
	Scalar const epsilon = std::pow(6 * std::sqrt(2.0) * min_V, 1.0 / 3.0);
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
	//#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_DEBUG
	std::string xs_str;
	for (auto const & x : MeshIters{mesh, attrs}.xs()) xs_str += fmt::format("{}\n", x);

	spdlog::info(xs_str);
	//	spdlog::debug(xs_str);
	//#else
	//	(void)mesh;
	//	(void)attrs;
	//#endif
}

template <typename T>
int sgn(T val)
{
	return (T(0) < val) - (val < T(0));
}
}  // namespace

void Matrix::solve()
{
	std::ofstream numbers("feltelements.tsv");
	Scalar epsilon = find_epsilon();

	Eigen::VectorXd const & mat_fixed_dof = ([&attrs = m_attrs,
											  rows = static_cast<Eigen::Index>(m_mesh.n_vertices()),
											  cols = static_cast<Eigen::Index>(Node::dim)]() {
		EigenMapTensorVerticesConst const & map{attrs.fixed_dof[Vtxh{0}].data(), rows, cols};
		// Copy to remove per-row alignment.
		Eigen::VectorXd vec = Eigen::Map<Eigen::VectorXd>{VerticesMatrix{map}.data(), rows * cols};
		return vec;
	})();
	Eigen::VectorXd const one_minus_fixed_dof =
		Eigen::VectorXd::Ones(mat_fixed_dof.size()) - mat_fixed_dof;

	auto const num_vertices = static_cast<Eigen::Index>(m_mesh.n_vertices());
	auto const num_dofs = static_cast<Eigen::Index>(Node::dim) * num_vertices;

	Eigen::VectorXd vec_F{num_dofs};
	Eigen::VectorXd vec_R{num_dofs};
	Eigen::MatrixXd mat_K{num_dofs, num_dofs};
	Eigen::VectorXd vec_u{num_dofs};
	Eigen::VectorXd vec_uF{num_dofs};
	Eigen::VectorXd vec_uR{num_dofs};
	Eigen::VectorXd vec_delta_x{num_dofs};
	Scalar delta_lambda = 0.0;
	Scalar lambda = 0.01;
	Scalar gamma = 0;
	Scalar max_norm;
	constexpr std::size_t step_target = 30;
	std::size_t step = step_target;
	Scalar const s20 = 1e-13;
	Scalar s2 = s20;
	Scalar const psi2 = 0;
	epsilon = 1e-6;
	Scalar s2_mult = 1e3;
	Scalar s2_min;
	Scalar s2_max = 0.01;
	Scalar const gamma_max = 0.01;
//	Scalar const gamma_max = 0.01;
	auto mat_x = as_matrix(m_attrs.x);

	constexpr Scalar s2_search = 1e-2;

	VerticesMatrix mat_X = as_matrix(m_attrs.X);

	for (std::size_t increment_num = 0; increment_num < m_params.num_force_increments;
		 increment_num++)
	{
		stats.force_increment_counter++;

		s2 = (step < step_target) ? s2 * (1.0 + s2_search) : s2 * (1.0 - s2_search);
		VerticesMatrix const mat_x0 = mat_x;
		delta_lambda = 0;
		s2_mult = 1e6;

		Scalar const lambda0 = lambda;

		// Pre-increment solution
		{
			update_elements_stiffness_and_residual(lambda);
			assemble(mat_K, vec_R, vec_F, one_minus_fixed_dof);
			vec_F /= lambda;
			auto mat_K_LU = mat_K.partialPivLu();
			vec_uR = mat_K_LU.solve(-vec_R);
			vec_uF = mat_K_LU.solve(vec_F);

//			s2 = gamma_max * gamma_max * vec_uF.squaredNorm();


			Eigen::VectorXd const n = vec_uF.normalized();
			Eigen::VectorXd const vec_s_min = vec_uR - n * n.dot(vec_uR);
			s2_min = vec_s_min.squaredNorm() * 1.001;

			auto const calc_s2 = [&]
			{
				s2 = vec_uR.squaredNorm() * s2_mult;
				gamma = sqrt(s2 / vec_uF.squaredNorm());
			};

//			calc_s2();
			gamma = sqrt(s2 / vec_uF.squaredNorm());
			while (gamma > gamma_max && s2 > s2_min)
			{
				s2_mult /= 2;
				calc_s2();
			}
			s2 = std::clamp(s2, s2_min, s2_max);

			//			auto const det_K = mat_K_LU.determinant();
//			const auto choose_correction = [&]
//			{
//				auto [gamma1, gamma2] =
//					arc_length(vec_uF, vec_uR, vec_F, Eigen::VectorXd::Zero(vec_delta_x.size()), 0, s2, psi2);
//				Scalar const sign = increment_num == 0 ? 1.0 : sgn(vec_delta_x.dot(vec_uF));
//				gamma = (gamma1 * sign > 0) ? gamma1 : gamma2;
//			};
//			choose_correction();
//			while (std::abs(gamma) > gamma_max && s2 / 2 > s2_min)
//			{
//				s2 /= 2;
//				choose_correction();
//			}

			delta_lambda = gamma;
			lambda += gamma;
			vec_delta_x = vec_uR + gamma * vec_uF;
			vec_delta_x.array() *= one_minus_fixed_dof.array();
			mat_x += as_matrix(vec_delta_x);
		}
//		Eigen::VectorXd vec_delta_x0 = vec_delta_x;

		spdlog::info(
			"increment = {}; total steps = {}; lambda = {}; gamma = {}; s2 = {}",
			increment_num,
			stats.step_counter.load(),
			lambda,
			gamma,
			s2);

		for (step = 0; step < m_params.num_steps; ++step)
		{
			if (step > step_target * 30)
			{
				lambda = lambda0;
				mat_x = mat_x0;
				break;
			}
			++stats.step_counter;
//			log_xs(m_mesh, m_attrs);

			update_elements_stiffness_and_residual(lambda);
			assemble(mat_K, vec_R, vec_F, one_minus_fixed_dof);
			vec_F /= lambda;
			max_norm = vec_R.squaredNorm();

//			spdlog::info(
//				"increment = {}; step = {}; lambda = {}; delta_lambda = {}; gamma = {}; s2 = {}; norm = {}",
//				increment_num,
//				step,
//				lambda,
//				delta_lambda,
//				gamma,
//				s2,
//				max_norm);

			stats.max_norm = max_norm;
			if (max_norm < epsilon)
				break;

			auto const & mat_K_LU = mat_K.partialPivLu();
			vec_uR = mat_K_LU.solve(-vec_R);
			vec_uF = mat_K_LU.solve(vec_F);

			Eigen::VectorXd const delta_xu = vec_delta_x + vec_uR;
			Eigen::VectorXd const n = vec_uF.normalized();
			Eigen::VectorXd const vec_s_min = delta_xu - n * n.dot(delta_xu);
			s2_min = vec_s_min.squaredNorm() * 1.001;

//			s2_min = (vec_delta_x + vec_uR).squaredNorm() * 1.001;
			if (s2 < s2_min)
			{
				s2 = s2_min;
//				lambda = lambda0;
//				mat_x = mat_x0;
//				step = 0;
//				break;
			}

			Scalar gamma1, gamma2;
			Eigen::VectorXd vec_delta_x1, vec_delta_x2;

			auto const choose_correction = [&]
			{
				std::tie(gamma1, gamma2) =
					arc_length(vec_uF, vec_uR, vec_F, vec_delta_x, delta_lambda, s2, psi2);

				auto const vec_u1 = vec_uR + gamma1 * vec_uF;
				auto const vec_u2 = vec_uR + gamma2 * vec_uF;
				vec_delta_x1 = vec_delta_x + vec_u1;
				vec_delta_x2 = vec_delta_x + vec_u2;

				auto vec_delta_x_norm = vec_delta_x.normalized();
				auto vec_delta_x1_norm = vec_delta_x1.normalized();
				auto vec_delta_x2_norm = vec_delta_x2.normalized();

				auto x1_proj = vec_delta_x_norm.dot(vec_delta_x1_norm);
				auto x2_proj = vec_delta_x_norm.dot(vec_delta_x2_norm);

				if (x1_proj > x2_proj)
					gamma = gamma1;
				else
					gamma = gamma2;
			};

			choose_correction();

			while (std::abs(gamma) > gamma_max && s2 / 2 > s2_min)
			{
				s2 /= 2;
				choose_correction();
			}
			if (std::abs(gamma) == std::numeric_limits<Scalar>::infinity())
			{
				s2 = (vec_delta_x + vec_uR).squaredNorm();
				choose_correction();
			}

			if (gamma == gamma1)
				vec_delta_x = vec_delta_x1;
			else
				vec_delta_x = vec_delta_x2;

//			gamma = -vec_delta_x0.dot(vec_uR) / vec_delta_x0.dot(vec_uF);
//			Eigen::VectorXd u = vec_uR + gamma * vec_uF;
//			u.array() *= one_minus_fixed_dof.array();
//			vec_delta_x += u;
//			Scalar const alpha = sqrt(s2 / vec_delta_x.squaredNorm());
//			vec_delta_x *= alpha;

			delta_lambda += gamma;
			lambda += gamma;

			mat_x = mat_x0 + as_matrix(vec_delta_x);
		}

		numbers << increment_num << "\t" << lambda << "\t" << mat_x(2, 1) << "\n" << std::flush;
//		if (max_norm < epsilon)
	}
}

void Matrix::assemble(
	Eigen::MatrixXd & mat_K,
	Eigen::VectorXd & vec_R,
	Eigen::VectorXd & vec_F,
	Eigen::VectorXd const & one_minus_fixed_dof) const
{
	vec_R.setZero();
	vec_F.setZero();
	mat_K.setZero();

	for (auto vtxh : boost::make_iterator_range(m_mesh.vertices()))
	{
		auto const vtx_idx = vtxh.idx();
		for (auto cellh : boost::make_iterator_range(m_mesh.vertex_cells(vtxh)))
		{
			auto const & cell_vtxhs = m_attrs.vtxhs[cellh];
			auto const & cell_vtx_idx = index_of<Eigen::Index>(cell_vtxhs, vtxh);
			auto const & cell_R = m_attrs.R[cellh];
			auto const & cell_F = m_attrs.F[cellh];

			// Residual force for node `a`.
			auto const & Ra = EigenConstTensorMap<4, 3>{cell_R.data()}.block<1, 3>(cell_vtx_idx, 0);
			vec_R.block<3, 1>(3 * vtx_idx, 0) += Ra;
			// External force for node `a`.
			auto const & Fa = EigenConstTensorMap<4, 3>{cell_F.data()}.block<1, 3>(cell_vtx_idx, 0);
			vec_F.block<3, 1>(3 * vtx_idx, 0) += Fa;
		}
	}
	// Compute sum to check equilibrium.
	//		Eigen::Vector3d sum; sum.setZero();
	//		Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor>>
	// mat_R_Nx3{mat_R.data(), mat_R.rows() / 3, 3}; 		for (Eigen::Index const node_idx
	// : boost::irange(mat_R_Nx3.rows())) 			sum += mat_R_Nx3.row(node_idx);

	auto const update_submatrix =
		[&mat_K, &attrs = m_attrs](auto const vtxh_src, auto const vtxh_dst, auto const cellh_)
	{
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

	// Zero-out penalised degrees of freedom.
	vec_R.array() *= one_minus_fixed_dof.array();
	vec_F.array() *= one_minus_fixed_dof.array();
	mat_K.array().colwise() *= one_minus_fixed_dof.array();
	mat_K.array().rowwise() *= one_minus_fixed_dof.transpose().array();
	Eigen::VectorXd const fixed_dof =
		Eigen::VectorXd::Constant(one_minus_fixed_dof.size(), 1.0) - one_minus_fixed_dof;
	mat_K += fixed_dof.asDiagonal();

	//	// Inspect matrix to calculate penalty of relative size.
	//	Scalar const penalty = pow(10, std::floor(std::log10(mat_K.lpNorm<Eigen::Infinity>()))) *
	// 1e1;
	//	// Set penalised degrees of freedom to penalty value.
	//	mat_K += penalty *
	//		(Eigen::VectorXd::Constant(one_minus_fixed_dof.size(), 1.0) - one_minus_fixed_dof)
	//			.asDiagonal();
}

std::tuple<Scalar, Scalar> Matrix::arc_length(
	Eigen::VectorXd const & vec_uF,
	Eigen::VectorXd const & vec_uR,
	Eigen::VectorXd const & vec_F,
	Eigen::VectorXd const & vec_delta_x,
	Scalar const delta_lambda,
	Scalar s2,
	Scalar const psi2)
{
	Scalar a, b, c, discriminant;
	a = vec_uF.squaredNorm() + psi2 * vec_F.squaredNorm();
	b = 2 * vec_uF.dot(vec_delta_x + vec_uR) + 2 * delta_lambda * psi2 * vec_F.squaredNorm();
	c = (vec_delta_x + vec_uR).squaredNorm() +
		psi2 * delta_lambda * delta_lambda * vec_F.squaredNorm() - s2;
	//	c = vec_uR.dot(2 * vec_delta_x + vec_uR) + vec_delta_x.squaredNorm() - s2;

	discriminant = b * b - 4 * a * c;

	if (discriminant < 0)
		return {std::numeric_limits<Scalar>::infinity(), std::numeric_limits<Scalar>::infinity()};

	//	while (discriminant < 0)
	//	{
	//		// TODO: decompose perpendicular component
	//		spdlog::debug(
	//			"uF = \n{}\nDx+uR = \n{}", vec_uF.normalized(), (vec_delta_x +
	// vec_uR).normalized()); 		c += s2; 		s2 *= 1.1; 		c -= s2; 		discriminant = b
	// *
	// b
	// -
	// 4
	// *
	// a
	// * c;
	//	}

	auto const gamma1 = (-b + sqrt(discriminant)) / (2 * a);
	auto const gamma2 = (-b - sqrt(discriminant)) / (2 * a);

	//	auto const error1 = (vec_delta_x + vec_uR + gamma1 * vec_uF).squaredNorm() +
	//		(delta_lambda + gamma1) * (delta_lambda + gamma1) * vec_F.squaredNorm();
	//	auto const error2 = (vec_delta_x + vec_uR + gamma2 * vec_uF).squaredNorm() +
	//		(delta_lambda + gamma2) * (delta_lambda + gamma2) * vec_F.squaredNorm();
	//	assert(abs(error1 - s2) < 0.000001);
	//	assert(abs(error2 - s2) < 0.000001);

	return {gamma1, gamma2};
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
