#include "Solver.hpp"

#include <cmath>
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

Scalar Base::find_approx_min_edge_length() const
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
	// Edge length assuming a regular tetrahedron, divided by a constant.
	Scalar const edge_length = scalar(std::pow(6 * std::sqrt(Scalar{2.0}) * min_V, 1.0 / 3.0));
	SPDLOG_DEBUG(
		"edge_length = {}, since total V = {}; mean V = {}; max V = {}; min V = {}; min X = \n{}",
		edge_length,
		total_V,
		total_V / scalar(m_mesh.n_cells()),
		max_V,
		min_V,
		min_X);
	return edge_length;
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

template <typename T>
int sgn(T val)
{
	return (T(0) < val) - (val < T(0));
}
}  // namespace

void Matrix::solve()
{
	Scalar residual_epsilon = find_approx_min_edge_length();

	VectorX const & mat_fixed_dof = ([&attrs = m_attrs,
											  rows = static_cast<Eigen::Index>(m_mesh.n_vertices()),
											  cols = static_cast<Eigen::Index>(Node::dim)]() {
		EigenMapTensorVerticesConst const & map{attrs.fixed_dof[Vtxh{0}].data(), rows, cols};
		// Copy to remove per-row alignment.
		VectorX vec = Eigen::Map<VectorX>{VerticesMatrix{map}.data(), rows * cols};
		return vec;
	})();
	VectorX const one_minus_fixed_dof = VectorX::Ones(mat_fixed_dof.size()) - mat_fixed_dof;

	auto const num_vertices = static_cast<Eigen::Index>(m_mesh.n_vertices());
	auto const num_dofs = static_cast<Eigen::Index>(Node::dim) * num_vertices;

	VectorX vec_F{num_dofs};
	VectorX vec_R{num_dofs};
	MatrixX mat_K{num_dofs, num_dofs};
	VectorX vec_u{num_dofs};
	VectorX vec_uF{num_dofs};
	VectorX vec_uR{num_dofs};
	VectorX vec_delta_x{num_dofs};
	Scalar delta_lambda = 0;
	Scalar lambda = scalar(1e-1);
	Scalar gamma = 0;
	constexpr std::size_t step_target = 5;
	std::size_t step = 0;
	Scalar s2 = scalar(1);
	Scalar const psi2 = 0;
	residual_epsilon = scalar(1e-6);
	constexpr Scalar s2_epsilon = 1e-8;
	auto mat_x = as_matrix(m_attrs.x);

	VerticesMatrix mat_X = as_matrix(m_attrs.X);

	for (std::size_t increment_num = 0; increment_num < m_params.num_force_increments;
		 increment_num++)
	{
		stats.force_increment_counter++;

		enum class State
		{
			arc,
			increment,
			done
		};

		State const state = [&]
		{
			lambda = std::min(lambda, scalar(1));

			update_elements_stiffness_and_residual(lambda);
			assemble(mat_K, vec_R, vec_F, one_minus_fixed_dof);
			vec_F /= lambda;

			s2 *= 2 * scalar(step_target) / scalar(step + step_target);
			step = 0;

			auto const mat_K_LU = mat_K.partialPivLu();
			vec_uR = mat_K_LU.solve(-vec_R);
			vec_uF = mat_K_LU.solve(vec_F);
			Scalar const uF2 = vec_uF.squaredNorm();

			gamma = std::sqrt(s2 / uF2);
			while (lambda + gamma > 1)
			{
				s2 *= 1e-1;
				gamma = std::sqrt(s2 / uF2);
			}
			gamma = 0;

			delta_lambda = gamma;
			lambda += gamma;
			vec_delta_x = vec_uR + gamma * vec_uF;
			vec_delta_x.array() *= one_minus_fixed_dof.array();
			mat_x += as_matrix(vec_delta_x);

			Scalar const rcond = mat_K_LU.rcond();
			Scalar const residual_norm = vec_R.squaredNorm();

			SPDLOG_DEBUG(
				"increment = {}; total steps = {}; lambda = {}; s2 = {}; max_norm = "
				"{}; rcond = {}; det(K) = {}",
				increment_num,
				stats.step_counter.load(),
				lambda,
				s2,
				residual_norm,
				mat_K_LU.rcond(),
				mat_K_LU.determinant());
			log_xs(m_mesh, m_attrs);

			stats.residual_norm = residual_norm;

			if (residual_norm > residual_epsilon)
			{
				if (s2 > s2_epsilon)
					return State::arc;
			}
			else if (lambda == 1)
				return State::done;

			gamma = std::sqrt(s2 / uF2);
			lambda = std::min(lambda + gamma, scalar(1));
			return State::increment;
		}();

		if (state == State::increment)
			continue;

		if (state == State::done)
			break;

		[&]
		{
			for (step = 0; step < m_params.num_steps; ++step)
			{
				++stats.step_counter;

				update_elements_stiffness_and_residual(lambda);
				assemble(mat_K, vec_R, vec_F, one_minus_fixed_dof);
				vec_F /= lambda;

				auto const mat_K_LU = mat_K.partialPivLu();
				vec_uR = mat_K_LU.solve(-vec_R);
				vec_uF = mat_K_LU.solve(vec_F);

				Scalar const residual_norm = vec_R.squaredNorm();

				SPDLOG_DEBUG(
					"increment = {}; step = {}; lambda = {}; delta_lambda = {}; gamma = {}; s2 = "
					"{}; norm = {}; rcond = {}; det(K) = {}; delta_x2 = {}",
					increment_num,
					step,
					lambda,
					delta_lambda,
					gamma,
					s2,
					residual_norm,
					mat_K_LU.rcond(),
					mat_K_LU.determinant(),
					vec_delta_x.squaredNorm());

				stats.residual_norm = residual_norm;
				if (residual_norm < residual_epsilon)
					return;

				{
					VectorX const delta_xu = vec_delta_x + vec_uR;
					VectorX const n = vec_uF.normalized();
					VectorX const vec_s_min = delta_xu - n * n.dot(delta_xu);
					Scalar const s2_min = vec_s_min.squaredNorm() * scalar(1.0 + 1e-5);
					if (s2 < s2_min)
						s2 = s2_min;
				}

				{
					auto [gamma1, gamma2] =
						arc_length(vec_uF, vec_uR, vec_F, vec_delta_x, delta_lambda, s2, psi2);

					auto const vec_u1 = vec_uR + gamma1 * vec_uF;
					auto const vec_u2 = vec_uR + gamma2 * vec_uF;
					auto const vec_delta_x1 = vec_delta_x + vec_u1;
					auto const vec_delta_x2 = vec_delta_x + vec_u2;
					auto const x1_proj = vec_delta_x.dot(vec_delta_x1);
					auto const x2_proj = vec_delta_x.dot(vec_delta_x2);

					if (x1_proj > x2_proj)
					{
						gamma = gamma1;
						vec_u = vec_u1;
						vec_delta_x = vec_delta_x1;
					}
					else
					{
						gamma = gamma2;
						vec_u = vec_u2;
						vec_delta_x = vec_delta_x2;
					}
				};

				delta_lambda += gamma;
				lambda += gamma;
				mat_x += as_matrix(vec_u);
			}
		}();
	}
}

void Matrix::assemble(
	MatrixX & mat_K, VectorX & vec_R, VectorX & vec_F, VectorX const & one_minus_fixed_dof) const
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
	VectorX const fixed_dof =
		VectorX::Constant(one_minus_fixed_dof.size(), 1.0) - one_minus_fixed_dof;
	mat_K += fixed_dof.asDiagonal();

	//	// Inspect matrix to calculate penalty of relative size.
	//	Scalar const penalty = pow(10, std::floor(std::log10(mat_K.lpNorm<Eigen::Infinity>()))) *
	// 1e1;
	//	// Set penalised degrees of freedom to penalty value.
	//	mat_K += penalty *
	//		(VectorX::Constant(one_minus_fixed_dof.size(), 1.0) - one_minus_fixed_dof)
	//			.asDiagonal();
}

std::tuple<Scalar, Scalar> Matrix::arc_length(
	VectorX const & vec_uF,
	VectorX const & vec_uR,
	VectorX const & vec_F,
	VectorX const & vec_delta_x,
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

	Scalar const sqrt_discriminant = std::sqrt(discriminant);

	auto const gamma1 = (-b + sqrt_discriminant) / (2 * a);
	auto const gamma2 = (-b - sqrt_discriminant) / (2 * a);

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
	Scalar const epsilon = find_approx_min_edge_length() * 1e-3;
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
			stats.residual_norm = max_norm;

			log_xs(m_mesh, m_attrs);

			if (max_norm < epsilon)
				break;
		}
	}
}
}  // namespace FeltElements::Solver
