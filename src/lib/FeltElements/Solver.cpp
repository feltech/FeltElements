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
#include "MeshIO.hpp"
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
	auto const num_cells = static_cast<int>(m_mesh_io.mesh.n_cells());

	FE_PARALLEL(shared(num_cells, m_mesh_io.attrs, lambda))
	for (int cell_idx = 0; cell_idx < num_cells; cell_idx++)
	{
		Cellh cellh{cell_idx};
		auto const & cell_vtxhs = m_mesh_io.attrs.vtxhs[cellh];
		auto const & boundary_faces_vtxh_idxs = m_mesh_io.attrs.boundary_faces_vtxh_idxs[cellh];

		auto const & [K, F, R] = Derivatives::KR(
			m_mesh_io.attrs.x.for_element(cell_vtxhs),
			boundary_faces_vtxh_idxs,
			m_mesh_io.attrs.x.for_elements(cell_vtxhs, boundary_faces_vtxh_idxs),
			m_mesh_io.attrs.dN_by_dX[cellh],
			*m_mesh_io.attrs.material,
			*m_mesh_io.attrs.forces,
			lambda);
		m_mesh_attrs.F[cellh] = F;
		m_mesh_attrs.R[cellh] = R;
		m_mesh_attrs.K[cellh] = K;
	}
}

Scalar Base::find_approx_min_edge_length() const
{
	Scalar total_V = 0;
	Scalar min_V = std::numeric_limits<Scalar>::max();
	Scalar max_V = 0;
	Element::NodePositions min_X;

	for (auto const & X : MeshIO{m_mesh_io.mesh, m_mesh_io.attrs}.Xs())
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
		total_V / scalar(m_mesh_io.mesh.n_cells()),
		max_V,
		min_V,
		min_X);
	return edge_length;
}

namespace
{
void log_xs(MeshIO const & mesh_io)
{
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_DEBUG
	std::string xs_str;
	for (auto const & x : mesh_io.xs()) { xs_str += fmt::format("{}\n", x); }
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
	auto const consts = [&]() -> Constants
	{
		constexpr std::size_t step_target = 5;
		constexpr Scalar psi2 = 0;
		Scalar const residual_epsilon = find_approx_min_edge_length() * 1e-3;
		Scalar const s2_epsilon = std::pow(residual_epsilon, 2) * 1e-1;

		return {step_target, psi2, residual_epsilon, s2_epsilon};
	}();

	auto soln = [&]() -> Solution
	{
		Scalar const lambda = scalar(1e-1);
		Scalar const s2 = std::pow(consts.residual_epsilon, 2);
		return {Dofs{m_mesh_io}, as_matrix(m_mesh_io.attrs.x), lambda, s2};
	}();

	std::size_t step = consts.step_target;

	for (std::size_t increment_num = 0; increment_num < m_params.num_force_increments;
		 increment_num++)
	{
		++stats.force_increment_counter;
		++stats.step_counter;

		IncrementState const state = increment(increment_num, step, consts, soln);

		if (state == IncrementState::increment_is_needed)
		{
			continue;
		}

		if (state == IncrementState::done)
			break;

		correction(increment_num, step, consts, soln);
	}
}

Matrix::IncrementState Matrix::increment(
	std::size_t const increment_num, std::size_t & step, Constants const & consts, Solution & soln)
{
	soln.lambda = std::min(soln.lambda, scalar(1));

	update_elements_stiffness_and_residual(soln.lambda);
	assemble(soln.mat_K, soln.vec_R, soln.vec_F, soln.dofs.one_minus_fixed_dof);
	soln.vec_F /= soln.lambda;

	auto const mat_K_LU = soln.mat_K.partialPivLu();
	soln.vec_uR = mat_K_LU.solve(-soln.vec_R);

	if (soln.lambda != 1)
	{
		// Increase/decrease arc length if below/above step target.
		soln.s2 *= 2 * scalar(consts.step_target) / scalar(step + consts.step_target);
		// Initial optimistic assumption that next increment needs no arc length refinement.
		step = consts.step_target;
		// Initialise displacement to uncorrected solution.
		soln.vec_delta_x = soln.vec_uR;
		soln.mat_x += as_matrix(soln.vec_delta_x);
	}
	else
	{
		// If lambda is 1 then revert to standard non-arc-length solution.
		soln.mat_x += as_matrix(soln.vec_uR);
	}

	Scalar const residual_norm = soln.vec_uR.squaredNorm();

	SPDLOG_DEBUG(
		"increment = {}; total steps = {}; lambda = {}; s2 = {}; max_norm = {}; rcond = {};"
		" det(K) = {}",
		increment_num,
		stats.step_counter.load(),
		soln.lambda,
		soln.s2,
		residual_norm,
		mat_K_LU.rcond(),
		mat_K_LU.determinant());
	log_xs(m_mesh_io);

	stats.residual_norm = residual_norm;

	if (residual_norm > consts.residual_epsilon)
	{
		if (soln.s2 > consts.s2_epsilon)
			// Try next arc direction.
			return IncrementState::correction_is_needed;
	}
	else if (soln.lambda == 1)
	{
		// Finished.
		return IncrementState::done;
	}

	// Next load increment.

	soln.vec_uF = mat_K_LU.solve(soln.vec_F);
	Scalar const uF2 = soln.vec_uF.squaredNorm();
	// Estimate next lambda.
	Scalar const gamma = std::sqrt(soln.s2 / uF2);
	soln.lambda = std::min(soln.lambda + gamma, scalar(1));

	return IncrementState::increment_is_needed;
}

void Matrix::correction(
	std::size_t const increment_num, size_t & step, Constants const & consts, Solution & soln)
{
	for (step = 0; step < m_params.num_steps; ++step)
	{
		++stats.step_counter;

		update_elements_stiffness_and_residual(soln.lambda);
		assemble(soln.mat_K, soln.vec_R, soln.vec_F, soln.dofs.one_minus_fixed_dof);
		soln.vec_F /= soln.lambda;

		auto const mat_K_LU = soln.mat_K.partialPivLu();
		soln.vec_uR = mat_K_LU.solve(-soln.vec_R);
		soln.vec_uF = mat_K_LU.solve(soln.vec_F);

		Scalar const residual_norm = soln.vec_R.squaredNorm();

		stats.residual_norm = residual_norm;
		if (residual_norm < consts.residual_epsilon)
			return;

		Scalar const gamma = arc_length(
			soln.vec_uF,
			soln.vec_uR,
			soln.vec_F,
			soln.delta_lambda,
			consts.psi2,
			soln.vec_delta_x,
			soln.vec_delta_delta_x,
			soln.s2);

		SPDLOG_DEBUG(
			"increment = {}; step = {}; lambda = {}; delta_lambda = {}; gamma = {}; s2 = "
			"{}; norm = {}; rcond = {}; det(K) = {}; delta_x2 = {}",
			increment_num,
			step,
			soln.lambda,
			soln.delta_lambda,
			gamma,
			soln.s2,
			residual_norm,
			mat_K_LU.rcond(),
			mat_K_LU.determinant(),
			soln.vec_delta_x.squaredNorm());

		soln.delta_lambda += gamma;
		soln.lambda += gamma;
		soln.mat_x += as_matrix(soln.vec_delta_delta_x);
	}
}

void Matrix::assemble(
	MatrixX & mat_K, VectorX & vec_R, VectorX & vec_F, VectorX const & one_minus_fixed_dof) const
{
	vec_R.setZero();
	vec_F.setZero();
	mat_K.setZero();

	for (auto vtxh : m_mesh_io.vertices())
	{
		auto const vtx_idx = vtxh.idx();
		for (auto cellh : m_mesh_io.vertex_cells(vtxh))
		{
			auto const & cell_vtxhs = m_mesh_io.attrs.vtxhs[cellh];
			auto const & cell_vtx_idx = index_of<Index>(cell_vtxhs, vtxh);
			auto const & cell_R = m_mesh_io.attrs.R[cellh];
			auto const & cell_F = m_mesh_io.attrs.F[cellh];

			// Residual force for node `a`.
			auto const & Ra = EigenConstTensorMap<4, 3>{cell_R.data()}.block<1, 3>(cell_vtx_idx, 0);
			vec_R.block<3, 1>(3 * vtx_idx, 0) += Ra;
			// External force for node `a`.
			auto const & Fa = EigenConstTensorMap<4, 3>{cell_F.data()}.block<1, 3>(cell_vtx_idx, 0);
			vec_F.block<3, 1>(3 * vtx_idx, 0) += Fa;
		}
	}

	auto const update_submatrix = [&mat_K, &attrs = m_mesh_io.attrs](
									  auto const vtxh_src, auto const vtxh_dst, auto const cellh_)
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

	for (auto vtxh : m_mesh_io.vertices())
	{
		for (auto cellh : m_mesh_io.vertex_cells(vtxh)) update_submatrix(vtxh, vtxh, cellh);
	}
	for (auto heh : m_mesh_io.halfedges())
	{
		auto const & halfedge = m_mesh_io.mesh.halfedge(heh);
		auto const & vtxh_src = halfedge.from_vertex();
		auto const & vtxh_dst = halfedge.to_vertex();
		for (auto cellh : m_mesh_io.halfedge_cells(heh))
			update_submatrix(vtxh_src, vtxh_dst, cellh);
	}

	// Zero-out penalised degrees of freedom.
	vec_R.array() *= one_minus_fixed_dof.array();
	vec_F.array() *= one_minus_fixed_dof.array();
	VectorX const fixed_dof =
		VectorX::Constant(one_minus_fixed_dof.size(), 1.0) - one_minus_fixed_dof;

	mat_K.array().colwise() *= one_minus_fixed_dof.array();
	mat_K.array().rowwise() *= one_minus_fixed_dof.transpose().array();
	// Set non-zero fixed DOFs on diagonal to absolute row sums, to hopefully ensure K is diagonally
	// dominant and thus non-singular.
	mat_K.diagonal() += fixed_dof.transpose() +
		(fixed_dof.asDiagonal() * mat_K).array().abs().colwise().sum().matrix();

	//	// Inspect matrix to calculate penalty of relative size.
	//	Scalar const penalty = pow(10, std::floor(std::log10(mat_K.lpNorm<Eigen::Infinity>()))) *
	// 1e1;
	//	// Set penalised degrees of freedom to penalty value.
	//	mat_K += penalty * fixed_dof.asDiagonal();
}

Scalar Matrix::arc_length(
	VectorX const & vec_uF,
	VectorX const & vec_uR,
	VectorX const & vec_F,
	Scalar delta_lambda,
	Scalar psi2,
	VectorX & vec_delta_x,
	VectorX & vec_delta_delta_x,
	Scalar & s2)
{
	// Clamp arc length to minimum allowed radius.
	VectorX const delta_xu = vec_delta_x + vec_uR;
	VectorX const n = vec_uF.normalized();
	VectorX const vec_s_min = delta_xu - n * n.dot(delta_xu);
	Scalar const s2_min = vec_s_min.squaredNorm() * scalar(2);
	if (s2 < s2_min)
		s2 = s2_min;

	// Solve quadratic, giving two potential solutions.

	auto [gamma1, gamma2] =
		arc_length_multipliers(vec_uF, vec_uR, vec_F, vec_delta_x, delta_lambda, s2, psi2);

	// Choose arc direction solution most aligned to current displacement direction.

	auto const vec_u1 = vec_uR + gamma1 * vec_uF;
	auto const vec_u2 = vec_uR + gamma2 * vec_uF;
	auto const vec_delta_x1 = vec_delta_x + vec_u1;
	auto const vec_delta_x2 = vec_delta_x + vec_u2;
	auto const x1_proj = vec_delta_x.dot(vec_delta_x1);
	auto const x2_proj = vec_delta_x.dot(vec_delta_x2);

	Scalar gamma;
	if (x1_proj > x2_proj)
	{
		gamma = gamma1;
		vec_delta_delta_x = vec_u1;
		vec_delta_x = vec_delta_x1;
	}
	else
	{
		gamma = gamma2;
		vec_delta_delta_x = vec_u2;
		vec_delta_x = vec_delta_x2;
	}
	return gamma;
}

std::tuple<Scalar, Scalar> Matrix::arc_length_multipliers(
	VectorX const & vec_uF,
	VectorX const & vec_uR,
	VectorX const & vec_F,
	VectorX const & vec_delta_x,
	Scalar const delta_lambda,
	Scalar const s2,
	Scalar const psi2)
{
	Scalar const uF2 = vec_uF.squaredNorm();
	Scalar const F2 = vec_F.squaredNorm();
	Scalar a, b, c, discriminant;

	a = uF2 + psi2 * vec_F.squaredNorm();
	b = 2 * vec_uF.dot(vec_delta_x + vec_uR) + 2 * delta_lambda * psi2 * F2;
	c = (vec_delta_x + vec_uR).squaredNorm() + psi2 * delta_lambda * delta_lambda * F2 - s2;
	//	c = vec_uR.dot(2 * vec_delta_x + vec_uR) + vec_delta_x.squaredNorm() - s2;

	discriminant = b * b - 4 * a * c;

	if (discriminant < 0)
		return {std::numeric_limits<Scalar>::infinity(), std::numeric_limits<Scalar>::infinity()};

	Scalar const sqrt_discriminant = std::sqrt(discriminant);

	auto const gamma1 = (-b + sqrt_discriminant) / (2 * a);
	auto const gamma2 = (-b - sqrt_discriminant) / (2 * a);

	return {gamma1, gamma2};
}

Matrix::Index Matrix::Dofs::calc_num_dofs(const MeshIO & mesh_io)
{
	auto const num_vertices = static_cast<Index>(mesh_io.mesh.n_vertices());
	return static_cast<Index>(Node::dim) * num_vertices;
}

Matrix::VectorX Matrix::Dofs::calc_fixed_dof(const MeshIO & mesh_io)
{
	auto const rows = static_cast<Index>(mesh_io.mesh.n_vertices());
	constexpr auto cols = static_cast<Index>(Node::dim);

	Matrix::EigenMapTensorVerticesConst const & map{
		mesh_io.attrs.fixed_dof[Vtxh{0}].data(), rows, cols};
	// Copy to remove per-row alignment.
	Matrix::VectorX vec =
		Eigen::Map<Matrix::VectorX>{Matrix::VerticesMatrix{map}.data(), rows * cols};
	return vec;
}

Matrix::Dofs::Dofs(const MeshIO & mesh_io)
	: num_dofs{calc_num_dofs(mesh_io)},
	  fixed_dof{calc_fixed_dof(mesh_io)},
	  one_minus_fixed_dof{VectorX::Ones(fixed_dof.size()) - fixed_dof}
{
}

void Gauss::solve()
{
	Scalar const epsilon = find_approx_min_edge_length() * 1e-5;
	using Tensor::Func::all;

	std::vector<Node::Pos> u(m_mesh_io.mesh.n_vertices(), {0, 0, 0});
	Scalar max_norm = 0;

	Node::Force const force_increment =
		m_mesh_io.attrs.forces->F_by_m / m_params.num_force_increments;
	m_mesh_attrs.forces->F_by_m = 0;

	for (std::size_t increment_num = 0; increment_num < m_params.num_force_increments;
		 increment_num++)
	{
		m_mesh_attrs.forces->F_by_m += force_increment;
		stats.force_increment_counter++;

		for (std::size_t step = 0; step < m_params.num_steps; step++)
		{
			stats.step_counter++;
			pauser.wait_while_paused();

			SPDLOG_DEBUG("Gauss iteration {}:{}", increment_num, step);

			update_elements_stiffness_and_residual();

			max_norm = 0;
			for (auto & u_a : u) u_a *= 0.99;  //.zeros();

			FE_PARALLEL(shared(m_mesh_io, m_mesh_attrs, u, force_increment) reduction(max
																					  : max_norm))
			for (Tensor::Index a = 0; a < m_mesh_io.mesh.n_vertices(); a++)
			{
				Vtxh vtxh_src{static_cast<int>(a)};
				auto const & fixed_dof = m_mesh_io.attrs.fixed_dof[vtxh_src];
				if (Tensor::Func::all_of(fixed_dof != 0))
					continue;

				Node::Force Ra = 0;
				Tensor::Matrix<3> Kaa = 0;
				Node::Force Ka_u = 0;
				for (auto cellh : m_mesh_io.vertex_cells(vtxh_src))
				{
					auto const & cell_vtxhs = m_mesh_io.attrs.vtxhs[cellh];
					auto const & cell_a = index_of(cell_vtxhs, vtxh_src);
					auto const & cell_R = m_mesh_io.attrs.R[cellh];
					auto const & cell_K = m_mesh_io.attrs.K[cellh];

					Ra += cell_R(cell_a, all);
					Kaa += cell_K(cell_a, all, cell_a, all);
				}
				//			diag(Kaa) += penalty * fixed_dof;

				for (auto heh : m_mesh_io.outgoing_halfedges(vtxh_src))
				{
					Tensor::Multi<Node::dim, Node::dim> Kab = 0;
					auto const & halfedge = m_mesh_io.mesh.halfedge(heh);
					auto const & vtxh_dst = halfedge.to_vertex();
					auto const b = static_cast<Tensor::Index>(vtxh_dst.idx());

					for (auto cellh : m_mesh_io.halfedge_cells(heh))
					{
						auto const & cell_K = m_mesh_io.attrs.K[cellh];
						auto const & cell_vtxhs = m_mesh_io.attrs.vtxhs[cellh];
						auto const cell_a = index_of(cell_vtxhs, vtxh_src);
						auto const cell_b = index_of(cell_vtxhs, vtxh_dst);
						Kab += cell_K(cell_a, all, cell_b, all);
					}
					Ka_u += Kab % u[b];
				}

				using Tensor::Func::inv;
				u[a] = inv(Kaa) % (-Ra - Ka_u) * (1.0 - fixed_dof);
				m_mesh_attrs.x[vtxh_src] += u[a];
				using Tensor::Func::abs;
				using Tensor::Func::max;
				// Note: double-reduction in case OpenMP disabled.
				max_norm = std::max(max_norm, max(abs(u[a])));
			}

			SPDLOG_DEBUG("Max norm: {}", max_norm);
			stats.residual_norm = max_norm;

			log_xs(m_mesh_io);

			if (max_norm < epsilon)
				break;
		}
	}
}
}  // namespace FeltElements::Solver
