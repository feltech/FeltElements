#pragma once
#include <atomic>
#include <mutex>

#include <Fastor/Fastor.h>
#include <boost/atomic.hpp>
#include <boost/atomic/atomic_flag.hpp>
#include <boost/preprocessor/cat.hpp>
#include <eigen3/Eigen/Core>

#include "MeshIO.hpp"
#include "Typedefs.hpp"

namespace FeltElements
{
struct Attributes;

namespace Solver
{
struct Params
{
	std::size_t num_steps{1};
	std::size_t num_force_increments{1};
};

struct Stats
{
	std::atomic_uint step_counter{0};
	std::atomic_uint force_increment_counter{0};
	std::atomic<Scalar> residual_norm{0};
};

class Base
{
public:
	Base(Mesh const & mesh, Attributes & attrs, Params params) noexcept
		: m_mesh_io{mesh, attrs}, m_mesh_attrs{attrs}, m_params{params}
	{
	}
	virtual ~Base() = default;
	virtual void solve() = 0;

protected:
	class Pauser
	{
	public:
		struct Unpauser
		{
			Unpauser(Unpauser const &) = delete;
			Unpauser(Unpauser &&) = delete;
			void operator=(Unpauser const &) = delete;
			void operator=(Unpauser &&) = delete;
			~Unpauser()
			{
				pause.flag.clear();
				pause.flag.notify_one();
				pause.lock.unlock();
			}
			Pauser & pause;
		};

		Unpauser scoped_pause() noexcept
		{
			flag.test_and_set();
			lock.lock();
			return Unpauser{*this};
		}

		void wait_while_paused() noexcept
		{
			if (!flag.test())
				return;
			lock.unlock();
			flag.wait(true);
			lock.lock();
		}

	private:
		std::mutex lock{};
		boost::atomic_flag flag{};
	};

	void update_elements_stiffness_and_residual(Scalar lambda = 1.0);
	[[nodiscard]] Scalar find_approx_min_edge_length() const;

public:
	using Unpauser = Pauser::Unpauser;
	Pauser pauser{};
	Stats stats{};

protected:
	MeshIO const m_mesh_io;
	Attributes & m_mesh_attrs;

	Params m_params;
};

class Matrix : public Base
{
#define EIGEN_FASTOR_ALIGN BOOST_PP_CAT(Eigen::Aligned, FASTOR_MEMORY_ALIGNMENT_VALUE)
	template <int rows, int cols>
	using EigenConstTensorMap =
		Eigen::Map<Eigen::Matrix<Scalar, rows, cols, Eigen::RowMajor> const, EIGEN_FASTOR_ALIGN>;
	using OvmVtxMatrix = Eigen::Matrix<OvmScalar, Eigen::Dynamic, 3, Eigen::RowMajor>;

public:
	using Index = Eigen::Index;
	using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
	using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
	using EigenMapOvmVertices = Eigen::Map<OvmVtxMatrix const>;
	using VerticesMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor>;
	using EigenMapTensorVerticesConst = Eigen::Map<
		VerticesMatrix const,
		Eigen::Unaligned,
		Eigen::Stride<(FASTOR_MEMORY_ALIGNMENT_VALUE / sizeof(Scalar)), 1>>;
	using EigenMapTensorVertices = Eigen::Map<
		VerticesMatrix,
		Eigen::Unaligned,
		Eigen::Stride<(FASTOR_MEMORY_ALIGNMENT_VALUE / sizeof(Scalar)), 1>>;

	using Base::Base;
	void solve() final;

private:
	enum class IncrementState
	{
		correction_is_needed,
		increment_is_needed,
		done
	};

	struct Constants
	{
		std::size_t const step_target;
		Scalar const psi2;
		Scalar const residual_epsilon;
		Scalar const s2_epsilon;
	};

	struct Dofs
	{
		explicit Dofs(MeshIO const & mesh_io);

		Index const num_dofs;
		VectorX const fixed_dof;
		VectorX const one_minus_fixed_dof;

	private:
		static Index calc_num_dofs(const MeshIO & mesh_io);
		static VectorX calc_fixed_dof(const MeshIO & mesh_io);
	};

	struct Solution
	{
		// Solution scope.
		Dofs dofs;
		EigenMapTensorVertices mat_x;
		Scalar lambda;
		Scalar s2;
		// Step scope.
		VectorX vec_F{dofs.num_dofs};
		VectorX vec_R{dofs.num_dofs};
		MatrixX mat_K{dofs.num_dofs, dofs.num_dofs};
		VectorX vec_delta_x{dofs.num_dofs};
		Scalar delta_lambda{0};
		// Transients.
		VectorX vec_uR{dofs.num_dofs};
		VectorX vec_uF{dofs.num_dofs};
		VectorX vec_delta_delta_x{dofs.num_dofs};
	};

	IncrementState increment(
		std::size_t increment_num,
		std::size_t & step,
		Constants const & consts,
		Solution& soln);

	void correction(
		std::size_t increment_num,
		std::size_t & step,
		Constants const & consts,
		Solution& soln);

	void assemble(
		MatrixX & mat_K,
		VectorX & vec_R,
		VectorX & vec_F,
		VectorX const & one_minus_fixed_dof) const;

	static Scalar arc_length(
		VectorX const & vec_uF,
		VectorX const & vec_uR,
		VectorX const & vec_F,
		Scalar delta_lambda,
		Scalar psi2,
		VectorX & vec_delta_x,
		VectorX & vec_delta_delta_x,
		Scalar & s2);

	static std::tuple<Scalar, Scalar> arc_length_multipliers(
		VectorX const & vec_uF,
		VectorX const & vec_uR,
		VectorX const & vec_F,
		VectorX const & vec_delta_x,
		Scalar delta_lambda,
		Scalar s2,
		Scalar psi2);

	template <typename T>
	static std::enable_if_t<
		std::is_member_function_pointer_v<decltype(&T::data_vector)>,
		EigenMapTensorVertices>
	as_matrix(T & prop)
	{
		auto & vec = prop.data_vector();
		return {vec[0].data(), static_cast<Index>(vec.size()), Node::dim};
	}

	static Eigen::Map<VerticesMatrix const> as_matrix(VectorX const & vec)
	{
		return {vec.data(), vec.size() / static_cast<Index>(Node::dim), Node::dim};
	}

};	// namespace Matrix

class Gauss : public Base
{
public:
	using Base::Base;
	void solve() final;
};

}  // namespace Solver
}  // namespace FeltElements