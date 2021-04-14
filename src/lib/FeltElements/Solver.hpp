#pragma once
#include <atomic>
#include <mutex>

#include <Fastor/Fastor.h>
#include <boost/atomic.hpp>
#include <boost/atomic/atomic_flag.hpp>
#include <boost/preprocessor/cat.hpp>
#include <eigen3/Eigen/Core>

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
	std::atomic<Scalar> max_norm{0};
};

class Base
{
public:
	Base(Mesh & mesh, Attributes & attrs, Params params) noexcept
		: m_mesh{mesh}, m_attrs{attrs}, m_params{std::move(params)}
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
				pause.flag.clear(boost::memory_order_release);
				pause.flag.notify_one();
				pause.lock.unlock();
			}
			Pauser & pause;
		};

		Unpauser scoped_pause() noexcept
		{
			flag.test_and_set(boost::memory_order_acquire);
			lock.lock();
			return Unpauser{*this};
		}

		void wait_while_paused() noexcept
		{
			if (!flag.test(boost::memory_order_relaxed))
				return;
			lock.unlock();
			flag.wait(true, boost::memory_order_relaxed);
			lock.lock();
		}

	private:
		std::mutex lock{};
		boost::atomic_flag flag{};
		std::scoped_lock<std::mutex> m_running{lock};
	};

	void update_elements_stiffness_and_residual(Scalar lambda = 1.0);
	Scalar find_epsilon() const;

public:
	using Unpauser = Pauser::Unpauser;
	Pauser pauser{};
	Stats stats{};

protected:
	Mesh & m_mesh;
	Attributes & m_attrs;
	Params m_params;
};

class Matrix : public Base
{
#define EIGEN_FASTOR_ALIGN BOOST_PP_CAT(Eigen::Aligned, FASTOR_MEMORY_ALIGNMENT_VALUE)
	template <int rows, int cols>
	using EigenConstTensorMap =
		Eigen::Map<Eigen::Matrix<Scalar, rows, cols, Eigen::RowMajor> const, EIGEN_FASTOR_ALIGN>;
	using VerticesMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor>;

public:
	using EigenMapOvmVertices = Eigen::Map<VerticesMatrix const>;
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
	void assemble(
		Eigen::MatrixXd & mat_K,
		Eigen::VectorXd & vec_R,
		Eigen::VectorXd & vec_F,
		Eigen::VectorXd const & one_minus_fixed_dof) const;

	std::tuple<Scalar, Scalar> arc_length(
		Eigen::VectorXd const & vec_uF,
		Eigen::VectorXd const & vec_uR,
		Eigen::VectorXd const & vec_delta_x,
		Eigen::VectorXd const & vec_F,
		Scalar const delta_lambda,
		Scalar const s2,
		Scalar const psi2);

	template <typename T>
	static std::enable_if_t<
		std::is_member_function_pointer_v<decltype(&T::data_vector)>,
		EigenMapTensorVertices>
	as_matrix(T & prop)
	{
		auto & vec = prop.data_vector();
		return {vec[0].data(), static_cast<Eigen::Index>(vec.size()), Node::dim};
	}

	static Eigen::Map<VerticesMatrix const> as_matrix(Eigen::VectorXd const & vec)
	{
		return {vec.data(), vec.size() / static_cast<Eigen::Index>(Node::dim), Node::dim};
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