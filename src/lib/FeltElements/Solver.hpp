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
	boost::atomic_flag pause{};
	std::mutex running{};
};

void update_elements_stiffness_and_residual(Mesh const & mesh, Attributes & attributes);

namespace Matrix
{
#define EIGEN_FASTOR_ALIGN BOOST_PP_CAT(Eigen::Aligned, FASTOR_MEMORY_ALIGNMENT_VALUE)
template <int rows, int cols>
using EigenConstTensorMap =
	Eigen::Map<Eigen::Matrix<Scalar, rows, cols, Eigen::RowMajor> const, EIGEN_FASTOR_ALIGN>;
using VerticesMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor>;
using EigenMapOvmVertices = Eigen::Map<VerticesMatrix const>;
using EigenMapTensorVertices = Eigen::Map<
	VerticesMatrix const,
	Eigen::Unaligned,
	Eigen::Stride<(FASTOR_MEMORY_ALIGNMENT_VALUE / sizeof(Scalar)), 1>>;
using EigenFixedDOFs = Eigen::VectorXd;

Scalar solve(Mesh & mesh, Attributes & attrib, Params params, Stats * stats = nullptr);
}  // namespace Matrix

namespace Gauss
{
Scalar solve(Mesh & mesh, Attributes & attrib, Params params, Stats * stats = nullptr);
}

}  // namespace Solver
}  // namespace FeltElements