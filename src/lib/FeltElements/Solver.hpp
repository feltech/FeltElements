#pragma once

#include <Fastor/Fastor.h>
#include <eigen3/Eigen/Core>
#include <boost/preprocessor/cat.hpp>

#include "Typedefs.hpp"

namespace FeltElements
{
class Attributes;

namespace Solver
{
void update_elements_stiffness_and_residual(Mesh const & mesh, Attributes & attributes);

namespace LDLT
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

std::size_t solve(Mesh & mesh, Attributes & attrib, std::size_t max_steps);
}  // namespace LDLT

namespace Gauss
{
std::size_t solve(Mesh & mesh, Attributes & attrib, std::size_t max_steps);
}

}  // namespace Solver
}  // namespace FeltElements