#pragma once
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

namespace FeltElements::internal
{
template <typename Matrix>
using is_matrix = typename std::enable_if<
	(Matrix::RowsAtCompileTime > 1 && Matrix::ColsAtCompileTime > 1), void*>::type;

template <typename Matrix>
using is_vector = typename std::enable_if<
	(Matrix::RowsAtCompileTime == 1 || Matrix::ColsAtCompileTime == 1), void*>::type;

template <typename Matrix, is_matrix<Matrix> = nullptr>
auto to_tensor(Matrix const & matrix)
{
	static_assert(
		!(Matrix::Options & Eigen::RowMajor),
		"Cannot map a row-major matrix to a tensor since row-major tensors are not yet"
		" supported by Eigen");

	return Eigen::TensorMap<
		   Eigen::Tensor<FeltElements::Tetrahedron::Scalar const, 2,
		Matrix::Options & Eigen::RowMajor ? Eigen::RowMajor : Eigen::ColMajor>>{
			matrix.data(), {Matrix::RowsAtCompileTime, Matrix::ColsAtCompileTime}
		};
}

template <typename Matrix, is_vector<Matrix> = nullptr>
auto to_tensor(Matrix const & matrix)
{
	return Eigen::TensorMap<Eigen::Tensor<FeltElements::Tetrahedron::Scalar const, 1>>{
		matrix.data(), {std::max(Matrix::RowsAtCompileTime, Matrix::ColsAtCompileTime)}
	};
}

}
