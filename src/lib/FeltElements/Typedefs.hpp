#pragma once
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <unsupported/Eigen/CXX11/Tensor>

namespace FeltElements
{
using Mesh = OpenVolumeMesh::GeometricTetrahedralMeshV3d;
using Scalar = Mesh::PointT::value_type;

namespace Tensor
{
template <Eigen::Index rows = 3, Eigen::Index cols = rows>
using Matrix = Eigen::TensorFixedSize<Scalar, Eigen::Sizes<rows, cols>>;
template <Eigen::Index dim = 3>
using Vector = Eigen::TensorFixedSize<Scalar, Eigen::Sizes<dim>>;

using Index = Eigen::Index;
using IndexPair = Eigen::IndexPair<Index>;
template <Eigen::Index num_pairs>
using IndexPairs = Eigen::array<IndexPair, num_pairs>;
}

namespace Node
{
static constexpr Tensor::Index count = 4;
static constexpr Tensor::Index dim = 3;

using Pos = Tensor::Vector<dim>;
using PosMap = Eigen::TensorMap<Pos>;
using Positions = Eigen::TensorFixedSize<Scalar, Eigen::Sizes<count, dim>>;
using Force = Tensor::Vector<dim>;
// Note: col-major map to matrix for solver requires unintuitive layout.
using Forces = Eigen::TensorFixedSize<Scalar, Eigen::Sizes<dim, count>>;
}  // namespace Node

namespace Element
{
using IsoCoordDerivative = Tensor::Matrix<Node::dim, Node::count>;
using ShapeDerivative = Tensor::Matrix<Node::count, Node::dim>;
using Elasticity =
	Eigen::TensorFixedSize<Scalar, Eigen::Sizes<Node::dim, Node::dim, Node::dim, Node::dim>>;
using ShapeCartesianTransform = Tensor::Matrix<4, 4>;
using CartesianDerivative = Tensor::Matrix<Node::dim, Node::count>;
using Gradient = Tensor::Matrix<Node::dim, Node::dim>;
using Stress = Gradient;
// Note: col-major map to matrix for solver requires unintuitive layout.
using Stiffness =
	Eigen::TensorFixedSize<Scalar, Eigen::Sizes<Node::dim, Node::count, Node::dim, Node::count>>;
using StiffnessAndForces = std::tuple<Stiffness, Node::Forces>;
}  // namespace Element
}  // namespace FeltElements
