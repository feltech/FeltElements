#pragma once
#include <Fastor/Fastor.h>

#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>

namespace FeltElements
{
using Mesh = OpenVolumeMesh::GeometricTetrahedralMeshV3d;
using Vtxh = OpenVolumeMesh::VertexHandle;
using Vtx = Mesh::PointT;
template <std::size_t dim>
using Vtxhs = std::array<Vtxh, dim>;
using Cellh = OpenVolumeMesh::CellHandle;
using Scalar = Mesh::PointT::value_type;

namespace Tensor
{
enum
{
	i,
	j,
	k,
	l,
	m,
	n,
	a,
	b
};
namespace Func = Fastor;
using Index = Fastor::FASTOR_INDEX;
template <Index... indices>
using Indices = Fastor::Index<indices...>;
template <Index... indices>
using Idxs = Indices<indices...>;
template <Index... indices>
using Order = Fastor::OIndex<indices...>;

template <typename value_type, Index... indices>
using Base = Fastor::Tensor<value_type, indices...>;

template <Index rows = 3, Index cols = rows>
using Matrix = Base<Scalar, rows, cols>;
template <Index dim = 3>
using Vector = Base<Scalar, dim>;
template <Index... indices>
using Multi = Base<Scalar, indices...>;

template <typename value_type, Index... indices>
using BaseMap = Fastor::TensorMap<value_type, indices...>;
template <Index... indices>
using Map = BaseMap<Scalar, indices...>;
template <Index... indices>
using ConstMap = BaseMap<Scalar const, indices...>;
}  // namespace Tensor

namespace Node
{
static constexpr Tensor::Index dim = 3;

using Pos = Tensor::Vector<dim>;
using Force = Tensor::Vector<dim>;
}  // namespace Node

namespace Element
{
static constexpr Tensor::Index count = 4;

using Vtxhs = FeltElements::Vtxhs<count>;
using IsoCoordDerivative = Tensor::Matrix<Node::dim, count>;
using ShapeDerivative = Tensor::Matrix<count, Node::dim>;
using SurfaceShapeDerivative = Tensor::Matrix<count - 1, Node::dim - 1>;
using ShapeDerivativeDeterminant = Tensor::Multi<count, count, count>;
using Elasticity = Tensor::Multi<Node::dim, Node::dim, Node::dim, Node::dim>;
using ShapeCartesianTransform = Tensor::Matrix<4, 4>;
using CartesianDerivative = Tensor::Matrix<Node::dim, count>;
using Gradient = Tensor::Matrix<Node::dim, Node::dim>;
using SurfaceGradient = Tensor::Matrix<Node::dim, Node::dim - 1>;
using Stress = Gradient;
using Stiffness = Tensor::Multi<count, Node::dim, count, Node::dim>;
using Positions = Fastor::Tensor<Scalar, count, Node::dim>;
using Forces = Fastor::Tensor<Scalar, count, Node::dim>;
using StiffnessResidual = std::tuple<Stiffness, Forces>;
}  // namespace Element

namespace SurfaceElement
{
static constexpr Tensor::Index count = Element::count - 1;
using Vtxhs = FeltElements::Vtxhs<count>;
using Positions = Fastor::Tensor<Scalar, count, Node::dim>;
}
}  // namespace FeltElements
