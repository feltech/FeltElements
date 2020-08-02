#pragma once
#include <Fastor/Fastor.h>

#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <boost/range/iterator_range_core.hpp>

namespace FeltElements
{
using Mesh = OpenVolumeMesh::GeometricTetrahedralMeshV3d;
using Vtxh = OpenVolumeMesh::VertexHandle;
using Cellh = OpenVolumeMesh::CellHandle;
using Halffaceh = OpenVolumeMesh::HalfFaceHandle;
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
static constexpr Tensor::Index count = 4;
static constexpr Tensor::Index dim = 3;

using Pos = Tensor::Vector<dim>;
using PosMap = Tensor::Map<dim>;
using Positions = Fastor::Tensor<Scalar, count, dim>;
using SurfacePositions = Fastor::Tensor<Scalar, count - 1, dim>;
using Force = Tensor::Vector<dim>;
using Forces = Fastor::Tensor<Scalar, count, dim>;
}  // namespace Node

namespace Element
{
using Vtxhs = std::array<OpenVolumeMesh::VertexHandle, Node::count>;
using IsoCoordDerivative = Tensor::Matrix<Node::dim, Node::count>;
using ShapeDerivative = Tensor::Matrix<Node::count, Node::dim>;
using SurfaceShapeDerivative = Tensor::Matrix<Node::count - 1, Node::dim - 1>;
using ShapeDerivativeDeterminant = Tensor::Multi<Node::count, Node::count, Node::count>;
using Elasticity = Tensor::Multi<Node::dim, Node::dim, Node::dim, Node::dim>;
using ShapeCartesianTransform = Tensor::Matrix<4, 4>;
using CartesianDerivative = Tensor::Matrix<Node::dim, Node::count>;
using Gradient = Tensor::Matrix<Node::dim, Node::dim>;
using SurfaceGradient = Tensor::Matrix<Node::dim, Node::dim - 1>;
using Stress = Gradient;
using Stiffness = Tensor::Multi<Node::count, Node::dim, Node::count, Node::dim>;
using StiffnessForcesVolume = std::tuple<Stiffness, Node::Forces, Scalar>;
}  // namespace Element

namespace SurfaceElement
{
using Vtxhs = std::array<OpenVolumeMesh::VertexHandle, Node::count - 1>;
}
}  // namespace FeltElements
