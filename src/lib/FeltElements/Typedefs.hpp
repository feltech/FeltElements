#pragma once
#include <Fastor/Fastor.h>

#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <boost/container/static_vector.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/range/iterator.hpp>

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
static constexpr Tensor::Index num_nodes = 4;
static constexpr Tensor::Index num_faces = 4;

using Vtxhs = FeltElements::Vtxhs<num_nodes>;
using IsoCoordDerivative = Tensor::Matrix<Node::dim, num_nodes>;
using ShapeDerivative = Tensor::Matrix<num_nodes, Node::dim>;
using SurfaceShapeDerivative = Tensor::Matrix<num_nodes - 1, Node::dim - 1>;
using ShapeDerivativeDeterminant = Tensor::Multi<num_nodes, num_nodes, num_nodes>;
using Elasticity = Tensor::Multi<Node::dim, Node::dim, Node::dim, Node::dim>;
using ShapeCartesianTransform = Tensor::Matrix<4, 4>;
using CartesianDerivative = Tensor::Matrix<Node::dim, num_nodes>;
using Gradient = Tensor::Matrix<Node::dim, Node::dim>;
using SurfaceGradient = Tensor::Matrix<Node::dim, Node::dim - 1>;
using Stress = Gradient;
using Stiffness = Tensor::Multi<num_nodes, Node::dim, num_nodes, Node::dim>;
using NodePositions = Fastor::Tensor<Scalar, num_nodes, Node::dim>;
using Forces = Fastor::Tensor<Scalar, num_nodes, Node::dim>;
using StiffnessResidual = std::tuple<Stiffness, Forces>;
}  // namespace Element

namespace BoundaryElement
{
static constexpr Tensor::Index num_nodes = Element::num_nodes - 1;

using Vtxhs = FeltElements::Vtxhs<num_nodes>;
using VtxhIdxs = std::array<Tensor::Index, BoundaryElement::num_nodes>;
using NodePositions = Tensor::Matrix<num_nodes, Node::dim>;
using Stiffness = Tensor::Multi<num_nodes, Node::dim, num_nodes, Node::dim>;
}  // namespace BoundaryElement

namespace Element
{
template <class T>
using PerFace = boost::container::static_vector<T, Element::num_faces>;
using BoundaryVtxhIdxs = PerFace<BoundaryElement::VtxhIdxs>;
using BoundaryNodePositions = PerFace<BoundaryElement::NodePositions>;
}  // namespace Element

template <typename Index = std::size_t, class Haystack, typename Needle>
constexpr Index index_of(Haystack && haystack, Needle && needle)
{
	auto const & it = boost::range::find(
		std::forward<decltype(haystack)>(haystack), std::forward<decltype(needle)>(needle));
	return static_cast<Index>(std::distance(haystack.begin(), it));
}
}  // namespace FeltElements
