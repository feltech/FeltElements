#pragma once
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>

#include "FeltElements/Tetrahedron.hpp"
#include "FeltElements/internal/Attribute.hpp"

namespace FeltElements::Attribute
{
namespace Node
{
class SpatialPosition;
class Force;
}
namespace Element
{
class VertexHandles;
class MaterialShapeDerivative;
}

template <> struct internal::Traits<Node::SpatialPosition> :
	public internal::VertexTraits<Tetrahedron::Node::Pos>
{
	static constexpr std::string_view prop_name = "spatial_position";
};

template <> struct internal::Traits<Node::Force> :
	public internal::VertexTraits<Tetrahedron::Node::Force>
{
	static constexpr std::string_view prop_name = "force";
};

template <> struct internal::Traits<Element::VertexHandles> :
	public internal::CellTraits<std::array<OpenVolumeMesh::VertexHandle, Tetrahedron::Node::count>>
{
	static constexpr std::string_view prop_name = "vertices";
};

template <> struct internal::Traits<Element::MaterialShapeDerivative> :
	public internal::CellTraits<Tetrahedron::ShapeDerivativeTensor>
{
	static constexpr std::string_view prop_name = "material_shape_derivative";
};

class Node::SpatialPosition : private internal::Vertex<SpatialPosition>
{
	using ThisBase = internal::Vertex<SpatialPosition>;
	using Mesh = typename ThisBase::Mesh;
public:
	explicit SpatialPosition(Mesh & mesh);
	using ThisBase::operator[];
};

class Node::Force : private internal::Vertex<Force>
{
	using ThisBase = internal::Vertex<Force>;
	using Mesh = typename ThisBase::Mesh;
public:
	explicit Force(Mesh & mesh);
	using ThisBase::operator[];
};

class Element::VertexHandles : private internal::Cell<VertexHandles>
{
	using ThisBase = internal::Cell<VertexHandles>;
	using Mesh = typename ThisBase::Mesh;
public:
	explicit VertexHandles(Mesh & mesh);
	using ThisBase::operator[];
};

class Element::MaterialShapeDerivative : private internal::Cell<MaterialShapeDerivative>
{
	using ThisBase = internal::Cell<MaterialShapeDerivative>;
	using Mesh = typename ThisBase::Mesh;
public:
	explicit MaterialShapeDerivative(Mesh & mesh, VertexHandles const & vtxhs);
	using ThisBase::operator[];
};

} // namespace FeltElements::Attribute
