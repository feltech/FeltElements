#pragma once

#include "FeltElements/Tetrahedron.hpp"
#include "FeltElements/internal/Attribute.hpp"

namespace FeltElements
{
namespace Attribute
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
} // namespace Attribute


namespace internal::Attribute
{
using namespace FeltElements::Attribute;

template <>
struct Traits<Node::SpatialPosition>
	: public VertexTraits<Tetrahedron::Node::Pos>
{
	static constexpr std::string_view prop_name = "spatial_position";
};

template <>
struct Traits<Node::Force>
	: public VertexTraits<Tetrahedron::Node::Force>
{
	static constexpr std::string_view prop_name = "force";
};

template <>
struct Traits<Element::VertexHandles>
	: public CellTraits<std::array<OpenVolumeMesh::VertexHandle, Tetrahedron::Node::count>>
{
	static constexpr std::string_view prop_name = "vertices";
};

template <>
struct Traits<Element::MaterialShapeDerivative>
	: public CellTraits<Tetrahedron::ShapeDerivativeTensor>
{
	static constexpr std::string_view prop_name = "material_shape_derivative";
};
} // namespace internal::Attribute


namespace Attribute
{
class Node::SpatialPosition : private internal::Attribute::Vertex<SpatialPosition>
{
	using ThisBase = internal::Attribute::Vertex<SpatialPosition>;
	using Mesh = typename ThisBase::Mesh;

public:
	explicit SpatialPosition(Mesh & mesh);
	using ThisBase::operator[];
};

class Node::Force : private internal::Attribute::Vertex<Force>
{
	using ThisBase = internal::Attribute::Vertex<Force>;
	using Mesh = typename ThisBase::Mesh;

public:
	explicit Force(Mesh & mesh);
	using ThisBase::operator[];
};

class Element::VertexHandles : private internal::Attribute::Cell<VertexHandles>
{
	using ThisBase = internal::Attribute::Cell<VertexHandles>;
	using Mesh = typename ThisBase::Mesh;

public:
	explicit VertexHandles(Mesh & mesh);
	using ThisBase::operator[];
};

class Element::MaterialShapeDerivative : private internal::Attribute::Cell<MaterialShapeDerivative>
{
	using ThisBase = internal::Attribute::Cell<MaterialShapeDerivative>;
	using Mesh = typename ThisBase::Mesh;

public:
	explicit MaterialShapeDerivative(Mesh & mesh, VertexHandles const & vtxhs);
	using ThisBase::operator[];
	static Tetrahedron::IsoCoordDerivativeTensor const dL_by_dN;
	static Tetrahedron::ShapeDerivativeTensor const dN_by_dL;
};

} // namespace Attribute
} // namespace FeltElements