#pragma once

#include "internal/Attributes.hpp"

namespace FeltElements
{
namespace Node::Attribute
{
class SpatialPosition : private internal::Attribute::Vertex<SpatialPosition>
{
	using ThisBase = internal::Attribute::Vertex<SpatialPosition>;

public:
	explicit SpatialPosition(Mesh& mesh);
	using ThisBase::operator[];
};

class Force : private internal::Attribute::Vertex<Force>
{
	using ThisBase = internal::Attribute::Vertex<Force>;

public:
	explicit Force(Mesh& mesh);
	using ThisBase::operator[];
};
}  // namespace Node::Attribute

namespace Element::Attribute
{
class VertexHandles : private internal::Attribute::Cell<VertexHandles>
{
	using ThisBase = internal::Attribute::Cell<VertexHandles>;

public:
	explicit VertexHandles(Mesh& mesh);
	using ThisBase::operator[];
};

class MaterialShapeDerivative : private internal::Attribute::Cell<MaterialShapeDerivative>
{
	using ThisBase = internal::Attribute::Cell<MaterialShapeDerivative>;

public:
	explicit MaterialShapeDerivative(Mesh& mesh, VertexHandles const& vtxhs);
	using ThisBase::operator[];
	static IsoCoordDerivative const dL_by_dN;
	static ShapeDerivative const dN_by_dL;
};
}  // namespace Element::Attribute
}  // namespace FeltElements
