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
	[[nodiscard]] Node::Positions for_element(Element::Vtxhs const & vtxhs) const;
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

class NodalForces : private internal::Attribute::Cell<NodalForces>
{
	using ThisBase = internal::Attribute::Cell<NodalForces>;

public:
	explicit NodalForces(Mesh& mesh);
	using ThisBase::operator[];
	static IsoCoordDerivative const dL_by_dN;
	static ShapeDerivative const dN_by_dL;
};

class Stiffness : private internal::Attribute::Cell<Stiffness>
{
	using ThisBase = internal::Attribute::Cell<Stiffness>;

public:
	explicit Stiffness(Mesh& mesh);
	using ThisBase::operator[];
	static IsoCoordDerivative const dL_by_dN;
	static ShapeDerivative const dN_by_dL;
};
}  // namespace Element::Attribute
}  // namespace FeltElements
