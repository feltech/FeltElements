#pragma once

#include "internal/Attributes.hpp"

namespace FeltElements
{
namespace Node::Attribute
{

class SpatialPosition final : private internal::Attribute::Position<SpatialPosition>
{
	using ThisBase = internal::Attribute::Position<SpatialPosition>;
public:
	using ThisBase::Position;
	using ThisBase::operator[];
	using ThisBase::for_element;
};

class MaterialPosition final : private internal::Attribute::Position<MaterialPosition>
{
	using ThisBase = internal::Attribute::Position<MaterialPosition>;
public:
	using ThisBase::Position;
	using ThisBase::operator[];
	using ThisBase::for_element;
};

}  // namespace Node::Attribute

namespace Element::Attribute
{
class VertexHandles final : private internal::Attribute::Cell<VertexHandles>
{
	using ThisBase = internal::Attribute::Cell<VertexHandles>;

public:
	explicit VertexHandles(Mesh& mesh);
	using ThisBase::operator[];
};

class MaterialShapeDerivative final : private internal::Attribute::Cell<MaterialShapeDerivative>
{
	using ThisBase = internal::Attribute::Cell<MaterialShapeDerivative>;

public:
	explicit MaterialShapeDerivative(
		Mesh& mesh, VertexHandles const & vtxhs, Node::Attribute::MaterialPosition const& X);
	using ThisBase::operator[];
	static Node::Positions const X;
	static IsoCoordDerivative const dL_by_dN;
	static ShapeDerivative const dN_by_dL;
};

class NodalForces final : private internal::Attribute::Cell<NodalForces>
{
	using ThisBase = internal::Attribute::Cell<NodalForces>;

public:
	explicit NodalForces(Mesh& mesh);
	using ThisBase::operator[];
};

class Stiffness final : private internal::Attribute::Cell<Stiffness>
{
	using ThisBase = internal::Attribute::Cell<Stiffness>;

public:
	explicit Stiffness(Mesh& mesh);
	using ThisBase::operator[];
};
}  // namespace Element::Attribute

struct Attributes
{
	explicit Attributes(Mesh& mesh);
	Node::Attribute::SpatialPosition x;
	Node::Attribute::MaterialPosition const X;
	Element::Attribute::VertexHandles const vtxh;
	Element::Attribute::MaterialShapeDerivative const dN_by_dX;
	Element::Attribute::NodalForces T;
	Element::Attribute::Stiffness K;
};

}  // namespace FeltElements
