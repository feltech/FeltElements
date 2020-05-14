#pragma once

#include "FeltElements/internal/Attributes.hpp"

namespace FeltElements::Attribute
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

} // namespace FeltElements::Attribute
