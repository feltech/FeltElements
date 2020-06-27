#pragma once

#include "internal/Attributes.hpp"

namespace FeltElements
{
namespace Attribute
{
/*
 * Forward declarations
 */
namespace Vertex
{
class MaterialPosition;
class SpatialPosition;
class FixedDOF;
}  // namespace Vertex
namespace Cell
{
class VertexHandles;
class MaterialShapeDerivative;
class NodalForces;
class Stiffness;
}  // namespace Cell

namespace internal
{
template <>
struct Traits<Vertex::MaterialPosition> : public VertexTraits<Node::Pos>
{
	static constexpr std::string_view prop_name = "material_position";
};
template <>
struct Traits<Vertex::SpatialPosition> : public VertexTraits<Node::Pos>
{
	static constexpr std::string_view prop_name = "spatial_position";
};
template <>
struct Traits<Vertex::FixedDOF> : public VertexTraits<Node::Pos>
{
	static constexpr std::string_view prop_name = "spatial_position";
};

using namespace Attribute::Cell;
template <>
struct Traits<VertexHandles> : public CellTraits<Element::Vtxhs>
{
	static constexpr std::string_view prop_name = "vertices";
};

template <>
struct Traits<NodalForces> : public CellTraits<Node::Forces>
{
	static constexpr std::string_view prop_name = "nodal_forces";
};

template <>
struct Traits<Stiffness> : public CellTraits<Element::Stiffness>
{
	static constexpr std::string_view prop_name = "stiffness";
};

template <>
struct Traits<MaterialShapeDerivative> : public CellTraits<Element::ShapeDerivative>
{
	static constexpr std::string_view prop_name = "material_shape_derivative";
};
}  // namespace internal

namespace Vertex
{
class MaterialPosition final : private internal::VertexPositionBase<MaterialPosition>
{
	using ThisBase = internal::VertexPositionBase<MaterialPosition>;

public:
	explicit MaterialPosition(Mesh& mesh);
	using ThisBase::operator[];
	using ThisBase::for_element;
};

class SpatialPosition final : private internal::VertexPositionBase<SpatialPosition>
{
	using ThisBase = internal::VertexPositionBase<SpatialPosition>;

public:
	explicit SpatialPosition(Mesh& mesh);
	using ThisBase::operator[];
	using ThisBase::for_element;
};

class FixedDOF final : private internal::VertexBase<FixedDOF>
{
	using ThisBase = internal::VertexBase<FixedDOF>;

public:
	explicit FixedDOF(Mesh& mesh);
	using ThisBase::operator[];
};
}  // namespace Vertex

namespace Cell
{
class VertexHandles final : private internal::CellBase<VertexHandles>
{
	using ThisBase = internal::CellBase<VertexHandles>;

public:
	explicit VertexHandles(Mesh& mesh);
	using ThisBase::operator[];
};

class MaterialShapeDerivative final : private internal::CellBase<MaterialShapeDerivative>
{
	using ThisBase = internal::CellBase<MaterialShapeDerivative>;

public:
	explicit MaterialShapeDerivative(
		Mesh& mesh, VertexHandles const& vtxhs, Vertex::MaterialPosition const& X);
	using ThisBase::operator[];
	static Node::Positions const X;
	static Element::IsoCoordDerivative const dL_by_dN;
	static Element::ShapeDerivative const dN_by_dL;
};

class NodalForces final : private internal::CellBase<NodalForces>
{
	using ThisBase = internal::CellBase<NodalForces>;

public:
	explicit NodalForces(Mesh& mesh);
	using ThisBase::operator[];
};

class Stiffness final : private internal::CellBase<Stiffness>
{
	using ThisBase = internal::CellBase<Stiffness>;

public:
	explicit Stiffness(Mesh& mesh);
	using ThisBase::operator[];
};
}  // namespace Cell
}  // namespace Attribute

struct Attributes final
{
	explicit Attributes(Mesh& mesh);
	Attribute::Vertex::SpatialPosition x;
	Attribute::Vertex::MaterialPosition const X;
	Attribute::Cell::VertexHandles const vtxh;
	Attribute::Cell::MaterialShapeDerivative const dN_by_dX;
	Attribute::Cell::NodalForces T;
	Attribute::Cell::Stiffness K;
};

}  // namespace FeltElements
