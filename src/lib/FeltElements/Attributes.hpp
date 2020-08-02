#pragma once

#include "internal/Attributes.hpp"

namespace FeltElements
{
namespace Attribute
{
struct MaterialProperties
{
	/// Density.
	Scalar rho;
	/// Lame's first parameter.
	Scalar lambda;
	/// Lame's second parameter.
	Scalar mu;
};
/*
 * Forward declarations
 */
namespace Body
{
class Properties;
class Force;
class Surface;
}  // namespace Body
namespace Vertex
{
class MaterialPosition;
class SpatialPosition;
class FixedDOF;
}  // namespace Vertex
namespace Cell
{
class MaterialVolume;
class SpatialVolume;
class VertexHandles;
class MaterialShapeDerivative;
class NodalForces;
class Stiffness;
}  // namespace Cell

namespace internal
{
template <>
struct Traits<Body::Properties> : public MeshTraits<MaterialProperties>
{
	static constexpr std::string_view prop_name = "material_properties";
};
template <>
struct Traits<Body::Force> : public MeshTraits<Node::Force>
{
	static constexpr std::string_view prop_name = "body_force";
};
template <>
struct Traits<Body::Surface> : public MeshTraits<std::vector<Halffaceh>>
{
	static constexpr std::string_view prop_name = "surface_halffaces";
};
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
	static constexpr std::string_view prop_name = "fixed_dof";
};
template <>
struct Traits<Cell::MaterialVolume> : public CellTraits<Scalar>
{
	static constexpr std::string_view prop_name = "material_volume";
};
template <>
struct Traits<Cell::SpatialVolume> : public CellTraits<Scalar>
{
	static constexpr std::string_view prop_name = "spatial_volume";
};
template <>
struct Traits<Cell::VertexHandles> : public CellTraits<Element::Vtxhs>
{
	static constexpr std::string_view prop_name = "vertices";
};
template <>
struct Traits<Cell::NodalForces> : public CellTraits<Node::Forces>
{
	static constexpr std::string_view prop_name = "nodal_forces";
};

template <>
struct Traits<Cell::Stiffness> : public CellTraits<Element::Stiffness>
{
	static constexpr std::string_view prop_name = "stiffness";
};

template <>
struct Traits<Cell::MaterialShapeDerivative> : public CellTraits<Element::ShapeDerivative>
{
	static constexpr std::string_view prop_name = "material_shape_derivative";
};
}  // namespace internal

namespace Body
{
class Properties final : private internal::MeshBase<Properties>
{
	using ThisBase = internal::MeshBase<Properties>;

public:
	explicit Properties(Mesh & mesh);
	using ThisBase::operator*;
	using ThisBase::operator->;
};

class Force final : private internal::MeshBase<Force>
{
	using ThisBase = internal::MeshBase<Force>;

public:
	explicit Force(Mesh & mesh);
	using ThisBase::operator*;
	using ThisBase::operator->;
};

class Surface final : private internal::MeshBase<Surface>
{
	using ThisBase = internal::MeshBase<Surface>;

public:
	explicit Surface(Mesh & mesh);
	using ThisBase::operator*;
	using ThisBase::operator->;
};
}  // namespace Body

namespace Vertex
{
class MaterialPosition final : private internal::VertexPositionBase<MaterialPosition>
{
	using ThisBase = internal::VertexPositionBase<MaterialPosition>;

public:
	explicit MaterialPosition(Mesh & mesh);
	using ThisBase::operator[];
	using ThisBase::for_element;
};

class SpatialPosition final : private internal::VertexPositionBase<SpatialPosition>
{
	using ThisBase = internal::VertexPositionBase<SpatialPosition>;

public:
	explicit SpatialPosition(Mesh & mesh);
	using ThisBase::operator[];
	using ThisBase::for_element;
};

class FixedDOF final : private internal::VertexBase<FixedDOF>
{
	using ThisBase = internal::VertexBase<FixedDOF>;

public:
	explicit FixedDOF(Mesh & mesh);
	using ThisBase::operator[];
};
}  // namespace Vertex

namespace Cell
{
class MaterialVolume final : private internal::CellBase<MaterialVolume>
{
	using ThisBase = internal::CellBase<MaterialVolume>;

public:
	explicit MaterialVolume(
		Mesh & mesh, Cell::VertexHandles const & vtxh, Vertex::MaterialPosition const & X);
	using ThisBase::operator[];
};

class SpatialVolume final : private internal::CellBase<SpatialVolume>
{
	using ThisBase = internal::CellBase<SpatialVolume>;

public:
	explicit SpatialVolume(Mesh & mesh);
	using ThisBase::operator[];
};

class VertexHandles final : private internal::CellBase<VertexHandles>
{
	using ThisBase = internal::CellBase<VertexHandles>;

public:
	explicit VertexHandles(Mesh & mesh);
	using ThisBase::operator[];
};

class MaterialShapeDerivative final : private internal::CellBase<MaterialShapeDerivative>
{
	using ThisBase = internal::CellBase<MaterialShapeDerivative>;

public:
	explicit MaterialShapeDerivative(
		Mesh & mesh, VertexHandles const & vtxhs, Vertex::MaterialPosition const & X);
	using ThisBase::operator[];
	static Node::Positions const X;
};

class NodalForces final : private internal::CellBase<NodalForces>
{
	using ThisBase = internal::CellBase<NodalForces>;

public:
	explicit NodalForces(Mesh & mesh);
	using ThisBase::operator[];
};

class Stiffness final : private internal::CellBase<Stiffness>
{
	using ThisBase = internal::CellBase<Stiffness>;

public:
	explicit Stiffness(Mesh & mesh);
	using ThisBase::operator[];
};
}  // namespace Cell
}  // namespace Attribute

struct Attributes final
{
	explicit Attributes(Mesh & mesh);
	Attribute::Body::Properties material;
	Attribute::Body::Force f;
	Attribute::Vertex::SpatialPosition x;
	Attribute::Vertex::MaterialPosition const X;
	Attribute::Vertex::FixedDOF fixed_dof;
	Attribute::Cell::VertexHandles const vtxh;
	Attribute::Cell::MaterialVolume const V;
	Attribute::Cell::SpatialVolume v;
	Attribute::Cell::MaterialShapeDerivative const dN_by_dX;
	Attribute::Cell::NodalForces T;
	Attribute::Cell::Stiffness K;
};

}  // namespace FeltElements
