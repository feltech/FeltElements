#pragma once

#include "Body.hpp"
#include "internal/Attributes.hpp"

namespace FeltElements
{
namespace Attribute
{
/*
 * Forward declarations
 */
namespace MeshBody
{
class MaterialProperties;
class Forces;
class Surface;
}  // namespace MeshBody
namespace Vertex
{
class MaterialPosition;
class SpatialPosition;
class FixedDOF;
}  // namespace Vertex
namespace Surface
{
class Traction;
}
namespace Cell
{
class VertexHandles;
class MaterialShapeDerivative;
class NodalForces;
class Stiffness;
class Boundary;
}  // namespace Cell

namespace internal
{
template <>
struct Traits<MeshBody::MaterialProperties> : public MeshTraits<Body::Material>
{
	static constexpr std::string_view prop_name = "material_properties";
};
template <>
struct Traits<MeshBody::Forces> : public MeshTraits<Body::Forces>
{
	static constexpr std::string_view prop_name = "body_forces";
};
template <>
struct Traits<MeshBody::Surface> : public MeshTraits<std::vector<BoundaryElement::Vtxhs>>
{
	static constexpr std::string_view prop_name = "surface_vertices";
};
template <>
struct Traits<Surface::Traction> : public SurfaceTraits<Node::Force>
{
	static constexpr std::string_view prop_name = "traction";
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
struct Traits<Cell::VertexHandles> : public CellTraits<Element::Vtxhs>
{
	static constexpr std::string_view prop_name = "vertices";
};
template <>
struct Traits<Cell::NodalForces> : public CellTraits<Element::Forces>
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
template <>
struct Traits<Cell::Boundary> : public CellTraits<Element::Boundary>
{
	static constexpr std::string_view prop_name = "element_boundary";
};
}  // namespace internal

namespace MeshBody
{
class MaterialProperties final : private internal::MeshBase<MaterialProperties>
{
	using ThisBase = internal::MeshBase<MaterialProperties>;

public:
	explicit MaterialProperties(Mesh & mesh);
	using ThisBase::operator*;
	using ThisBase::operator->;
};

class Forces final : private internal::MeshBase<Forces>
{
	using ThisBase = internal::MeshBase<Forces>;

public:
	explicit Forces(Mesh & mesh);
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
}  // namespace MeshBody

namespace Surface
{
class Traction final : private internal::SurfaceBase<Traction>
{
	using ThisBase = internal::SurfaceBase<Traction>;

public:
	explicit Traction(Mesh & mesh);
	using ThisBase::operator[];
};
}  // namespace Surface

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
		Mesh & mesh, VertexHandles const & vtxhs, Vertex::MaterialPosition const & X_);
	using ThisBase::operator[];
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

class Boundary final : private internal::CellBase<Boundary>
{
	using ThisBase = internal::CellBase<Boundary>;

public:
	explicit Boundary(Mesh & mesh);
	using ThisBase::operator[];
};
}  // namespace Cell
}  // namespace Attribute

struct Attributes final
{
	explicit Attributes(Mesh & mesh);
	Attribute::MeshBody::MaterialProperties material;
	Attribute::MeshBody::Forces forces;
	Attribute::MeshBody::Surface surface_vtxh;
	Attribute::Vertex::SpatialPosition x;
	Attribute::Vertex::MaterialPosition const X;
	Attribute::Vertex::FixedDOF fixed_dof;
	Attribute::Cell::VertexHandles const vtxh;
	Attribute::Cell::MaterialShapeDerivative const dN_by_dX;
	Attribute::Cell::NodalForces R;
	Attribute::Cell::Stiffness K;
};

}  // namespace FeltElements
