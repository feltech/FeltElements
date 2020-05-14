#include <string_view>

#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <unsupported/Eigen/CXX11/Tensor>

#include "FeltElements/Tetrahedron.hpp"

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

template <class Derived> struct Traits
{
	static_assert(sizeof(Derived) != sizeof(Derived), "Traits are required");
};

template <class TData> struct BaseTraits
{
	using Mesh = OpenVolumeMesh::GeometricTetrahedralMeshV3d;
	using Data = TData;
};

template <class Derived> class Base
{
protected:
	using ThisTraits = Traits<Derived>;
	using Mesh = typename ThisTraits::Mesh;
	using Data = typename ThisTraits::Data;
	using Prop = typename ThisTraits::Prop;
	using Handle = typename ThisTraits::Handle;
	static constexpr std::string_view prop_name = ThisTraits::prop_name;

	explicit Base(Prop && prop) : m_prop{prop} {}

	Data const & operator[](Handle const & handle) const
	{
		return m_prop[handle];
	}

	Data & operator[](Handle const & handle)
	{
		return m_prop[handle];
	}

	Prop m_prop;
};

template <typename TData> struct VertexTraits : public BaseTraits<TData>
{
	using Handle = OpenVolumeMesh::VertexHandle;
	using Prop = OpenVolumeMesh::VertexPropertyT<TData>;
};

template <class Derived> class Vertex : protected Base<Derived>
{
	using ThisBase = Base<Derived>;

protected:
	using Data = typename ThisBase::Data;
	explicit Vertex(typename ThisBase::Mesh & mesh)
		: ThisBase{mesh.template request_vertex_property<Data>(ThisBase::prop_name.data())}
	{
	}
};

template <typename TData> struct CellTraits : public BaseTraits<TData>
{
	using Handle = OpenVolumeMesh::CellHandle;
	using Prop = OpenVolumeMesh::CellPropertyT<TData>;
};

template <class Derived> class Cell : protected Base<Derived>
{
	using ThisBase = Base<Derived>;
	using Data = typename ThisBase::Data;

protected:
	explicit Cell(typename ThisBase::Mesh & mesh)
		: ThisBase{mesh.template request_cell_property<Data>(ThisBase::prop_name.data())}
	{
	}
};

using namespace FeltElements::Attribute;

template <> struct Traits<Node::SpatialPosition> : public VertexTraits<Tetrahedron::Node::Pos>
{
	static constexpr std::string_view prop_name = "spatial_position";
};

template <> struct Traits<Node::Force> : public VertexTraits<Tetrahedron::Node::Force>
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
} // namespace FeltElements