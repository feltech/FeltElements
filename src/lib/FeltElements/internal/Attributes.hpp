#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <string_view>

#include "FeltElements/Derivatives.hpp"
#include "FeltElements/Typedefs.hpp"

namespace FeltElements
{
namespace internal::Attribute
{
/*
 * Base traits
 */

template <typename T>
struct always_false : std::false_type
{
};
template <class Derived>
struct Traits
{
	static_assert(always_false<Derived>::value, "Traits are required");
};

template <class TData>
struct BaseTraits
{
	using Data = TData;
};

template <class Derived>
class Base
{
protected:
	using ThisTraits = Traits<Derived>;
	using Data = typename ThisTraits::Data;
	using Prop = typename ThisTraits::Prop;
	using Handle = typename ThisTraits::Handle;
	static constexpr std::string_view prop_name = ThisTraits::prop_name;

	explicit Base(Prop&& prop) : m_prop{prop} {}

	Data & operator[](Handle const & handle)
	{
		return m_prop[handle];
	}

	Data const & operator[](Handle const & handle) const
	{
		return m_prop[handle];
	}

	Prop m_prop;
};

template <typename TData>
struct VertexTraits : public BaseTraits<TData>
{
	using Handle = OpenVolumeMesh::VertexHandle;
	using Prop = OpenVolumeMesh::VertexPropertyT<TData>;
};

template <class Derived>
class Vertex : protected Base<Derived>
{
	using ThisBase = Base<Derived>;

protected:
	using Data = typename ThisBase::Data;
	explicit Vertex(Mesh& mesh)
		: ThisBase{mesh.request_vertex_property<Data>(ThisBase::prop_name.data())}
	{
	}
};

template <typename TData>
struct CellTraits : public BaseTraits<TData>
{
	using Handle = OpenVolumeMesh::CellHandle;
	using Prop = OpenVolumeMesh::CellPropertyT<TData>;
};

template <class Derived>
class Cell : protected Base<Derived>
{
	using ThisBase = Base<Derived>;
	using Data = typename ThisBase::Data;

protected:
	explicit Cell(Mesh& mesh)
		: ThisBase{mesh.request_cell_property<Data>(ThisBase::prop_name.data())}
	{
	}
};
}  // namespace internal::Attribute

/*
 * Forward declarations
 */
namespace Node::Attribute
{
class SpatialPosition;
class Force;
}  // namespace Node::Attribute
namespace Element::Attribute
{
class VertexHandles;
class MaterialShapeDerivative;
}  // namespace Element::Attribute

/*
 * Specific traits.
 */
namespace internal::Attribute
{
using namespace Node::Attribute;

template <>
struct Traits<SpatialPosition> : public VertexTraits<Node::Pos>
{
	static constexpr std::string_view prop_name = "spatial_position";
};

template <>
struct Traits<Force> : public VertexTraits<Node::Force>
{
	static constexpr std::string_view prop_name = "force";
};

using namespace Element::Attribute;
template <>
struct Traits<VertexHandles>
	: public CellTraits<std::array<OpenVolumeMesh::VertexHandle, Node::count>>
{
	static constexpr std::string_view prop_name = "vertices";
};

template <>
struct Traits<MaterialShapeDerivative> : public CellTraits<Element::ShapeDerivative>
{
	static constexpr std::string_view prop_name = "material_shape_derivative";
};
}  // namespace internal::Attribute
}  // namespace FeltElements