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
struct always_false : std::false_type {};
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


template <class Derived>
class Position : protected Vertex<Derived>
{
	using ThisBase = Vertex<Derived>;
public:
	explicit Position(Mesh& mesh) : ThisBase{mesh}
	{
		for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
		{
			Mesh::PointT vtx = mesh.vertex(*itvtxh);
			(*this)[*itvtxh] =
				Tensor::BaseMap<Mesh::PointT::value_type, Mesh::PointT::size()>{vtx.data()};
		}
	}

	[[nodiscard]] Node::Positions for_element(Element::Vtxhs const & vtxhs) const
	{
		using Tensor::Func::all;
		Node::Positions x;
		for (Tensor::Index node_idx = 0; node_idx < Node::Positions::dimension(0); node_idx++)
			x(node_idx, all) = ThisBase::m_prop[vtxhs[node_idx]];

		return x;
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

protected:
	using Data = typename ThisBase::Data;
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
class MaterialPosition;
class SpatialPosition;
}  // namespace Node::Attribute
namespace Element::Attribute
{
class VertexHandles;
class MaterialShapeDerivative;
class NodalForces;
class Stiffness;
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
struct Traits<MaterialPosition> : public VertexTraits<Node::Pos>
{
	static constexpr std::string_view prop_name = "material_position";
};

using namespace Element::Attribute;
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
}  // namespace internal::Attribute
}  // namespace FeltElements