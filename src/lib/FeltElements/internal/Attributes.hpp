#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <string_view>

#include "FeltElements/Derivatives.hpp"
#include "FeltElements/Typedefs.hpp"

namespace FeltElements::Attribute::internal
{
template <class Derived>
struct Traits
{
	// Enforce that this must be specialised.
	template <typename T>
	struct missing_traits : std::false_type
	{
	};
	static_assert(missing_traits<Derived>::value, "Traits are required");
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

	Data& operator[](Handle const& handle)
	{
		return m_prop[handle];
	}

	Data const& operator[](Handle const& handle) const
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
class VertexBase : protected Base<Derived>
{
	using ThisBase = Base<Derived>;

protected:
	using Data = typename ThisBase::Data;
	explicit VertexBase(Mesh& mesh)
		: ThisBase{mesh.request_vertex_property<Data>(ThisBase::prop_name.data())}
	{
	}
};

template <class Derived>
class VertexPositionBase : protected VertexBase<Derived>
{
	using ThisBase = VertexBase<Derived>;

protected:
	explicit VertexPositionBase(Mesh& mesh) : ThisBase{mesh}
	{
		for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
		{
			Mesh::PointT vtx = mesh.vertex(*itvtxh);
			(*this)[*itvtxh] =
				Tensor::BaseMap<Mesh::PointT::value_type, Mesh::PointT::size()>{vtx.data()};
		}
	}

	[[nodiscard]] Node::Positions for_element(Element::Vtxhs const& vtxhs) const
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
class CellBase : protected Base<Derived>
{
	using ThisBase = Base<Derived>;

protected:
	using Data = typename ThisBase::Data;
	explicit CellBase(Mesh& mesh)
		: ThisBase{mesh.request_cell_property<Data>(ThisBase::prop_name.data())}
	{
	}
};
}  // namespace FeltElements::Attribute::internal