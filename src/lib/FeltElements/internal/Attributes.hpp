#include <string_view>

#include <boost/range/iterator_range.hpp>
#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>

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

	explicit Base(Prop && prop) : m_prop{prop} {}

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
struct MeshTraits : public BaseTraits<TData>
{
	using Handle = OpenVolumeMesh::MeshHandle;
	using Prop = OpenVolumeMesh::MeshPropertyT<TData>;
};

template <class Derived>
class MeshBase : protected Base<Derived>
{
	using ThisBase = Base<Derived>;
	using Handle = typename ThisBase::Handle;

protected:
	using Data = typename ThisBase::Data;
	explicit MeshBase(Mesh & mesh)
		: ThisBase{mesh.request_mesh_property<Data>(ThisBase::prop_name.data())}
	{
	}

	Data const & operator*() const
	{
		return (*this)[Handle{0}];
	}

	Data & operator*()
	{
		return (*this)[Handle{0}];
	}

	Data const * operator->() const
	{
		return &(*this)[Handle{0}];
	}

	Data * operator->()
	{
		return &(*this)[Handle{0}];
	}
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
	explicit VertexBase(Mesh & mesh)
		: ThisBase{mesh.request_vertex_property<Data>(ThisBase::prop_name.data())}
	{
	}
};

template <class Derived>
class VertexPositionBase : protected VertexBase<Derived>
{
	using ThisBase = VertexBase<Derived>;

protected:
	explicit VertexPositionBase(Mesh & mesh) : ThisBase{mesh}
	{
		for (auto vtxh : boost::make_iterator_range(mesh.vertices()))
		{
			Vtx const & vtx = mesh.vertex(vtxh);
			(*this)[vtxh] = Tensor::BaseMap<Vtx::value_type const, Vtx::size()>{vtx.data()};
		}
	}

	template <std::size_t count>
	[[nodiscard]] Tensor::Matrix<count, Node::dim> for_element(Vtxhs<count> const & vtxhs) const
	{
		using Tensor::Func::all;
		using Positions = Tensor::Matrix<count, Node::dim>;
		Positions x;
		for (Tensor::Index node_idx = 0; node_idx < Positions::dimension(0); node_idx++)
			x(node_idx, all) = (*this)[vtxhs[node_idx]];

		return x;
	}
};

template <typename TData>
struct SurfaceTraits : public BaseTraits<TData>
{
	using Handle = OpenVolumeMesh::HalfFaceHandle;
	using Prop = OpenVolumeMesh::HalfFacePropertyT<TData>;
};

template <class Derived>
class SurfaceBase : protected Base<Derived>
{
	using ThisBase = Base<Derived>;

protected:
	using Data = typename ThisBase::Data;
	explicit SurfaceBase(Mesh & mesh)
		: ThisBase{mesh.request_halfface_property<Data>(ThisBase::prop_name.data())}
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
class CellBase : protected Base<Derived>
{
	using ThisBase = Base<Derived>;

protected:
	using Data = typename ThisBase::Data;
	explicit CellBase(Mesh & mesh)
		: ThisBase{mesh.request_cell_property<Data>(ThisBase::prop_name.data())}
	{
	}
};
}  // namespace FeltElements::Attribute::internal