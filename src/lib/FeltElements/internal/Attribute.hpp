#pragma once
#include <string_view>

#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <unsupported/Eigen/CXX11/Tensor>

namespace FeltElements::internal::Attribute
{

template <class Derived>
struct Traits
{
	static_assert(sizeof(Derived) != sizeof(Derived), "Traits are required");
};

template <class TData>
struct BaseTraits
{
	using Mesh = OpenVolumeMesh::GeometricTetrahedralMeshV3d;
	using Data = TData;
};

template <class Derived>
class Base
{
protected:
	using ThisTraits = Traits<Derived>;
	using Mesh = typename ThisTraits::Mesh;
	using Data = typename ThisTraits::Data;
	using Prop = typename ThisTraits::Prop;
	using Handle = typename ThisTraits::Handle;
	static constexpr std::string_view prop_name = ThisTraits::prop_name;

	explicit Base(Prop&& prop) : m_prop{prop} {}

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


template <typename TData>
struct VertexTraits : public BaseTraits<TData>
{
	using Handle = OpenVolumeMesh::VertexHandle;
	using Prop = OpenVolumeMesh::VertexPropertyT<TData>;
};


template <class Derived>
class Vertex: protected Base<Derived>
{
	using ThisBase = Base<Derived>;
protected:
	using Data = typename ThisBase::Data;
	explicit Vertex(typename ThisBase::Mesh & mesh) :
		ThisBase{mesh.template request_vertex_property<Data>(ThisBase::prop_name.data())}
	{}
};


template <typename TData>
struct CellTraits : public BaseTraits<TData>
{
	using Handle = OpenVolumeMesh::CellHandle;
	using Prop = OpenVolumeMesh::CellPropertyT<TData>;
};


template <class Derived>
class Cell: protected Base<Derived>
{
	using ThisBase = Base<Derived>;
	using Data = typename ThisBase::Data;
protected:
	explicit Cell(typename ThisBase::Mesh & mesh) :
		ThisBase{mesh.template request_cell_property<Data>(ThisBase::prop_name.data())}
	{}
};
}