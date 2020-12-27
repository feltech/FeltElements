#include "Attributes.hpp"

#include <boost/range.hpp>

#include "Derivatives.hpp"

namespace FeltElements
{
namespace Attribute
{
namespace MeshBody
{
MaterialProperties::MaterialProperties(Mesh & mesh) : ThisBase(mesh)
{
	(*(*this)) = Body::Material{0, 0, 0};
}

Surface::Surface(Mesh & mesh) : ThisBase(mesh)
{
	for (auto ithfh = mesh.bhf_iter(); ithfh.valid(); ithfh++)
	{
		SurfaceElement::Vtxhs vtxhs;
		Tensor::Index vtx_idx = 0;
		for (auto halfedgeh : boost::make_iterator_range(mesh.halfface_halfedges(*ithfh)))
			vtxhs[vtx_idx++] = mesh.halfedge(halfedgeh).from_vertex();
		(*this)->push_back(vtxhs);
	}
}

Forces::Forces(Mesh & mesh) : MeshBase(mesh)
{
	(*(*this)) = Body::Forces{0, {0, 0, 0}};
}
}  // namespace MeshBody

namespace Surface
{
Traction::Traction(Mesh & mesh) : ThisBase{mesh}
{
	for (auto hfh : boost::make_iterator_range(mesh.halffaces())) (*this)[hfh] = 0;
}
}  // namespace Surface

namespace Vertex
{
MaterialPosition::MaterialPosition(Mesh & mesh) : ThisBase(mesh) {}
SpatialPosition::SpatialPosition(Mesh & mesh) : ThisBase(mesh) {}

FixedDOF::FixedDOF(Mesh & mesh) : ThisBase(mesh)
{
	for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
		(*this)[*itvtxh] = 0;
}
}  // namespace Vertex

namespace Cell
{
VertexHandles::VertexHandles(Mesh & mesh) : ThisBase{mesh}
{
	for (auto cellh : boost::make_iterator_range(mesh.cells()))
	{
		auto const & vtxhs = mesh.get_cell_vertices(cellh);
		std::copy_n(vtxhs.begin(), vtxhs.size(), (*this)[cellh].begin());
	}
}

MaterialShapeDerivative::MaterialShapeDerivative(
	Mesh & mesh, VertexHandles const & vtxhs, Vertex::MaterialPosition const & X_)
	: ThisBase{mesh}
{
	for (auto cellh : boost::make_iterator_range(mesh.cells()))
		(*this)[cellh] = Derivatives::dN_by_dX(X_.for_element(vtxhs[cellh]));
}

NodalForces::NodalForces(Mesh & mesh) : ThisBase{mesh}
{
	for (auto cellh : boost::make_iterator_range(mesh.cells()))
		(*this)[cellh] = 0;
}

Stiffness::Stiffness(Mesh & mesh) : ThisBase(mesh)
{
	for (auto cellh : boost::make_iterator_range(mesh.cells()))
		(*this)[cellh] = 0;
}
}  // namespace Cell
}  // namespace Attribute

Attributes::Attributes(Mesh & mesh)
	: material{mesh},
	  forces{mesh},
	  surface_vtxh{mesh},
	  x{mesh},
	  X{mesh},
	  fixed_dof(mesh),
	  vtxh{mesh},
	  dN_by_dX{mesh, vtxh, X},
	  R{mesh},
	  K{mesh}
{
}
}  // namespace FeltElements