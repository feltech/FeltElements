#include "Attributes.hpp"

#include <boost/range/algorithm/find.hpp>

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
		BoundaryElement::Vtxhs vtxhs;
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

Boundary::Boundary(Mesh & mesh, VertexHandles const & vtxhs) : CellBase(mesh)
{
	for (auto cellh : boost::make_iterator_range(mesh.cells()))
	{
		auto const& cell_vtxhs = vtxhs[cellh];

		for (auto hfh : boost::make_iterator_range(mesh.cell_halffaces(cellh)))
		{
			if (!mesh.is_boundary(Mesh::face_handle(hfh)))
				continue;

			BoundaryElement::VtxhIdxs vtxhidxs;
			Tensor::Index vtx_idx = 0;
			// Get opposite half-face so winding order is on the outside (i.e. normal calculated
			// from cross product points outward).
			Mesh::Face const boundary = mesh.opposite_halfface(hfh);

			for (auto const & heh : boundary.halfedges())
			{
				Vtxh const & vtxh = mesh.halfedge(heh).from_vertex();
				// Get index of handle in cell.
				auto const vtxhidx = static_cast<Tensor::Index>(
					std::distance(cell_vtxhs.begin(), boost::range::find(cell_vtxhs, vtxh)));
				vtxhidxs[vtx_idx++] = vtxhidx;
			}
			assert(vtx_idx == 3);

			(*this)[cellh].push_back(vtxhidxs);
		}
	}
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
	  boundary{mesh, vtxh},
	  dN_by_dX{mesh, vtxh, X},
	  R{mesh},
	  K{mesh}
{
}
}  // namespace FeltElements