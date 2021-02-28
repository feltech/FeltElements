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
	for (auto cellh : boost::make_iterator_range(mesh.cells())) (*this)[cellh] = 0;
}

Stiffness::Stiffness(Mesh & mesh) : ThisBase(mesh)
{
	for (auto cellh : boost::make_iterator_range(mesh.cells())) (*this)[cellh] = 0;
}

Boundary::Boundary(Mesh & mesh, VertexHandles const & vtxhs) : CellBase(mesh)
{
	for (auto cellh : boost::make_iterator_range(mesh.cells()))
	{
		auto const & cell_vtxhs = vtxhs[cellh];

		for (auto hfh : boost::make_iterator_range(mesh.cell_halffaces(cellh)))
		{
			if (!mesh.is_boundary(Mesh::face_handle(hfh)))
				continue;

			// Indices of vertex handles for this face within this cell's list of vertex handles.
			BoundaryElement::VtxhIdxs vtxh_idxs;
			// Get opposite half-face so winding order is on the outside (i.e. normal calculated
			// from cross product points outward).
			Mesh::Face const & boundary_face = mesh.opposite_halfface(hfh);

			for (auto const & [vtx_idx, heh] : boost::adaptors::index(boundary_face.halfedges()))
			{
				// Handle to vertex at start endpoint of this halfedge.
				Vtxh const & vtxh = mesh.halfedge(heh).from_vertex();
				using VtxhIdx = decltype(vtxh_idxs)::size_type;
				// Face vertex's index within cell's list of vertex handles.
				vtxh_idxs[static_cast<VtxhIdx>(vtx_idx)] = index_of(cell_vtxhs, vtxh);
			}

			(*this)[cellh].push_back(vtxh_idxs);
		}
	}
}
}  // namespace Cell
}  // namespace Attribute

Attributes::Attributes(Mesh & mesh)
	: material{mesh},
	  forces{mesh},
	  x{mesh},
	  X{mesh},
	  fixed_dof(mesh),
	  vtxhs{mesh},
	  boundary_faces_vtxh_idxs{mesh, vtxhs},
	  dN_by_dX{mesh, vtxhs, X},
	  R{mesh},
	  K{mesh}
{
}
}  // namespace FeltElements