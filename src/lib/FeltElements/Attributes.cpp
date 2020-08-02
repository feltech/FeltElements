#include "Attributes.hpp"

#include <boost/range.hpp>

namespace FeltElements
{
namespace Attribute
{
namespace Body
{
Properties::Properties(Mesh &mesh) : ThisBase(mesh)
{
	(*(*this)) = MaterialProperties{0, 0, 0};
}
Force::Force(Mesh & mesh) : ThisBase(mesh)
{
	(*(*this)).zeros();
}
Surface::Surface(Mesh & mesh) : ThisBase(mesh)
{
	for (auto halffaceh : boost::make_iterator_range(mesh.halffaces()))
		if (mesh.is_boundary(halffaceh))
			(*this)->push_back(halffaceh);
}
}  // namespace Body

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
MaterialVolume::MaterialVolume(
	Mesh & mesh, Cell::VertexHandles const & vtxh, Vertex::MaterialPosition const & X)
	: ThisBase{mesh}
{
	for (auto cellh : boost::make_iterator_range(mesh.cells()))
		(*this)[cellh] = Derivatives::V(X.for_element(vtxh[cellh]));
}

SpatialVolume::SpatialVolume(Mesh & mesh) : ThisBase{mesh}
{
	for (auto cellh : boost::make_iterator_range(mesh.cells())) (*this)[cellh] = 0;
}

VertexHandles::VertexHandles(Mesh & mesh) : ThisBase{mesh}
{
	for (auto cellh : boost::make_iterator_range(mesh.cells()))
	{
		auto const & vtxhs = mesh.get_cell_vertices(cellh);
		std::copy_n(vtxhs.begin(), vtxhs.size(), (*this)[cellh].begin());
	}
}

MaterialShapeDerivative::MaterialShapeDerivative(
	Mesh & mesh, VertexHandles const & vtxhs, Vertex::MaterialPosition const & X)
	: ThisBase{mesh}
{
	for (auto itcellh = mesh.cells_begin(); itcellh != mesh.cells_end(); itcellh++)
	{ (*this)[*itcellh] = Derivatives::dN_by_dX(X.for_element(vtxhs[*itcellh])); }
}

NodalForces::NodalForces(Mesh & mesh) : ThisBase{mesh}
{
	for (auto itcellh = mesh.cells_begin(); itcellh != mesh.cells_begin(); itcellh++)
		(*this)[*itcellh] = 0;
}

Stiffness::Stiffness(Mesh & mesh) : ThisBase(mesh)
{
	for (auto itcellh = mesh.cells_begin(); itcellh != mesh.cells_begin(); itcellh++)
		(*this)[*itcellh] = 0;
}
}  // namespace Cell
}  // namespace Attribute

Attributes::Attributes(Mesh & mesh)
	: material{mesh},
	  f{mesh},
	  x{mesh},
	  X{mesh},
	  fixed_dof(mesh),
	  vtxh{mesh},
	  V{mesh, vtxh, X},
	  v{mesh},
	  dN_by_dX{mesh, vtxh, X},
	  T{mesh},
	  K{mesh}
{
}
}  // namespace FeltElements