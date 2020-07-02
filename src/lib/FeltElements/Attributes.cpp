#include "Attributes.hpp"

#include <boost/range.hpp>

namespace FeltElements
{
namespace Attribute
{
namespace Global
{
BodyForce::BodyForce(Mesh & mesh) : ThisBase(mesh)
{
	(*(*this)).zeros();
}

BodyForce::Data const & BodyForce::operator*() const
{
	return (*this)[Handle{0}];
}

BodyForce::Data & BodyForce::operator*()
{
	return (*this)[Handle{0}];
}
}  // namespace Global

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
	: f{mesh},
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