#include "Attributes.hpp"

namespace FeltElements
{
namespace Attribute
{
namespace Vertex
{
MaterialPosition::MaterialPosition(Mesh& mesh) : ThisBase(mesh) {}
SpatialPosition::SpatialPosition(Mesh& mesh) : ThisBase(mesh) {}

FixedDOF::FixedDOF(Mesh& mesh) : ThisBase(mesh)
{
	for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
		(*this)[*itvtxh] = 0;
}
}  // namespace Vertex

namespace Cell
{
VertexHandles::VertexHandles(Mesh& mesh) : ThisBase{mesh}
{
	for (auto itcellh = mesh.cells_begin(); itcellh != mesh.cells_end(); itcellh++)
	{
		auto const& vtxhs = mesh.get_cell_vertices(*itcellh);
		std::copy_n(vtxhs.begin(), vtxhs.size(), (*this)[*itcellh].begin());
	}
}

MaterialShapeDerivative::MaterialShapeDerivative(
	Mesh& mesh, VertexHandles const& vtxhs, Vertex::MaterialPosition const& X)
	: ThisBase{mesh}
{
	for (auto itcellh = mesh.cells_begin(); itcellh != mesh.cells_end(); itcellh++)
	{ (*this)[*itcellh] = Derivatives::dN_by_dX(X.for_element(vtxhs[*itcellh])); }
}

// clang-format off
Element::IsoCoordDerivative const MaterialShapeDerivative::dL_by_dN = // NOLINT(cert-err58-cpp)
	Tensor::Matrix<4, 4>{
		{1, 1, 1, 1},
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1}}(Fastor::fseq<1, 4>(), Fastor::all);

Element::ShapeDerivative const MaterialShapeDerivative::dN_by_dL = // NOLINT(cert-err58-cpp)
	Fastor::evaluate(Fastor::inv(Tensor::Matrix<4, 4>{
		{1, 1, 1, 1},
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1}}))(Fastor::all, Fastor::fseq<1, 4>());
// clang-format on

NodalForces::NodalForces(Mesh& mesh) : ThisBase{mesh}
{
	for (auto itcellh = mesh.cells_begin(); itcellh != mesh.cells_begin(); itcellh++)
		(*this)[*itcellh] = 0;
}

Stiffness::Stiffness(Mesh& mesh) : ThisBase(mesh)
{
	for (auto itcellh = mesh.cells_begin(); itcellh != mesh.cells_begin(); itcellh++)
		(*this)[*itcellh] = 0;
}
}  // namespace Cell
}  // namespace Attribute

Attributes::Attributes(Mesh& mesh)
	: x{mesh}, X{mesh}, fixed_dof(mesh),
	  vtxh{mesh}, dN_by_dX{mesh, vtxh, X}, T{mesh}, K{mesh}
{
}
}  // namespace FeltElements