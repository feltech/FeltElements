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

Element::ShapeDerivativeDeterminant const MaterialShapeDerivative::det_dN_by_dL = ([](){
	using Tensor::Index;
	using Tensor::Func::all;
	using Tensor::Func::det;
	Element::ShapeDerivativeDeterminant det_dN_by_dL_;
	for (Index i = 0; i < Node::count; i++)
		for (Index j = 0; j < Node::count; j++)
			for (Index k = 0; k < Node::count; k++)
			{
				Tensor::Matrix<3> dN_by_dL_ijk;
				dN_by_dL_ijk(0, all) = dN_by_dL(i, all);
				dN_by_dL_ijk(1, all) = dN_by_dL(j, all);
				dN_by_dL_ijk(2, all) = dN_by_dL(k, all);
				det_dN_by_dL_(i, j, k) = det(dN_by_dL_ijk);
			}
	return det_dN_by_dL_;
}());

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