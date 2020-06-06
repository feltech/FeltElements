#include "Attributes.hpp"
#include "internal/Conversions.hpp"

namespace FeltElements
{

namespace Node::Attribute
{
SpatialPosition::SpatialPosition(Mesh& mesh) : ThisBase{mesh}
{
	for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
	{
		using OvmVtx = std::decay<decltype(mesh.vertex(*itvtxh))>::type;
		OvmVtx vtx = mesh.vertex(*itvtxh);
		Data & position = (*this)[*itvtxh];
		position = Tensor::BaseMap<OvmVtx::value_type, OvmVtx::size()>{vtx.data()};
	}
}

Force::Force(Mesh& mesh) : ThisBase{mesh}
{
	for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
		(*this)[*itvtxh] = 0;
}
}

namespace Element::Attribute
{
VertexHandles::VertexHandles(Mesh& mesh) : ThisBase{mesh}
{
	for (auto itcellh = mesh.cells_begin(); itcellh != mesh.cells_end(); itcellh++)
	{
		auto const& vtxhs = mesh.get_cell_vertices(*itcellh);
		std::copy_n(vtxhs.begin(), vtxhs.size(), (*this)[*itcellh].begin());
	}
}

MaterialShapeDerivative::MaterialShapeDerivative(Mesh& mesh, VertexHandles const& vtxhs)
	: ThisBase{mesh}
{
	for (auto itcellh = mesh.cells_begin(); itcellh != mesh.cells_end(); itcellh++)
		(*this)[*itcellh] = Derivatives::dN_by_dX(Derivatives::X(mesh, vtxhs[*itcellh]));
}

// clang-format off
IsoCoordDerivative const
MaterialShapeDerivative::dL_by_dN = // NOLINT(cert-err58-cpp)
	Tensor::Matrix<4, 4>{
		{1, 1, 1, 1},
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1}}(Fastor::fseq<1, 4>(), Fastor::all);

ShapeDerivative const
MaterialShapeDerivative::dN_by_dL = // NOLINT(cert-err58-cpp)
	Fastor::evaluate(Fastor::inv(Tensor::Matrix<4, 4>{
		{1, 1, 1, 1},
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1}}))(Fastor::all, Fastor::fseq<1, 4>());
// clang-format on

} // namespace FeltElements::Attribute
}