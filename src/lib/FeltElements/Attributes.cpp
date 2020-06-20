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

Node::Positions SpatialPosition::for_element(Element::Vtxhs const & vtxhs) const
{
	using Tensor::Func::all;
	Node::Positions x;
	for (Tensor::Index node_idx = 0; node_idx < Node::Positions::dimension(0); node_idx++)
		x(node_idx, all) = m_prop[vtxhs[node_idx]];

	return x;
}
} // namespace Node::Attribute

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

NodalForces::NodalForces(Mesh &mesh) : ThisBase{mesh}
{
	for (auto itcellh = mesh.cells_begin(); itcellh != mesh.cells_begin(); itcellh++)
		(*this)[*itcellh] = 0;
}

Stiffness::Stiffness(Mesh& mesh) : ThisBase(mesh)
{
	for (auto itcellh = mesh.cells_begin(); itcellh != mesh.cells_begin(); itcellh++)
		(*this)[*itcellh] = 0;
}
} // namespace FeltElements::Attribute
}