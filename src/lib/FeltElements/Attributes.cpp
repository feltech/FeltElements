#include "Attributes.hpp"

namespace FeltElements::Attribute
{

Node::SpatialPosition::SpatialPosition(Node::SpatialPosition::Mesh & mesh) : ThisBase{mesh}
{
	for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
	{
		using OvmVtx = std::decay<decltype(mesh.vertex(*itvtxh))>::type;
		using OvmVtxTensor = Eigen::Tensor<OvmVtx::value_type const, 1>;
		auto const & vtx = mesh.vertex(*itvtxh);
		(*this)[*itvtxh] = Eigen::TensorMap<OvmVtxTensor>{vtx.data(), OvmVtx::size()};
	}
}

Node::Force::Force(Mesh & mesh) : ThisBase{mesh}
{
	for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
		(*this)[*itvtxh].setZero();
}

Element::VertexHandles::VertexHandles(Mesh & mesh) : ThisBase{mesh}
{
	for (auto itcellh = mesh.cells_begin(); itcellh != mesh.cells_end(); itcellh++)
	{
		auto const & vtxhs = mesh.get_cell_vertices(*itcellh);
		std::copy_n(vtxhs.begin(), vtxhs.size(), (*this)[*itcellh].begin());
	}
}

Element::MaterialShapeDerivative::MaterialShapeDerivative(Mesh & mesh, VertexHandles const & vtxhs)
	: ThisBase{mesh}
{
	for (auto itcellh = mesh.cells_begin(); itcellh != mesh.cells_end(); itcellh++)
		(*this)[*itcellh] = Tetrahedron::dN_by_dX(Tetrahedron::X(mesh, vtxhs[*itcellh]));
}
} // namespace FeltElements::Attribute