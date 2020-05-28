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
		using OvmVtxTensor = Eigen::Tensor<OvmVtx::value_type const, 1>;
		OvmVtx const& vtx = mesh.vertex(*itvtxh);
		(*this)[*itvtxh] = Eigen::TensorMap<OvmVtxTensor>{vtx.data(), OvmVtx::size()};
	}
}

Force::Force(Mesh& mesh) : ThisBase{mesh}
{
	for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
		(*this)[*itvtxh].setZero();
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
	internal::to_tensor(Derivatives::IsoCoordDerivativeMatrix{(Eigen::Matrix4d{} <<
		// (1, L) = A * N
		1, 1, 1, 1,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1).finished().block<3, 4>(1, 0)});

ShapeDerivative const
MaterialShapeDerivative::dN_by_dL = // NOLINT(cert-err58-cpp)
	internal::to_tensor(Derivatives::ShapeDerivativeMatrix{(Eigen::Matrix4d{} <<
		1, 1, 1, 1,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1).finished().inverse().block<4, 3>(0, 1)});
// clang-format on

} // namespace FeltElements::Attribute
}