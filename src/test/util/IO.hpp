#pragma once
#include <FeltElements/Attributes.hpp>
#include <FeltElements/Typedefs.hpp>
#include <string_view>
#include <tuple>

namespace FeltElements::Test
{
Mesh load_ovm_mesh(std::string_view const& file_name);

std::tuple<Node::Positions, Node::Positions> load_tet(std::string_view const& file_name);

void write_ovm_mesh(
	OpenVolumeMesh::GeometricTetrahedralMeshV3d const& mesh_src,
	Attribute::Vertex::SpatialPosition const& attrib_x,
	std::string_view const& file_name);
}  // namespace FeltElements::Test