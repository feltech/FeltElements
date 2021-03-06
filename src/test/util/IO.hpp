#pragma once
#include <string_view>
#include <tuple>

#include <FeltElements/Attributes.hpp>
#include <FeltElements/Typedefs.hpp>

namespace FeltElements::Test
{
Mesh load_ovm_mesh(std::string_view const & file_name);

std::tuple<Element::NodePositions, Element::NodePositions> load_tet(
	std::string_view const & file_name);

void write_ovm_mesh(
	OpenVolumeMesh::GeometricTetrahedralMeshV3d const & mesh_src,
	Attribute::Vertex::SpatialPosition const & attrib_x,
	std::string_view const & file_name);
}  // namespace FeltElements::Test