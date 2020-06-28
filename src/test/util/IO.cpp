#include "IO.hpp"
#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <fmt/format.h>

namespace FeltElements::Test
{
Mesh load_ovm_mesh(std::string_view const& file_name)
{
	using namespace FeltElements;
	FeltElements::Mesh mesh;
	OpenVolumeMesh::IO::FileManager file_manager{};
	file_manager.readFile(file_name.data(), mesh);
	return mesh;
}

std::tuple<Node::Positions, Node::Positions> load_tet(std::string_view const& file_name)
{
	using namespace FeltElements;
	FeltElements::Mesh mesh = load_ovm_mesh(file_name);
	FeltElements::Attributes attrib{mesh};
	auto const& cell_vtxhs = attrib.vtxh[FeltElements::Cellh{0}];

	return std::tuple(attrib.X.for_element(cell_vtxhs), attrib.x.for_element(cell_vtxhs));
}

void write_ovm_mesh(
	FeltElements::Mesh const& mesh_src,
	FeltElements::Attribute::Vertex::SpatialPosition const& attrib_x,
	std::string_view const& file_name)
{
	FeltElements::Mesh mesh_dst{mesh_src};
	for (auto const& vtxh : boost::make_iterator_range(mesh_src.vertices()))
	{
		auto const& x = attrib_x[vtxh];
		mesh_dst.set_vertex(vtxh, {x(0), x(1), x(2)});
	}
	OpenVolumeMesh::IO::FileManager{}.writeFile(
		fmt::format("artefacts/{}.ovm", file_name), mesh_dst);
}
}
