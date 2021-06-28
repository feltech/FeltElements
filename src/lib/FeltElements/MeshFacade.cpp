#include "MeshFacade.hpp"

#include <filesystem>

#include <fmt/format.h>
#include <OpenVolumeMesh/FileManager/FileManager.hh>

namespace FeltElements
{
Mesh MeshIO::fromFile(std::string const & file_path, Scalar scale)
{
	Mesh mesh{};
	OpenVolumeMesh::IO::FileManager file_manager{};
	if (!file_manager.readFile(file_path, mesh))
		throw std::filesystem::filesystem_error{
			fmt::format("Unable to read mesh file '{}'", file_path),
			std::make_error_code(std::errc::no_such_file_or_directory)};

	for (auto const & vtxh : boost::make_iterator_range(mesh.vertices()))
		mesh.set_vertex(vtxh, mesh.vertex(vtxh) * scale);

	return mesh;
}

void MeshIO::toFile(std::string const & file_path) const
{
	Mesh out{mesh};
	OpenVolumeMesh::IO::FileManager file_manager{};

	for (auto const & vtxh : MeshIters{mesh, attrs}.vertices)
	{
		auto const & x = attrs.x[vtxh];
		out.set_vertex(vtxh, {x(0), x(1), x(2)});
	}
	file_manager.writeFile(file_path, out);
}
}  // namespace FeltElements
