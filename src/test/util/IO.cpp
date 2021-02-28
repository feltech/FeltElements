#include "IO.hpp"

#include <fmt/format.h>

#include <OpenVolumeMesh/FileManager/FileManager.hh>

#include <FeltElements/MeshFacade.hpp>

namespace FeltElements::Test
{
std::tuple<Element::NodePositions, Element::NodePositions> load_tet(std::string_view const & file_name)
{
	using namespace FeltElements;
	FeltElements::Mesh mesh = MeshIO::fromFile(std::string{file_name});
	FeltElements::Attributes attrib{mesh};
	auto const & cell_vtxhs = attrib.vtxhs[FeltElements::Cellh{0}];

	return std::tuple(attrib.X.for_element(cell_vtxhs), attrib.x.for_element(cell_vtxhs));
}

}  // namespace FeltElements::Test
