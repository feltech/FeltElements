#pragma once
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/iterator_range.hpp>

#include "Attributes.hpp"
#include "Typedefs.hpp"

namespace FeltElements
{
struct MeshDecorator
{
	using CellRng = typename boost::iterator_range<typename boost::range_iterator<
		const std::pair<OpenVolumeMesh::CellIter, OpenVolumeMesh::CellIter>>::type>;

	using VtxRng = typename boost::iterator_range<typename boost::range_iterator<
		const std::pair<OpenVolumeMesh::VertexIter, OpenVolumeMesh::VertexIter>>::type>;

	Mesh const & mesh;
	Attributes const & attrs;

	VtxRng vertices{boost::make_iterator_range(mesh.vertices())};

	CellRng cells{boost::make_iterator_range(mesh.cells())};

	[[nodiscard]] inline auto Xs() const
	{
		return cells | boost::adaptors::transformed([this](auto const & cellh) {
				   return attrs.X.for_element(attrs.vtxhs[cellh]);
			   });
	}

	[[nodiscard]] inline auto xs() const
	{
		return cells | boost::adaptors::transformed([this](auto const & cellh) {
				   return attrs.x.for_element(attrs.vtxhs[cellh]);
			   });
	}

	static Mesh fromFile(std::string const & file_path);

	void toFile(std::string const & file_path);

};
}  // namespace FeltElements