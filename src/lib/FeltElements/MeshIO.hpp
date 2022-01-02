#pragma once
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/iterator_range.hpp>

#include "Attributes.hpp"
#include "Typedefs.hpp"

namespace FeltElements
{

struct MeshIO
{
	using CellRng = typename boost::iterator_range<typename boost::range_iterator<
		const std::pair<OpenVolumeMesh::CellIter, OpenVolumeMesh::CellIter>>::type>;

	using VtxRng = typename boost::iterator_range<typename boost::range_iterator<
		const std::pair<OpenVolumeMesh::VertexIter, OpenVolumeMesh::VertexIter>>::type>;

	using VtxCellRng = typename boost::iterator_range<typename boost::range_iterator<
		const std::pair<OpenVolumeMesh::VertexCellIter, OpenVolumeMesh::VertexCellIter>>::type>;

	using HalfEdgeRng = typename boost::iterator_range<typename boost::range_iterator<
		const std::pair<OpenVolumeMesh::HalfEdgeIter, OpenVolumeMesh::HalfEdgeIter>>::type>;

	using HalfEdgeCellRng = typename boost::iterator_range<typename boost::range_iterator<
		const std::pair<OpenVolumeMesh::HalfEdgeCellIter, OpenVolumeMesh::HalfEdgeCellIter>>::type>;

	using VtxHalfEdgeRng = typename boost::iterator_range<typename boost::range_iterator<
		const std::pair<OpenVolumeMesh::VertexOHalfEdgeIter, OpenVolumeMesh::VertexOHalfEdgeIter>>::
															  type>;

	using Heh = OpenVolumeMesh::HalfEdgeHandle;

	[[nodiscard]] static Mesh fromFile(std::string const & file_path, Scalar scale = 1);

	void toFile(std::string const & file_path) const;

	[[nodiscard]] inline auto Xs() const
	{
		return cells() |
			boost::adaptors::transformed([&attrs = attrs](auto const & cellh)
										 { return attrs.X.for_element(attrs.vtxhs[cellh]); });
	}

	[[nodiscard]] inline auto xs() const
	{
		return cells() |
			boost::adaptors::transformed([&attrs = attrs](auto const & cellh)
										 { return attrs.x.for_element(attrs.vtxhs[cellh]); });
	}

	[[nodiscard]] inline VtxRng vertices() const
	{
		return boost::make_iterator_range(mesh.vertices());
	}

	[[nodiscard]] inline CellRng cells() const
	{
		return boost::make_iterator_range(mesh.cells());
	}

	[[nodiscard]] inline VtxCellRng vertex_cells(Vtxh const & vtxh) const
	{
		return boost::make_iterator_range(mesh.vertex_cells(vtxh));
	}

	[[nodiscard]] inline HalfEdgeRng halfedges() const
	{
		return boost::make_iterator_range(mesh.halfedges());
	}

	[[nodiscard]] inline HalfEdgeCellRng halfedge_cells(Heh const & heh) const
	{
		return boost::make_iterator_range(mesh.halfedge_cells(heh));
	}

	[[nodiscard]] inline VtxHalfEdgeRng outgoing_halfedges(Vtxh const & vtxh) const
	{
		return boost::make_iterator_range(mesh.outgoing_halfedges(vtxh));
	}

	Mesh const & mesh;
	Attributes const & attrs;
};
}  // namespace FeltElements