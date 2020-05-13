//
// Created by dave on 04/11/2019.
//

#include "TetGenIO.hpp"

#include <boost/range/irange.hpp>

using namespace FeltElements;

TetGenIO::TetGenIO(std::string&& file_name) : m_file_name{file_name}, m_io{}
{
	m_io.load_node(const_cast<char *>(m_file_name.c_str()));
	m_io.load_tet(const_cast<char *>(m_file_name.c_str()));
	m_io.load_face(const_cast<char *>(m_file_name.c_str()));
	m_io.numberofpointattributes = 3;
	m_io.pointattributelist = new REAL[m_io.numberofpointattributes * m_io.numberofpoints];
	memset(
		m_io.pointattributelist, 0,
		m_io.numberofpointattributes * m_io.numberofpoints * sizeof(REAL));
}
TetGenIO::~TetGenIO() = default;


OpenVolumeMesh::GeometricTetrahedralMeshV3d TetGenIO::to_mesh() const
{
	using VtxhList = std::vector<OpenVolumeMesh::VertexHandle>;

	VtxhList hvtxs(num_points());
	VtxhList tet_hvtxs(num_corners());
	OpenVolumeMesh::GeometricTetrahedralMeshV3d mesh;

	for (std::size_t vtx_idx : boost::irange(num_points()))
	{
		auto const [x, y, z] = vertex(vtx_idx);
		hvtxs[vtx_idx] = mesh.add_vertex(OpenVolumeMesh::Vec3d{x, y, z});
	}

	for (std::size_t tet_idx : boost::irange(num_simplexes()))
	{
		for (std::size_t corner_idx : boost::irange(num_corners()))
		{
			std::size_t const vtx_idx = tet_vertex_idx(tet_idx, corner_idx);
			tet_hvtxs[corner_idx] = hvtxs[vtx_idx];
		}
		mesh.add_cell(tet_hvtxs);
	}
	return mesh;
}

std::size_t TetGenIO::num_points() const
{
	return m_io.numberofpoints;
}

std::size_t TetGenIO::num_simplexes() const
{
	return m_io.numberoftetrahedra;
}

std::size_t TetGenIO::num_corners() const
{
	return m_io.numberofcorners;
}

double * TetGenIO::points() const
{
	return m_io.pointlist;
}

double * TetGenIO::displacements() const
{
	return m_io.pointattributelist;
}

int * TetGenIO::corners() const
{
	return m_io.tetrahedronlist;
}

std::array<double, 3> TetGenIO::vertex(std::size_t const vertex_idx) const
{
	std::size_t const idx = point_idx(vertex_idx);
	return std::array<double, 3>{points()[idx], points()[idx + 1], points()[idx + 2]};
}

std::size_t TetGenIO::tet_vertex_idx(std::size_t const tet_idx, std::size_t const corner_idx) const
{
	return corners()[tet_idx * num_corners() + corner_idx];
}

std::size_t TetGenIO::point_idx(std::size_t vertex_idx)
{
	return 3 * vertex_idx;
}
