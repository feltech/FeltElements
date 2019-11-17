//
// Created by dave on 04/11/2019.
//

#include "MeshFile.hpp"
#include <tetgen/tetgen.h>

using namespace FeltElements;

MeshFile::MeshFile(std::string file_name) : m_file_name{std::move(file_name)}
{
	m_io = std::make_unique<tetgenio>();
	m_io->load_node(const_cast<char *>(m_file_name.c_str()));
	m_io->load_tet(const_cast<char *>(m_file_name.c_str()));
	m_io->load_face(const_cast<char *>(m_file_name.c_str()));
	m_io->numberofpointattributes = 3;
	m_io->pointattributelist = new REAL[m_io->numberofpointattributes * m_io->numberofpoints];
	memset(
		m_io->pointattributelist, 0,
		m_io->numberofpointattributes * m_io->numberofpoints * sizeof(REAL));
}
MeshFile::~MeshFile() = default;

std::size_t MeshFile::num_simplexes() const
{
	return m_io->numberoftetrahedra;
}

std::size_t MeshFile::num_corners() const
{
	return m_io->numberofcorners;
}

std::size_t MeshFile::num_trifaces() const
{
	return m_io->numberoftrifaces;
}

double * MeshFile::points() const
{
	return m_io->pointlist;
}

double * MeshFile::displacements() const
{
	return m_io->pointattributelist;
}

int * MeshFile::corners() const
{
	return m_io->tetrahedronlist;
}

double * MeshFile::tet_corner_point(std::size_t tet_idx, std::size_t corner_idx) const
{
	return &points()[3 * corners()[tet_idx * num_corners() + corner_idx]];
}

double * MeshFile::tet_corner_displacement(std::size_t tet_idx, std::size_t corner_idx) const
{
	return &displacements()[3 * corners()[tet_idx * num_corners() + corner_idx]];
}
