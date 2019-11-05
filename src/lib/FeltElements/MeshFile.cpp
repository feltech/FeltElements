//
// Created by dave on 04/11/2019.
//

#include "MeshFile.hpp"
#include <tetgen/tetgen.h>

using namespace FeltElements;

MeshFile::MeshFile(std::string file_name)
: m_file_name{std::move(file_name)}
{
	m_io = std::make_unique<tetgenio>();
	m_io->load_node(const_cast<char*>(m_file_name.c_str()));
	m_io->load_tet(const_cast<char*>(m_file_name.c_str()));
}
MeshFile::~MeshFile() = default;


std::size_t MeshFile::num_simplexes() const
{
	return m_io->numberoftetrahedra;
}

