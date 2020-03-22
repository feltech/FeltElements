//
// Created by dave on 04/11/2019.
//
#pragma once

#include <memory>
#include <string>

#include <tetgen/tetgen.h>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>

namespace FeltElements
{
class TetGenIO
{
  public:
	explicit TetGenIO(std::string&& file_name);
	~TetGenIO();

	[[nodiscard]] OpenVolumeMesh::GeometricTetrahedralMeshV3d to_mesh() const;

	[[nodiscard]] std::size_t num_points() const;
	[[nodiscard]] std::size_t num_simplexes() const;
	[[nodiscard]] std::size_t num_corners() const;
	[[nodiscard]] double * points() const;
	[[nodiscard]] double * displacements() const;
	[[nodiscard]] int * corners() const;

	[[nodiscard]] double * tet_corner_point(std::size_t tet_idx, std::size_t corner_idx) const;
	[[nodiscard]] double *
	tet_corner_displacement(std::size_t tet_idx, std::size_t corner_idx) const;
	[[nodiscard]] std::size_t tet_vertex_idx(std::size_t tet_idx, std::size_t corner_idx) const;
	[[nodiscard]] std::array<double, 3> vertex(std::size_t vertex_idx) const;
	[[nodiscard]] std::size_t point_idx(std::size_t vertex_idx) const;

  private:
	std::string m_file_name;
	tetgenio m_io;
};
} // namespace FeltElements
