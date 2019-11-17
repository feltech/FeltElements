//
// Created by dave on 07/11/2019.
//

#include "Tetrahedron.hpp"
#include "MeshFile.hpp"

namespace FeltElements
{
Tetrahedron::Tetrahedron(FeltElements::MeshFile const & mesh, std::size_t const tet_idx)
	:
	m_vertices{
		  Node::Pos{mesh.tet_corner_point(tet_idx, 0)},
		  Node::Pos{mesh.tet_corner_point(tet_idx, 1)},
		  Node::Pos{mesh.tet_corner_point(tet_idx, 2)},
		  Node::Pos{mesh.tet_corner_point(tet_idx, 3)},
	},
	m_displacements{
		Node::Pos{mesh.tet_corner_displacement(tet_idx, 0)},
		Node::Pos{mesh.tet_corner_displacement(tet_idx, 1)},
		Node::Pos{mesh.tet_corner_displacement(tet_idx, 2)},
		Node::Pos{mesh.tet_corner_displacement(tet_idx, 3)},
	}
{
	Eigen::Matrix4d mat_interp;
	// Interpolation: (1, x, y, z)^T = mat_interp * N, where N is 4x natural coordinates (corners).
	// clang-format off
	mat_interp <<
		1,	1,	1,	1,
		m_vertices[0][0], m_vertices[1][0], m_vertices[2][0], m_vertices[3][0],
		m_vertices[0][1], m_vertices[1][1], m_vertices[2][1], m_vertices[3][1],
		m_vertices[0][2], m_vertices[1][2], m_vertices[2][2], m_vertices[3][2];
	// clang-format on

	// dN_i/dX_j, where i is row, j is column
	m_dN_by_dX = mat_interp.inverse().block<4, 3>(0, 1);
}

Tetrahedron::Node::Pos const & Tetrahedron::vertex(Node::Index idx) const
{
	return m_vertices[idx];
}

Tetrahedron::Node::Pos & Tetrahedron::displacement(Tetrahedron::Node::Index idx)
{
	return m_displacements[idx];
}

Tetrahedron::GradientMatrix Tetrahedron::deformation_gradient() const
{
	Eigen::Matrix<Node::Coord, 3, 4> mat_spatial;
	mat_spatial.col(0) = m_vertices[0] + m_displacements[0];
	mat_spatial.col(1) = m_vertices[1] + m_displacements[1];
	mat_spatial.col(2) = m_vertices[2] + m_displacements[2];
	mat_spatial.col(3) = m_vertices[3] + m_displacements[3];

	GradientMatrix mat_deformation = mat_spatial * m_dN_by_dX;

	return mat_deformation;
}
} // namespace FeltElements
