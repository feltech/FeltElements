//
// Created by dave on 07/11/2019.
//

#ifndef FELTELEMENTS_TETRAHEDRON_HPP
#define FELTELEMENTS_TETRAHEDRON_HPP

#include <cstddef>

#include <Eigen/Dense>

namespace FeltElements
{
class MeshFile;

class Tetrahedron
{
	struct Node
	{
		using Coord = double;
		using Pos = Eigen::Map<Eigen::Vector3d>;
		using List = std::array<Pos, 4>;
		using Index = Eigen::Index;
	};
	using GradientMatrix = Eigen::Matrix<Node::Coord, 3, 3>;
public:
	Tetrahedron(MeshFile const& mesh, std::size_t tet_idx);

	[[nodiscard]] Node::Pos const& vertex(Node::Index idx) const;
	[[nodiscard]] Node::Pos & displacement(Node::Index idx);

	[[nodiscard]] GradientMatrix deformation_gradient() const;

private:
	Node::List m_vertices;
	Node::List m_displacements;
	Eigen::Matrix<Node::Coord, 4, 3> m_dN_by_dX;

};
}
#endif // FELTELEMENTS_TETRAHEDRON_HPP
