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
		using Pos = Eigen::Vector3d;
		using PosMap = Eigen::Map<Pos>;
		using List = std::array<PosMap, 4>;
		using Index = Eigen::Index;
	};
	using GradientMatrix = Eigen::Matrix<Node::Coord, 3, 3>;
	using ShapeDerivativeMatrix = Eigen::Matrix<Node::Coord, 4, 3>;
	using IsoCoordDerivativeMatrix = Eigen::Matrix<Node::Coord, 3, 4>;

  public:
	Tetrahedron(MeshFile const & mesh, std::size_t tet_idx);

	[[nodiscard]] IsoCoordDerivativeMatrix dx_by_dN() const;
	[[nodiscard]] Node::Pos x(Node::Index const idx) const;
	[[nodiscard]] Node::PosMap const & X(Node::Index const idx) const;
	[[nodiscard]] Node::PosMap const & u(Node::Index const idx) const;
	[[nodiscard]] Node::PosMap & u(Node::Index const idx);

	[[nodiscard]] Tetrahedron::GradientMatrix
	dx_by_dX(ShapeDerivativeMatrix const & dN_by_dX) const;
	[[nodiscard]] Tetrahedron::ShapeDerivativeMatrix dN_by_dX();
	[[nodiscard]] Tetrahedron::ShapeDerivativeMatrix
	dN_by_dx(ShapeDerivativeMatrix const & dN_by_dX, GradientMatrix const & F);

	static IsoCoordDerivativeMatrix const dL_by_dN;
	static ShapeDerivativeMatrix const dN_by_dL;

  private:
	Node::List m_vertices;
	Node::List m_displacements;
};
} // namespace FeltElements
#endif // FELTELEMENTS_TETRAHEDRON_HPP
