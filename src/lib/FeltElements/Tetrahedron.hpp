//
// Created by dave on 07/11/2019.
//

#ifndef FELTELEMENTS_TETRAHEDRON_HPP
#define FELTELEMENTS_TETRAHEDRON_HPP

#include <cstddef>

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

namespace FeltElements
{
class MeshFile;

template <typename T>
double delta(T const i, T const j) {
	return i == j;
};

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
	using ElasticityTensor = Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3, 3>>;

	Tetrahedron(MeshFile const & mesh, std::size_t tet_idx);

	[[nodiscard]] IsoCoordDerivativeMatrix dx_by_dN() const;
	[[nodiscard]] Node::Pos x(Node::Index idx) const;
	[[nodiscard]] Node::PosMap const & X(Node::Index idx) const;
	[[nodiscard]] Node::PosMap const & u(Node::Index idx) const;
	[[nodiscard]] Node::PosMap & u(Node::Index idx);
	[[nodiscard]] double V() const;

	[[nodiscard]] Tetrahedron::ShapeDerivativeMatrix
	dN_by_dx(ShapeDerivativeMatrix const & dN_by_dX);
	[[nodiscard]] Tetrahedron::GradientMatrix
	dx_by_dX(ShapeDerivativeMatrix const & dN_by_dX) const;
	[[nodiscard]] Tetrahedron::ShapeDerivativeMatrix dN_by_dX();

	static IsoCoordDerivativeMatrix const dL_by_dN;
	static ShapeDerivativeMatrix const dN_by_dL;
	ElasticityTensor
	neo_hookian_elasticity(GradientMatrix const & F, double const lambda, double const mu);

  private:
	Node::List const m_vertices;
	Node::List m_displacements;
	double const m_volume;
};
} // namespace FeltElements
#endif // FELTELEMENTS_TETRAHEDRON_HPP
