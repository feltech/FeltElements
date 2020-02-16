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
	using ShapeDerivativeMatrix = Eigen::Matrix<Node::Coord, 4, 3, Eigen::RowMajor>;
	using IsoCoordDerivativeMatrix = Eigen::Matrix<Node::Coord, 3, 4>;
  public:
	using GradientMatrix = Eigen::Matrix<Node::Coord, 3, 3>;
	using ElasticityTensor = Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3, 3>>;

	Tetrahedron(MeshFile const & mesh, std::size_t tet_idx);

	static IsoCoordDerivativeMatrix const dL_by_dN;
	static ShapeDerivativeMatrix const dN_by_dL;
	[[nodiscard]] Tetrahedron::ShapeDerivativeMatrix dN_by_dX() const;
	[[nodiscard]] IsoCoordDerivativeMatrix dx_by_dN() const;
	[[nodiscard]] Node::Pos x(Node::Index idx) const;
	[[nodiscard]] Node::PosMap const & X(Node::Index idx) const;
	[[nodiscard]] Node::PosMap const & u(Node::Index idx) const;
	[[nodiscard]] Node::PosMap & u(Node::Index idx);
	[[nodiscard]] double V() const;

	[[nodiscard]] static Tetrahedron::GradientMatrix dx_by_dX(
		IsoCoordDerivativeMatrix const & dx_by_dN, ShapeDerivativeMatrix const & dN_by_dX);
	[[nodiscard]] static Tetrahedron::ShapeDerivativeMatrix dN_by_dx(
		ShapeDerivativeMatrix const & dN_by_dX, GradientMatrix const & F);
	[[nodiscard]] static double J(GradientMatrix const & F);
	[[nodiscard]] static Tetrahedron::GradientMatrix b(GradientMatrix const & F);
	[[nodiscard]] static ElasticityTensor
	neo_hookian_elasticity(double J, double lambda, double mu);
	[[nodiscard]] static GradientMatrix
	neo_hookian_stress(double J, GradientMatrix const & b, double lambda, double mu);

	static Tetrahedron::GradientMatrix Kcab(
		Tetrahedron::ShapeDerivativeMatrix const & dN_by_dx,
		Tetrahedron::ElasticityTensor const & c,
		Tetrahedron::Node::Index a,
		Tetrahedron::Node::Index b);

private:
	Node::List const m_vertices;
	Node::List m_displacements;
	double const m_volume;
};
} // namespace FeltElements
#endif // FELTELEMENTS_TETRAHEDRON_HPP
