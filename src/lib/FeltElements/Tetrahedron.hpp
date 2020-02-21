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
class TetGenIO;

class Tetrahedron
{
public:
	using Scalar = double;
	struct Node
	{
		using Coord = Scalar;
		using Pos = Eigen::Vector3d;
		using PosMap = Eigen::Map<Pos>;
		using List = std::array<PosMap, 4>;
		using Index = Eigen::Index;
	};
	using ShapeDerivativeMatrix = Eigen::Matrix<Node::Coord, 4, 3, Eigen::RowMajor>;
	using IsoCoordDerivativeMatrix = Eigen::Matrix<Node::Coord, 3, 4>;
	using GradientMatrix = Eigen::Matrix<Node::Coord, 3, 3>;
	using StressMatrix = GradientMatrix;
	using StiffnessMatrix = GradientMatrix;
	using ElasticityTensor = Eigen::TensorFixedSize<Scalar, Eigen::Sizes<3, 3, 3, 3>>;

	Tetrahedron(TetGenIO const & mesh, std::size_t tet_idx);

	static IsoCoordDerivativeMatrix const dL_by_dN;
	static ShapeDerivativeMatrix const dN_by_dL;
	[[nodiscard]] Tetrahedron::ShapeDerivativeMatrix dN_by_dX() const;
	[[nodiscard]] IsoCoordDerivativeMatrix dx_by_dN() const;
	[[nodiscard]] Node::Pos x(Node::Index idx) const;
	[[nodiscard]] Node::PosMap const & X(Node::Index idx) const;
	[[nodiscard]] Node::PosMap const & u(Node::Index idx) const;
	[[nodiscard]] Node::PosMap & u(Node::Index idx);
	[[nodiscard]] Scalar V() const;
	[[nodiscard]] Scalar v() const;

	[[nodiscard]] static Tetrahedron::GradientMatrix dx_by_dX(
		IsoCoordDerivativeMatrix const & dx_by_dN, ShapeDerivativeMatrix const & dN_by_dX);
	[[nodiscard]] static Tetrahedron::ShapeDerivativeMatrix dN_by_dx(
		ShapeDerivativeMatrix const & dN_by_dX, GradientMatrix const & F);
	[[nodiscard]] static Scalar J(GradientMatrix const & F);
	[[nodiscard]] static Tetrahedron::GradientMatrix b(GradientMatrix const & F);
	[[nodiscard]] static ElasticityTensor
	neo_hookian_elasticity(Scalar J, Scalar lambda, Scalar mu);
	[[nodiscard]] static StressMatrix
	neo_hookian_stress(Scalar J, GradientMatrix const & b, Scalar lambda, Scalar mu);

	static StiffnessMatrix Kcab(
		ShapeDerivativeMatrix const & dN_by_dx,
		Scalar v,
		ElasticityTensor const & c,
		Node::Index a,
		Node::Index b);

	static StiffnessMatrix Ksab(
		ShapeDerivativeMatrix const & dN_by_dx,
		Scalar v,
		GradientMatrix const & sigma,
		Node::Index a,
		Node::Index b);

private:
	Node::List const m_vertices;
	Node::List m_displacements;
	Scalar const m_material_volume;
};
} // namespace FeltElements
#endif // FELTELEMENTS_TETRAHEDRON_HPP
