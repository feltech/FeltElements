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
	using Mesh = OpenVolumeMesh::GeometricTetrahedralMeshV3d;
	using CellHandle = OpenVolumeMesh::CellHandle;

	struct Node
	{
		using Coord = Scalar;
		using Pos = Eigen::Vector3d;
		using PosMap = Eigen::Map<Pos>;
		using List = std::array<PosMap, 4>;
		using Index = Eigen::Index;
		using Positions = Eigen::TensorFixedSize<Scalar, Eigen::Sizes<4, 3>>;
		using PosProperty = OpenVolumeMesh::VertexPropertyT<Pos>;
	};
	using ShapeDerivativeTensor = Eigen::TensorFixedSize<Node::Coord, Eigen::Sizes<4, 3>>;
	using ShapeDerivativeMatrix = Eigen::Matrix<Node::Coord, 4, 3, Eigen::RowMajor>;
	using IsoCoordDerivativeMatrix = Eigen::Matrix<Node::Coord, 3, 4>;
	using GradientMatrix = Eigen::Matrix<Node::Coord, 3, 3>;
	using StressMatrix = GradientMatrix;
	using StiffnessMatrix = GradientMatrix;
	using ElasticityTensor = Eigen::TensorFixedSize<Scalar, Eigen::Sizes<3, 3, 3, 3>>;
	template <Eigen::Index... dim>
	using MatrixTensor = Eigen::TensorFixedSize<Scalar, Eigen::Sizes<dim...>>;
	template <Eigen::Index dim = 3>
	using VectorTensor = Eigen::TensorFixedSize<Scalar, Eigen::Sizes<dim>>;

	Tetrahedron() = default;
	Tetrahedron(TetGenIO const & mesh, std::size_t tet_idx);

	static IsoCoordDerivativeMatrix const dL_by_dN;
	static ShapeDerivativeMatrix const dN_by_dL;
	[[nodiscard]] static ShapeDerivativeTensor dN_by_dX(Node::Positions const & X);
	[[nodiscard]] ShapeDerivativeMatrix dN_by_dX() const;
	[[nodiscard]] IsoCoordDerivativeMatrix dx_by_dN() const;
	[[nodiscard]] static Node::Positions x(Node::Positions const & X, Node::Positions const & u);
	[[nodiscard]] Node::Pos x(Node::Index idx) const;
	[[nodiscard]] static Node::Positions X(Mesh const & mesh, CellHandle const & cellh);
	[[nodiscard]] Node::PosMap const & X(Node::Index idx) const;
	[[nodiscard]] static Node::Positions u(
		Mesh const & mesh, Node::PosProperty const & displacements, CellHandle const & cellh);
	[[nodiscard]] Node::PosMap const & u(Node::Index idx) const;
	[[nodiscard]] Node::PosMap & u(Node::Index idx);
	[[nodiscard]] Scalar V() const;
	[[nodiscard]] Scalar v() const;

	[[nodiscard]] static GradientMatrix dx_by_dX(
		IsoCoordDerivativeMatrix const & dx_by_dN, ShapeDerivativeMatrix const & dN_by_dX);
	[[nodiscard]] static ShapeDerivativeMatrix dN_by_dx(
		ShapeDerivativeMatrix const & dN_by_dX, GradientMatrix const & F);
	[[nodiscard]] static Scalar J(GradientMatrix const & F);
	[[nodiscard]] static GradientMatrix b(GradientMatrix const & F);
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
	static Node::Pos null_pos;
	Node::List m_vertices = {
		Node::PosMap{null_pos.data(), null_pos.size()},
		Node::PosMap{null_pos.data(), null_pos.size()},
		Node::PosMap{null_pos.data(), null_pos.size()},
		Node::PosMap{null_pos.data(), null_pos.size()}};
	Node::List m_displacements = {
		Node::PosMap{null_pos.data(), null_pos.size()},
			Node::PosMap{null_pos.data(), null_pos.size()},
			Node::PosMap{null_pos.data(), null_pos.size()},
			Node::PosMap{null_pos.data(), null_pos.size()}};
	Scalar m_material_volume;
};
} // namespace FeltElements
#endif // FELTELEMENTS_TETRAHEDRON_HPP
