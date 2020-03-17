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
		using SpatialCoordProp = OpenVolumeMesh::VertexPropertyT<OpenVolumeMesh::Vec3d>;
	};
	template <Eigen::Index rows = 3, Eigen::Index cols = rows, int options = 0>
	using Matrix = Eigen::Matrix<Scalar, rows, cols, options>;
	using GradientMatrix = Matrix<3, 3>;
	using StressMatrix = GradientMatrix;
	using StiffnessMatrix = GradientMatrix;
	using ShapeDerivativeMatrix = Matrix<4, 3>;
	using IsoCoordDerivativeMatrix = Matrix<3, 4>;
	template <Eigen::Index rows = 3, Eigen::Index cols = 3>
	using MatrixTensor = Eigen::TensorFixedSize<Scalar, Eigen::Sizes<rows, cols>>;
	template <Eigen::Index dim = 3>
	using VectorTensor = Eigen::TensorFixedSize<Scalar, Eigen::Sizes<dim>>;
	using IsoCoordDerivativeTensor = MatrixTensor<3, 4>;
	using ShapeDerivativeTensor = MatrixTensor<4, 3>;
	using ElasticityTensor = Eigen::TensorFixedSize<Scalar, Eigen::Sizes<3, 3, 3, 3>>;
	using ShapeCartesianTransform = MatrixTensor<4, 4>;
	using CartesianDerivativeTensor = MatrixTensor<3, 4>;
	using GradientTensor = MatrixTensor<3, 3>;
	using IndexPair = Eigen::IndexPair<Eigen::Index>;
	template <Eigen::Index num_pairs>
	using IndexPairs = Eigen::array<IndexPair, num_pairs>;

	static IsoCoordDerivativeTensor const dL_by_dN;
	static ShapeDerivativeTensor const dN_by_dL;



	[[nodiscard]] static ElasticityTensor
	neo_hookian_elasticity(Scalar J, Scalar lambda, Scalar mu);
	[[nodiscard]] static StressMatrix
	neo_hookian_stress(Scalar J, GradientMatrix const & b, Scalar lambda, Scalar mu);

	[[nodiscard]] static Scalar J(GradientTensor const & F);
	[[nodiscard]] static GradientTensor b(GradientTensor const & F);

	[[nodiscard]] static GradientTensor dx_by_dX(
		Node::Positions const & x, ShapeDerivativeTensor const & dN_by_dX);
	[[nodiscard]] static GradientTensor dx_by_dX(
		GradientTensor const & dx_by_dL, GradientTensor const & dL_by_dX);
	[[nodiscard]] static GradientMatrix dx_by_dX(
		IsoCoordDerivativeMatrix const & dx_by_dN, ShapeDerivativeMatrix const & dN_by_dX);
	[[nodiscard]] static ShapeDerivativeMatrix dN_by_dx(
		ShapeDerivativeMatrix const & dN_by_dX, GradientMatrix const & F);

	[[nodiscard]] static Scalar V(Node::Positions const & x);

	[[nodiscard]] static GradientTensor dX_by_dL(Node::Positions const & X);
	[[nodiscard]] static GradientTensor dL_by_dX(GradientTensor const & dX_by_dL);
	[[nodiscard]] static CartesianDerivativeTensor dx_by_dN(ShapeCartesianTransform const & N_to_x);
	[[nodiscard]] static ShapeDerivativeTensor dN_by_dX(
		Tetrahedron::GradientTensor const & dL_by_dX);

	[[nodiscard]] static ShapeDerivativeTensor dN_by_dX(ShapeCartesianTransform const & N_to_x);
	[[nodiscard]] static ShapeDerivativeTensor dN_by_dX(Node::Positions const & X);

	[[nodiscard]] static ShapeCartesianTransform N_to_x(Node::Positions const & X);

	[[nodiscard]] static Node::Positions X(Mesh const & mesh, CellHandle const & cellh);
	[[nodiscard]] static Node::SpatialCoordProp x(Mesh & mesh);
	[[nodiscard]] static Node::Positions x(
		Mesh const & mesh, CellHandle const & cellh, Node::SpatialCoordProp const & x_prop);


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
};
} // namespace FeltElements
#endif // FELTELEMENTS_TETRAHEDRON_HPP
