#pragma once

#include <Eigen/Dense>
#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <cstddef>

#include "Typedefs.hpp"

namespace FeltElements
{
class Derivatives
{
public:
	using Vtxh = OpenVolumeMesh::VertexHandle;
	using Vtxhs = std::array<Vtxh, Node::count>;
	using PosProperty = OpenVolumeMesh::VertexPropertyT<OpenVolumeMesh::Vec3d>;
	using SpatialCoordProp = OpenVolumeMesh::VertexPropertyT<OpenVolumeMesh::Vec3d>;

	template <Eigen::Index rows = 3, Eigen::Index cols = rows, int options = 0>
	using Matrix = Eigen::Matrix<Scalar, rows, cols, options>;
	using GradientMatrix = Matrix<3, 3>;
	using ShapeDerivativeMatrix = Matrix<4, 3>;
	using IsoCoordDerivativeMatrix = Matrix<3, 4>;

	[[nodiscard]] static Element::StiffnessAndForces KT(
		Node::Positions const& x,
		Element::ShapeDerivative const& dN_by_dX,
		Scalar lambda,
		Scalar mu);

	[[nodiscard]] static Element::Stiffness Kc(
		Element::ShapeDerivative const& dN_by_dx,
		Scalar v,
		Element::Elasticity const& c);
	[[nodiscard]] static Element::Stiffness Ks(
		Element::ShapeDerivative const& dN_by_dx,
		Scalar v,
		Element::Stress const& s);

	[[nodiscard]] static Element::Elasticity c(Scalar J, Scalar lambda, Scalar mu);

	[[nodiscard]] static Node::Forces T(
		Element::ShapeDerivative const& dN_by_dx,
		Scalar v,
		Element::Stress const& sigma);
	[[nodiscard]] static Element::Stress sigma(
		Scalar J, Element::Gradient const& b, Scalar lambda, Scalar mu);

	[[nodiscard]] static Scalar J(Element::Gradient const& F);
	[[nodiscard]] static Element::Gradient b(Element::Gradient const& F);

	[[nodiscard]] static Element::Gradient dx_by_dX(
		Node::Positions const& x, Element::ShapeDerivative const& dN_by_dX);
	[[nodiscard]] static Element::Gradient dx_by_dX(
		Element::Gradient const& dx_by_dL, Element::Gradient const& dL_by_dX);

	[[nodiscard]] static Scalar V(Node::Positions const& x);

	[[nodiscard]] static Element::Gradient dX_by_dL(Node::Positions const& X);
	[[nodiscard]] static Element::Gradient dL_by_dX(
		Element::Gradient const& dX_by_dL);
	[[nodiscard]] static Element::CartesianDerivative dx_by_dN(
		Element::ShapeCartesianTransform const& N_to_x);

	[[nodiscard]] static Element::ShapeDerivative dN_by_dX(
		Element::Gradient const& dL_by_dx);
	[[nodiscard]] static Element::ShapeDerivative dN_by_dX(
		Element::ShapeCartesianTransform const& N_to_x);
	[[nodiscard]] static Element::ShapeDerivative dN_by_dX(
		Node::Positions const& X);

	[[nodiscard]] static Element::ShapeCartesianTransform N_to_x(
		Node::Positions const& X);
};
}  // namespace FeltElements
