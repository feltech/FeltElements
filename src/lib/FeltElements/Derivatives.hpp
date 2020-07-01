#pragma once

#include <Eigen/Dense>
#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <cstddef>

#include "Typedefs.hpp"

namespace FeltElements::Derivatives
{
[[nodiscard]] Element::StiffnessAndForces KT(
	Node::Positions const& x,
	Element::ShapeDerivative const& dN_by_dX,
	Scalar lambda,
	Scalar mu);

[[nodiscard]] Element::Stiffness Kc(
	Element::ShapeDerivative const& dN_by_dx,
	Scalar v,
	Element::Elasticity const& c);
[[nodiscard]] Element::Stiffness Ks(
	Element::ShapeDerivative const& dN_by_dx,
	Scalar v,
	Element::Stress const& s);

[[nodiscard]] Element::Elasticity c(Scalar J, Scalar lambda, Scalar mu);

[[nodiscard]] Node::Forces T(
	Element::ShapeDerivative const& dN_by_dx,
	Scalar v,
	Element::Stress const& sigma);
[[nodiscard]] Element::Stress sigma(
	Scalar J, Element::Gradient const& b, Scalar lambda, Scalar mu);

[[nodiscard]] Scalar J(Element::Gradient const& F);
[[nodiscard]] Element::Gradient b(Element::Gradient const& F);

[[nodiscard]] Element::Gradient dx_by_dX(
	Node::Positions const& x, Element::ShapeDerivative const& dN_by_dX);
[[nodiscard]] Element::Gradient dx_by_dX(
	Element::Gradient const& dx_by_dL, Element::Gradient const& dL_by_dX);

[[nodiscard]] Scalar V(Node::Positions const& x);

[[nodiscard]] Element::Gradient dX_by_dL(Node::Positions const& X);
[[nodiscard]] Element::Gradient dL_by_dX(
	Element::Gradient const& dX_by_dL);
[[nodiscard]] Element::CartesianDerivative dx_by_dN(
	Element::ShapeCartesianTransform const& N_to_x);

[[nodiscard]] Element::ShapeDerivative dN_by_dX(
	Element::Gradient const& dL_by_dx);
[[nodiscard]] Element::ShapeDerivative dN_by_dX(
	Element::ShapeCartesianTransform const& N_to_x);
[[nodiscard]] Element::ShapeDerivative dN_by_dX(
	Node::Positions const& X);

[[nodiscard]] Element::ShapeCartesianTransform N_to_x(
	Node::Positions const& X);

[[nodiscard]] Scalar det_dx_by_dL(Node::Positions const& x);
}  // namespace FeltElements
