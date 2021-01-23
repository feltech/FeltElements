#pragma once

#include <Eigen/Dense>
#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <cstddef>

#include "Typedefs.hpp"

namespace FeltElements
{
// Forward declarations.
namespace Body
{
struct Material;
struct Forces;
}  // namespace Body

namespace Derivatives
{
[[nodiscard]] Element::StiffnessResidual KR(
	Element::NodePositions const & x,
	Element::BoundaryVtxhIdxs const & boundary_faces_idxs,
	Element::BoundaryNodePositions const & boundary_faces_x,
	Element::ShapeDerivative const & dN_by_dX,
	Body::Material const & material,
	Body::Forces const & forces);

[[nodiscard]] Element::Stiffness Kc(
	Element::ShapeDerivative const & dN_by_dx, Scalar v, Element::Elasticity const & c);
[[nodiscard]] Element::Stiffness Ks(
	Element::ShapeDerivative const & dN_by_dx, Scalar v, Element::Stress const & s);
[[nodiscard]] Element::Stiffness Kp(
	Element::BoundaryNodePositions const & xs, Element::BoundaryVtxhIdxs const & S_to_Vs, Scalar p);

[[nodiscard]] Element::Elasticity c(Scalar J, Scalar lambda, Scalar mu);
[[nodiscard]] Node::Force t(Scalar p, Element::SurfaceGradient const & dX_by_dS);
[[nodiscard]] Element::Forces T(
	Element::ShapeDerivative const & dN_by_dx, Scalar v, Element::Stress const & sigma);
[[nodiscard]] Element::Stress sigma(
	Scalar J, Element::Gradient const & b, Scalar lambda, Scalar mu);

[[nodiscard]] Scalar det_dx_by_dX(Element::Gradient const & F);
[[nodiscard]] Element::Gradient b(Element::Gradient const & F);

[[nodiscard]] Element::Gradient dx_by_dX(
	Element::NodePositions const & x, Element::ShapeDerivative const & dN_by_dX);
[[nodiscard]] Element::Gradient dx_by_dX(
	Element::Gradient const & dx_by_dL, Element::Gradient const & dL_by_dX);

[[nodiscard]] Scalar V(Element::NodePositions const & x);
[[nodiscard]] Scalar v(Element::NodePositions const & x);
[[nodiscard]] Scalar A(BoundaryElement::NodePositions const & s);

[[nodiscard]] Element::Gradient dX_by_dL(Element::NodePositions const & X);
[[nodiscard]] Element::Gradient dL_by_dX(Element::Gradient const & dX_by_dL);
[[nodiscard]] Element::CartesianDerivative dx_by_dN(
	Element::ShapeCartesianTransform const & N_to_x);

[[nodiscard]] Element::SurfaceGradient dX_by_dS(BoundaryElement::NodePositions const & X);

[[nodiscard]] Element::ShapeDerivative dN_by_dX(Element::Gradient const & dL_by_dx);
[[nodiscard]] Element::ShapeDerivative dN_by_dX(Element::ShapeCartesianTransform const & N_to_x);
[[nodiscard]] Element::ShapeDerivative dN_by_dX(Element::NodePositions const & X);

[[nodiscard]] Element::ShapeCartesianTransform N_to_x(Element::NodePositions const & X);

[[nodiscard]] Scalar det_dx_by_dL(Element::NodePositions const & x);

extern Element::IsoCoordDerivative const dL_by_dN;
extern Element::ShapeDerivative const dN_by_dL;
extern Element::SurfaceShapeDerivative const dN_by_dS;
extern Element::ShapeDerivativeDeterminant const det_dN_by_dL;
extern Tensor::Multi<3, 3, 3> const levi_civita;
}  // namespace Derivatives
}  // namespace FeltElements