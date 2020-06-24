#include "Derivatives.hpp"

#include "Attributes.hpp"
#include "TetGenIO.hpp"

namespace
{
using namespace FeltElements;

auto constexpr delta = [](auto const i, auto const j)
{
	return i == j;
};

template <typename Matrix>
using is_ovm_vector = typename std::enable_if<
	(std::is_same<Matrix, typename Matrix::vector_type>::value), void*>::type;

Element::Elasticity const c_lambda = ([]() { // NOLINT(cert-err58-cpp)
  std::size_t constexpr N = 3;
  Element::Elasticity c{};
  for (std::size_t i = 0; i < N; i++)
	  for (std::size_t j = 0; j < N; j++)
		  for (std::size_t k = 0; k < N; k++)
			  for (std::size_t l = 0; l < N; l++)
				  c(i, j, k, l) = delta(i, j) * delta(k, l);

  return c;
}());

// TODO: constexpr - requires constexpr std::copy (C++20)
Element::Elasticity const c_mu = ([]() { // NOLINT(cert-err58-cpp)
  std::size_t constexpr N = 3;
  Element::Elasticity c{};
  for (std::size_t i = 0; i < N; i++)
	  for (std::size_t j = 0; j < N; j++)
		  for (std::size_t k = 0; k < N; k++)
			  for (std::size_t l = 0; l < N; l++)
				  c(i, j, k, l) = delta(i, k) * delta(j, l) + delta(i, l) * delta(j, k);

  return c;
}());

// TODO: constexpr - requries simd_vector_type to satisfy literal type requirements
Tensor::Matrix<3> const I = ([]() {
	Tensor::Matrix<3> mat{};
	mat.eye2();
	return mat;
}());
} // End anon namespace

namespace FeltElements::ex
{
using namespace Tensor;

auto const dX_by_dL = [](auto const & X) {
	return Func::einsum<Idxs<k, i>, Idxs<k, j>>(
		X, Element::Attribute::MaterialShapeDerivative::dN_by_dL);
};

auto const dL_by_dX = [](auto const & dX_by_dL_) {
	return Func::inv(dX_by_dL_);
};

auto const dN_by_dX = [](auto const & dL_by_dX_) {
	// dN/dX^T = dX/dL^(-T) * dN/dL^T => dN/dX = dN/dL * dX/dL^(-1) = dN/dL * dL/dX
	return Func::einsum<Idxs<i, k>, Idxs<k, j>>(
		Element::Attribute::MaterialShapeDerivative::dN_by_dL, dL_by_dX_);
};

auto const dx_by_dX = [](auto const & x, auto const & dN_by_dX_) {
	return Func::einsum<Idxs<k, i>, Idxs<k, j>>(x, dN_by_dX_);
};

auto const finger = [](auto const & F) {
	return Func::einsum<Idxs<i, k>, Idxs<j, k>>(F, F);
};

auto const sigma = [](
	Scalar const J, auto const & b,
	Scalar const lambda, Scalar const mu) {
	return (mu / J) * (b - I) + (lambda / J) * log(J) * I;
};

auto const T = [](auto const & dN_by_dx, Scalar const v, auto const & sigma_) {
	// T = v * sigma * dN/dx^T
	return v * Func::einsum<Idxs<a, k>, Idxs<i, k>>(dN_by_dx, sigma_);
};

auto const c = [](Scalar J, Scalar lambda, Scalar mu) {
	Scalar const lambda_prime = lambda / J;
	Scalar const mu_prime = (mu - lambda * std::log(J)) / J;

	return lambda_prime * c_lambda + mu_prime * c_mu;
};


auto const Kc = [](auto const & dN_by_dx, auto const & c_) {

  // Kc = v * dN_a/dx_k * c_ikjl * dN_b/dx_l
	return Func::einsum<Idxs<a, k>, Idxs<i, k, j, l>, Idxs<b, l>, Order<a, i, b, j>>(
	  dN_by_dx, c_, dN_by_dx);
};


auto const Ks = [](auto const & dN_by_dx, auto const & s)
{
	// Ks = v * dN_a/dx_k * sigma_kl * dN_b/dx_l * delta_ij
	return Func::einsum<Idxs<a, k>, Idxs<k, l>, Idxs<b, l>, Idxs<i, j>, Order<a, i, b, j>>(
		dN_by_dx, s, dN_by_dx, I);
};
}


namespace FeltElements::Derivatives
{
Element::StiffnessAndForces KT(
	Node::Positions const & x, Element::ShapeDerivative const & dN_by_dX,
	Scalar const lambda, Scalar const mu)
{
	Scalar const v = V(x);

	Element::Gradient const F = ex::dx_by_dX(x, dN_by_dX);
	auto const b = ex::finger(F);
	Scalar const J = Derivatives::J(F);

	auto const & dx_by_dL = ex::dX_by_dL(x);
	auto const & dL_by_dx = ex::dL_by_dX(dx_by_dL);
	Element::ShapeDerivative const dN_by_dx = ex::dN_by_dX(dL_by_dx);

	Element::Stress const sigma = ex::sigma(J, b, lambda, mu);

	auto const & c = ex::c(J, lambda, mu);
	auto const & Kc = ex::Kc(dN_by_dx, c);
	auto const & Ks = ex::Ks(dN_by_dx, sigma);

	Element::Stiffness K = v * (Kc + Ks);
	Node::Forces T = ex::T(dN_by_dx, v, sigma);

	return Element::StiffnessAndForces(K, T);
}

Element::Stiffness Kc(
	Element::ShapeDerivative const & dN_by_dx, Scalar const v, Element::Elasticity const & c)
{
	return v * ex::Kc(dN_by_dx, c);
}

Element::Stiffness Ks(
	Element::ShapeDerivative const & dN_by_dx, Scalar const v, Element::Stress const & s)
{
	return v * ex::Ks(dN_by_dx, s);
}

Element::Elasticity c(Scalar J, Scalar lambda, Scalar mu)
{
	return ex::c(J, lambda, mu);
}

Node::Forces T(
	Element::ShapeDerivative const & dN_by_dx, Scalar const v, Element::Stress const & sigma)
{
	return ex::T(dN_by_dx, v, sigma);
}

Element::Stress sigma(
	Scalar const J, Element::Gradient const & b,
	Scalar const lambda, Scalar const mu)
{
	return ex::sigma(J, b, lambda, mu);
}

Scalar J(Element::Gradient const & dx_by_dX)
{
	return Fastor::det(dx_by_dX);
}

Element::Gradient b(Element::Gradient const & F)
{
	return ex::finger(F);
}

Element::Gradient dx_by_dX(
	Node::Positions const & x, Element::ShapeDerivative const & dN_by_dX)
{
	return ex::dx_by_dX(x, dN_by_dX);
}

Element::Gradient dx_by_dX(
	Element::Gradient const & dx_by_dL, Element::Gradient const & dL_by_dX)
{
	using namespace Tensor;
	return Func::einsum<Idxs<k, i>, Idxs<j, k>>(dx_by_dL, dL_by_dX);
}

Element::ShapeDerivative dN_by_dX(Element::Gradient const& dL_by_dx)
{
	return ex::dN_by_dX(dL_by_dx);
}

Element::ShapeDerivative dN_by_dX(Node::Positions const & X)
{
	auto const & dX_by_dL = ex::dX_by_dL(X);
	auto const & dL_by_dX = ex::dL_by_dX(dX_by_dL);
	return ex::dN_by_dX(dL_by_dX);
}

Element::CartesianDerivative dx_by_dN(
	Element::ShapeCartesianTransform const & N_to_x)
{
	using namespace Tensor::Func;
	return N_to_x(fseq<1, last>{}, all);
}

Element::Gradient dL_by_dX(Element::Gradient const & dX_by_dL)
{
	return ex::dL_by_dX(dX_by_dL);
}

Element::Gradient dX_by_dL(Node::Positions const & X)
{
	return ex::dX_by_dL(X);
}

Element::ShapeDerivative dN_by_dX(Element::ShapeCartesianTransform const & N_to_x)
{
	using namespace Tensor::Func;
	// Interpolation: (1, x, y, z)^T = N_to_x * N, where N is 4x natural coordinates (corners).
	// Invert then strip constant terms, leaving just coefficients, i.e. the derivative.
	return evaluate(inv(N_to_x))(all, fseq<1, last>{});
}

Element::ShapeCartesianTransform N_to_x(Node::Positions const & X)
{
	Element::ShapeCartesianTransform mat;
	mat(Fastor::ffirst, Fastor::all) = 1.0;
	mat(Fastor::fseq<1, Fastor::last>{}, Fastor::all) = Fastor::transpose(X);
	return mat;
}

Scalar V(Node::Positions const & x)
{
	using namespace Tensor::Func;
	auto const & start_3x3 = x(fseq<0, 3>{}, all);
//	std::clog << x(fix<3>, all).self().size() << "\n";
	Tensor::Vector<3> end_1x3 = x(fix<3>, all);
	Tensor::Matrix<3> end_3x3;
	end_3x3(fix<0>, all) = end_1x3;
	end_3x3(fix<1>, all) = end_1x3;
	end_3x3(fix<2>, all) = end_1x3;
	auto const & delta = start_3x3 - end_3x3;
	return std::abs(Fastor::det(delta) / 6.0);
}

} // namespace FeltElements
