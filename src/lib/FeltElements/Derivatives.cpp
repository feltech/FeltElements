#include "Derivatives.hpp"

#include "Body.hpp"

namespace
{
using namespace FeltElements;

template <typename T>
int sgn(T val)
{
	return (T(0) < val) - (val < T(0));
}
auto constexpr delta = [](auto const i, auto const j) { return i == j; };

Element::Elasticity const c_lambda = ([]() {  // NOLINT(cert-err58-cpp)
	std::size_t constexpr N = 3;
	Element::Elasticity c{};
	for (std::size_t i = 0; i < N; i++)
		for (std::size_t j = 0; j < N; j++)
			for (std::size_t k = 0; k < N; k++)
				for (std::size_t l = 0; l < N; l++) c(i, j, k, l) = delta(i, j) * delta(k, l);

	return c;
}());

// TODO: constexpr - requires constexpr std::copy (C++20)
Element::Elasticity const c_mu = ([]() {  // NOLINT(cert-err58-cpp)
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
Tensor::Matrix<3> const I = ([]() {	 // NOLINT(cert-err58-cpp)
	Tensor::Matrix<3> mat{};
	mat.eye2();
	return mat;
}());
}  // namespace

namespace FeltElements::ex
{
using namespace Tensor;

constexpr auto dX_by_dL = [](auto const & X) {
	return Func::einsum<Idxs<k, i>, Idxs<k, j>>(X, Derivatives::dN_by_dL);
};

constexpr auto dL_by_dX = [](auto const & dX_by_dL_) { return Func::inv(dX_by_dL_); };

constexpr auto dX_by_dS = [](auto const & X) {
	return Func::einsum<Idxs<k, i>, Idxs<k, j>>(X, Derivatives::dN_by_dS);
};

constexpr auto dN_by_dX = [](auto const & dL_by_dX_) {
	// dN/dX^T = dX/dL^(-T) * dN/dL^T => dN/dX = dN/dL * dX/dL^(-1) = dN/dL * dL/dX
	return Func::einsum<Idxs<i, k>, Idxs<k, j>>(Derivatives::dN_by_dL, dL_by_dX_);
};

constexpr auto dx_by_dX = [](auto const & x, auto const & dN_by_dX_) {
	return Func::einsum<Idxs<k, i>, Idxs<k, j>>(x, dN_by_dX_);
};

constexpr auto finger = [](auto const & F) {
	return Func::einsum<Idxs<i, k>, Idxs<j, k>>(F, F);
};

constexpr auto sigma =
	[](Scalar const J, auto const & b, Scalar const lambda, Scalar const mu) {
		return (mu / J) * (b - I) + (lambda / J) * log(J) * I;
	};

constexpr auto t = [](Scalar const p, auto const & dX_by_dS_) {
	using Tensor::Func::all;
	using Tensor::Func::fix;
	Node::Force const & dX1_by_dS = dX_by_dS_(all, fix<0>);
	Node::Force const & dX2_by_dS = dX_by_dS_(all, fix<1>);
	return (1.0 / 2.0) * p * cross(dX1_by_dS, dX2_by_dS);
};

constexpr auto T = [](auto const & dN_by_dx, auto const & sigma_) {
	// T = v * sigma * dN/dx^T
	return Func::einsum<Idxs<a, k>, Idxs<i, k>>(dN_by_dx, sigma_);
};

constexpr auto c = [](Scalar J, Scalar lambda, Scalar mu) {
	Scalar const lambda_prime = lambda / J;
	Scalar const mu_prime = (mu - lambda * std::log(J)) / J;

	return lambda_prime * c_lambda + mu_prime * c_mu;
};

constexpr auto Kc = [](auto const & dN_by_dx, auto const & c_) {
	// Kc_ij = v * dN_a/dx_k * c_ikjl * dN_b/dx_l
	return Func::einsum<Idxs<a, k>, Idxs<i, k, j, l>, Idxs<b, l>, Order<a, i, b, j>>(
		dN_by_dx, c_, dN_by_dx);
};

constexpr auto Ks = [](auto const & dN_by_dx, auto const & s) {
	// Ks_ij = v * dN_a/dx_k * sigma_kl * dN_b/dx_l * delta_ij
	return Func::einsum<Idxs<a, k>, Idxs<k, l>, Idxs<b, l>, Idxs<i, j>, Order<a, i, b, j>>(
		dN_by_dx, s, dN_by_dx, I);
};
}  // namespace FeltElements::ex

namespace FeltElements::Derivatives
{
Element::StiffnessResidual KR(
	Element::Positions const & x,
	Element::ShapeDerivative const & dN_by_dX,
	Body::Material const & material, Body::Forces const & forces)
{
	Element::Gradient const F = ex::dx_by_dX(x, dN_by_dX);
	auto const b = ex::finger(F);
	Scalar const J = Derivatives::det_dx_by_dX(F);
	Scalar const v = Derivatives::v(x); // TODO: but also v = J*V

	auto const & dx_by_dL = ex::dX_by_dL(x);
	auto const & dL_by_dx = ex::dL_by_dX(dx_by_dL);
	Element::ShapeDerivative const dN_by_dx = ex::dN_by_dX(dL_by_dx);

	Element::Stress const sigma = ex::sigma(J, b, material.lambda, material.mu);

	auto const & c = ex::c(J, material.lambda, material.mu);
	auto const & Kc = ex::Kc(dN_by_dx, c);
	auto const & Ks = ex::Ks(dN_by_dx, sigma);

	// Tangent stiffness matrix.
	Element::Stiffness K = v * (Kc + Ks);
	// Internal forces.
	auto const & T_by_v = ex::T(dN_by_dx, sigma);
	// External forces
	auto const & F_by_V = forces.F_by_m * material.rho;
	Node::Force const F_by_v_per_node = 1.0 / 4.0 * F_by_V / J;
	Element::Forces F_by_v;
	using Tensor::Func::all;
	F_by_v(0, all) = F_by_v_per_node;
	F_by_v(1, all) = F_by_v_per_node;
	F_by_v(2, all) = F_by_v_per_node;
	F_by_v(3, all) = F_by_v_per_node;

	Element::Forces const R = v * (T_by_v - F_by_v);

	assert(Tensor::Func::all_of(R == R));  // Assert no NaNs

	return {K, R};
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

Node::Force t(Scalar const p, Element::SurfaceGradient const & dX_by_dS)
{
	return ex::t(p, dX_by_dS);
}

Element::Forces T(
	Element::ShapeDerivative const & dN_by_dx, Scalar const v, Element::Stress const & sigma)
{
	return v * ex::T(dN_by_dx, sigma);
}

Element::Stress sigma(
	Scalar const J, Element::Gradient const & b, Scalar const lambda, Scalar const mu)
{
	return ex::sigma(J, b, lambda, mu);
}

Scalar det_dx_by_dX(Element::Gradient const & dx_by_dX)
{
	return Fastor::det(dx_by_dX);
}

Element::Gradient b(Element::Gradient const & F)
{
	return ex::finger(F);
}

Element::Gradient dx_by_dX(Element::Positions const & x, Element::ShapeDerivative const & dN_by_dX)
{
	return ex::dx_by_dX(x, dN_by_dX);
}

Element::Gradient dx_by_dX(Element::Gradient const & dx_by_dL, Element::Gradient const & dL_by_dX)
{
	using namespace Tensor;
	return Func::einsum<Idxs<k, i>, Idxs<j, k>>(dx_by_dL, dL_by_dX);
}

Element::ShapeDerivative dN_by_dX(Element::Gradient const & dL_by_dx)
{
	return ex::dN_by_dX(dL_by_dx);
}

Element::ShapeDerivative dN_by_dX(Element::Positions const & X)
{
	auto const & dX_by_dL = ex::dX_by_dL(X);
	auto const & dL_by_dX = ex::dL_by_dX(dX_by_dL);
	return ex::dN_by_dX(dL_by_dX);
}

Element::CartesianDerivative dx_by_dN(Element::ShapeCartesianTransform const & N_to_x)
{
	using namespace Tensor::Func;
	return N_to_x(fseq<1, last>{}, all);
}

Element::Gradient dL_by_dX(Element::Gradient const & dX_by_dL)
{
	return ex::dL_by_dX(dX_by_dL);
}

Element::Gradient dX_by_dL(Element::Positions const & X)
{
	return ex::dX_by_dL(X);
}

Element::SurfaceGradient dX_by_dS(BoundaryElement::Positions const & X)
{
	return ex::dX_by_dS(X);
}

Element::ShapeDerivative dN_by_dX(Element::ShapeCartesianTransform const & N_to_x)
{
	using namespace Tensor::Func;
	// Interpolation: (1, x, y, z)^T = N_to_x * N, where N is 4x natural coordinates (corners).
	// Invert then strip constant terms, leaving just coefficients, i.e. the derivative.
	return evaluate(inv(N_to_x))(all, fseq<1, last>{});
}

Element::ShapeCartesianTransform N_to_x(Element::Positions const & X)
{
	Element::ShapeCartesianTransform mat;
	mat(Fastor::ffirst, Fastor::all) = 1.0;
	mat(Fastor::fseq<1, Fastor::last>{}, Fastor::all) = Fastor::transpose(X);
	return mat;
}

Scalar V(Element::Positions const & x)
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

Scalar v(Element::Positions const & x)
{
	// Volume in local coords (constant)
	constexpr Scalar v_wrt_L = 1.0 / 6.0;
	return v_wrt_L * det_dx_by_dL(x);
}

Scalar det_dx_by_dL(Element::Positions const & x)
{
	using namespace Tensor;
	using Func::all;
	using Func::einsum;
	using Func::fix;

	//	Scalar det_ = 0;
	//	for (Index i = 0; i < Element::count; i++)
	//		for (Index j = 0; j < Element::count; j++)
	//			for (Index k = 0; k < Element::count; k++)
	//				det_ += det_dN_by_dL(i, j, k) * x(i, 0) * x(j, 1) * x(k, 2);
	//
	//	return det_;

	//	Element::ShapeDerivativeDeterminant xs = ([&x](){
	//		using Tensor::Index;
	//		using Tensor::Func::all;
	//		Element::ShapeDerivativeDeterminant xs_;
	//		for (Index i = 0; i < Element::count; i++)
	//			for (Index j = 0; j < Element::count; j++)
	//				for (Index k = 0; k < Element::count; k++)
	//					xs_(i, j, k) = x(i, 0) * x(j, 1) * x(k, 2);
	//		return xs_;
	//	}());
	//	Tensor::Vector<1> sum = einsum<Idxs<i, j, k>, Idxs<i, j, k>>(
	//		det_dN_by_dL, xs);

	// "Consideration of Body Forces within Finite Element Analysis", Glenk et al., 2018
	Tensor::Vector<4> x0 = x(all, fix<0>);
	Tensor::Vector<4> x1 = x(all, fix<1>);
	Tensor::Vector<4> x2 = x(all, fix<2>);
	Tensor::Vector<1> sum =
		einsum<Idxs<i, j, k>, Idxs<i>, Idxs<j>, Idxs<k>>(Derivatives::det_dN_by_dL, x0, x1, x2);

	return sum(0);
}

// clang-format off
Element::IsoCoordDerivative const dL_by_dN = // NOLINT(cert-err58-cpp)
	Tensor::Matrix<4, 4>{
		{1, 1, 1, 1},
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1}}(Fastor::fseq<1, 4>(), Fastor::all);

Element::ShapeDerivative const dN_by_dL = // NOLINT(cert-err58-cpp)
	Fastor::evaluate(Fastor::inv(Tensor::Matrix<4, 4>{
		{1, 1, 1, 1},
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1}}))(Fastor::all, Fastor::fseq<1, 4>());

Element::SurfaceShapeDerivative const dN_by_dS = // NOLINT(cert-err58-cpp)
	Fastor::evaluate(Fastor::inv(Tensor::Matrix<3, 3>{
		{1, 1, 1},
		{0, 1, 0},
		{0, 0, 1},
		}))(Fastor::all, Fastor::fseq<1, 3>());
// clang-format on

Element::ShapeDerivativeDeterminant const det_dN_by_dL = ([]() {  // NOLINT(cert-err58-cpp)
	using Tensor::Index;
	using Tensor::Func::all;
	using Tensor::Func::det;
	Element::ShapeDerivativeDeterminant det_dN_by_dL_;
	for (Index i = 0; i < Element::count; i++)
		for (Index j = 0; j < Element::count; j++)
			for (Index k = 0; k < Element::count; k++)
			{
				Tensor::Matrix<3> dN_by_dL_ijk;
				dN_by_dL_ijk(0, all) = dN_by_dL(i, all);
				dN_by_dL_ijk(1, all) = dN_by_dL(j, all);
				dN_by_dL_ijk(2, all) = dN_by_dL(k, all);
				det_dN_by_dL_(i, j, k) = det(dN_by_dL_ijk);
			}
	return det_dN_by_dL_;
}());

Tensor::Multi<Node::dim, Node::dim, Node::dim> const levi_civita = // NOLINT(cert-err58-cpp)
	([]() {
		using LeviCivita = Tensor::Multi<Node::dim, Node::dim, Node::dim>;
		LeviCivita E;
		int constexpr const count = static_cast<int>(Node::dim);
		for (int i = 0; i < count; i++)
			for (int j = 0; j < count; j++)
				for (int k = 0; k < count; k++) E(i, j, k) = sgn(j - i) * sgn(k - i) * sgn(k - j);
		return E;
	}());
}  // namespace FeltElements::Derivatives
