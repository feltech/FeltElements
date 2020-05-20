#include "TetGenIO.hpp"
#include "Tetrahedron.hpp"
#include "internal/Conversions.hpp"
#include "Attributes.hpp"

namespace
{
using namespace FeltElements;

auto const delta = [](auto const i, auto const j)
{
	return i == j;
};

template <typename Matrix>
using is_ovm_vector = typename std::enable_if<
	(std::is_same<Matrix, typename Matrix::vector_type>::value), void*>::type;


template <typename Matrix, is_ovm_vector<Matrix> = nullptr >
auto to_tensor(Matrix const & matrix)
{
	return Eigen::TensorMap<Eigen::Tensor<Scalar const, 1>>{
		matrix.data(), {static_cast<Eigen::Index>(matrix.size())}
	};
}

template <class Tensor, std::size_t idx>
Eigen::Index const size = Eigen::internal::get<idx, typename Tensor::Dimensions::Base>::value;

template <typename Tensor>
auto to_matrix(Tensor const & tensor)
{
	return Eigen::Map<
		Eigen::Matrix<double, size<Tensor, 0>, size<Tensor, 1>, Tensor::Layout> const>{
		tensor.data(), size<Tensor, 0>, size<Tensor, 1>};
}

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

template <std::size_t dim = Node::dim>
const Tensor::Matrix<dim, dim> I = internal::to_tensor(
	Tetrahedron::Matrix<dim, dim>{Tetrahedron::Matrix<dim, dim>::Identity()});
} // End anon namespace

namespace FeltElements
{
}

namespace FeltElements::ex
{

auto const dX_by_dL = [](auto const & X) {
	constexpr Tensor::IndexPairs<1> X_L{{{0, 0}}};
	return X.contract(Element::Attribute::MaterialShapeDerivative::dN_by_dL, X_L);
};

auto const dL_by_dX = [](Element::Gradient const & dX_by_dL)
	-> Element::Gradient {
	// Note:
	// * Can't accept an expression because may not have `Dimensions` for `to_matrix()` to use.
	// * Can't return expression because `TensorMap`ing a temporary.
	return internal::to_tensor(Tetrahedron::GradientMatrix{to_matrix(dX_by_dL).inverse()});
};

auto const dN_by_dX = [](auto const & dL_by_dX) {
  // dN/dX^T = dX/dL^(-T) * dN/dL^T => dN/dX = dN/dL * dX/dL^(-1) = dN/dL * dL/dX
  constexpr Tensor::IndexPairs<1> dN_dL{{{1, 0}}};
  return Element::Attribute::MaterialShapeDerivative::dN_by_dL.contract(dL_by_dX, dN_dL);
};

auto const dx_by_dX = [](auto const & x, auto const & dN_by_dX) {
	constexpr Tensor::IndexPairs<1> x_dN{{{0, 0}}};
	return x.contract(dN_by_dX, x_dN);
};

auto const b = [](auto const & F) {
	constexpr Tensor::IndexPairs<1> F_F{{{1, 1}}};
	return F.contract(F, F_F);
};

auto const sigma = [](
	Scalar const J, auto const & b,
	Scalar const lambda, Scalar const mu) {
	return (mu / J) * (b - I<>) + (lambda / J) * log(J) * I<>;
};

auto const T = [](auto const & dN_by_dx, Scalar const v, auto const & sigma) {
	// T = v * sigma * dN/dx^T
	constexpr Tensor::IndexPairs<1> sigma_dN{{{1, 1}}};
	return v * sigma.contract(dN_by_dx, sigma_dN);
};

auto const c = [](Scalar J, Scalar lambda, Scalar mu) {
	Scalar const lambda_prime = lambda / J;
	Scalar const mu_prime = (mu - lambda * std::log(J)) / J;

	return lambda_prime * c_lambda + mu_prime * c_mu;
};

auto const Kc = [](auto const & dN_by_dx, Scalar const v, auto const & c) {
	constexpr Tensor::IndexPairs<1> c_a{{{1, 1}}};
	constexpr Tensor::IndexPairs<1> c_b{{{3, 1}}};
	using Shuffle = Eigen::array<Eigen::Index, 4>;
	constexpr Shuffle shuffle{1, 0, 2, 3};
	// Kc = v * dN_a/dx_k * c_ikjl * dN_b/dx_l
	return v * dN_by_dx.contract(c, c_a).contract(dN_by_dx, c_b).shuffle(shuffle);
};

auto const Ks = [](auto const & dN_by_dx, Scalar const v, auto const & s)
{
	constexpr Tensor::IndexPairs<1> s_a{{{1, 0}}};
	constexpr Tensor::IndexPairs<1> s_b{{{1, 1}}};
	constexpr Tensor::IndexPairs<0> s_I{};
	using Shuffle = Eigen::array<Eigen::Index, 4>;
	constexpr Shuffle shuffle{2, 0, 3, 1};
	// Ks = v * dN_a/dx_k * sigma_kl * dN_b/dx_l * delta_ij
	return v * dN_by_dx.contract(s, s_a).contract(dN_by_dx, s_b).contract(I<>, s_I).shuffle(shuffle);
};
}


namespace FeltElements
{
Element::StiffnessAndForces Tetrahedron::KT(
	Node::Positions const & x, Element::ShapeDerivative const & dN_by_dX,
	Scalar const lambda, Scalar const mu)
{
	Scalar const v = Tetrahedron::V(x);

	Element::Gradient const & F = ex::dx_by_dX(x, dN_by_dX);  // non-expression for calc'ing Jacobian
	auto const & b = ex::b(F);
	Scalar const J = Tetrahedron::J(F);

	auto const & dx_by_dL = ex::dX_by_dL(x);
	auto const & dL_by_dx = ex::dL_by_dX(dx_by_dL);
	auto const & dN_by_dx = ex::dN_by_dX(dL_by_dx);

	auto const & sigma = ex::sigma(J, b, lambda, mu); // Evaluate since used twice?

	auto const & c = ex::c(J, lambda, mu);
	auto const & Kc = ex::Kc(dN_by_dx, v, c);

	auto const & Ks = ex::Ks(dN_by_dx, v, sigma);

	auto const & K = Kc + Ks;
	auto const & T = ex::T(dN_by_dx, v, sigma);

	return Element::StiffnessAndForces(K, T);
}

Element::Stiffness Tetrahedron::Kc(
	Element::ShapeDerivative const & dN_by_dx, Scalar const v, Element::Elasticity const & c)
{
	return ex::Kc(dN_by_dx, v, c);
}

Element::Stiffness Tetrahedron::Ks(
	Element::ShapeDerivative const & dN_by_dx, Scalar const v, Element::Stress const & s)
{
	return ex::Ks(dN_by_dx, v, s);
}

Element::Elasticity Tetrahedron::c(Scalar J, Scalar lambda, Scalar mu)
{
	return ex::c(J, lambda, mu);
}

Node::Forces Tetrahedron::T(
	Element::ShapeDerivative const & dN_by_dx, Scalar const v, Element::Stress const & sigma)
{
	return ex::T(dN_by_dx, v, sigma);
}

Element::Stress Tetrahedron::sigma(
	Scalar const J, Element::Gradient const & b,
	Scalar const lambda, Scalar const mu)
{
	return ex::sigma(J, b, lambda, mu);
}

Scalar Tetrahedron::J(Element::Gradient const & dx_by_dX)
{
	return to_matrix(dx_by_dX).determinant();
}

Element::Gradient Tetrahedron::b(Element::Gradient const & F)
{
	return ex::b(F);
}

Element::Gradient Tetrahedron::dx_by_dX(
	Node::Positions const & x, Element::ShapeDerivative const & dN_by_dX)
{
	return ex::dx_by_dX(x, dN_by_dX);
}

Element::Gradient Tetrahedron::dx_by_dX(
	Element::Gradient const & dx_by_dL, Element::Gradient const & dL_by_dX)
{
	constexpr Tensor::IndexPairs<1> x_X{{{0, 1}}};
	return dx_by_dL.contract(dL_by_dX, x_X);
}

Element::ShapeDerivative Tetrahedron::dN_by_dX(Node::Positions const & X)
{
	auto const & dX_by_dL = ex::dX_by_dL(X);
	auto const & dL_by_dX = ex::dL_by_dX(dX_by_dL);
	return ex::dN_by_dX(dL_by_dX);
}

Element::CartesianDerivative Tetrahedron::dx_by_dN(
	Element::ShapeCartesianTransform const & N_to_x)
{
	Matrix<3, 4> dx_by_dN = to_matrix(N_to_x).block<3, 4>(1, 0);
	return internal::to_tensor(dx_by_dN);
}

Element::Gradient Tetrahedron::dL_by_dX(Element::Gradient const & dX_by_dL)
{
	return ex::dL_by_dX(dX_by_dL);
}

Element::Gradient Tetrahedron::dX_by_dL(Node::Positions const & X)
{
	return ex::dX_by_dL(X);
}

Element::ShapeDerivative Tetrahedron::dN_by_dX(Element::ShapeCartesianTransform const & N_to_x)
{
	// Interpolation: (1, x, y, z)^T = N_to_x * N, where N is 4x natural coordinates (corners).
	// Invert then strip constant terms, leaving just coefficients, i.e. the derivative.
	Matrix<4, 3> dN_by_dX = to_matrix(N_to_x).inverse().block<4, 3>(0, 1);
	return internal::to_tensor(dN_by_dX);
}

Element::ShapeCartesianTransform Tetrahedron::N_to_x(Node::Positions const & X)
{
	using Shuffle = Eigen::array<Eigen::Index, 2>;
	using Padding = Eigen::array<std::pair<Eigen::Index, Eigen::Index>, 2>;
	constexpr Shuffle transpose{{1, 0}};
	constexpr Padding padding{{{1, 0}, {0, 0}}};
	return X.shuffle(transpose).pad(padding, 1.0);
}

Scalar Tetrahedron::V(Node::Positions const & x)
{
	using Indices = Eigen::array<Eigen::Index, 2>;
	auto const & start_3x3 = x.slice(Indices{0, 0}, Indices{3, 3});
	auto const & end_1x3 = x.slice(Indices{3, 0}, Indices{1, 3});
	auto const & end_3x3 = end_1x3.broadcast(Indices{3, 1});

	Tensor::Matrix<3, 3> const delta = start_3x3 - end_3x3;
	auto const & mat_delta = to_matrix(delta);
	return std::abs(mat_delta.determinant() / 6.0);
}

Tetrahedron::SpatialCoordProp Tetrahedron::x(Mesh & mesh)
{
	SpatialCoordProp x_prop =
		mesh.request_vertex_property<OpenVolumeMesh::Vec3d>("x");
	x_prop->set_persistent(true);

	for (
		OpenVolumeMesh::VertexIter it_vtx = mesh.vertices_begin(); it_vtx != mesh.vertices_end();
		it_vtx++)
		x_prop[*it_vtx] = mesh.vertex(*it_vtx);

	return x_prop;
}

Node::Positions Tetrahedron::x(
	Vtxhs const & vtxhs, SpatialCoordProp const & x_prop)
{
	Node::Positions p;
	std::size_t node_idx = 0;
	for (Vtxh const & vtxh : vtxhs)
	{
		OpenVolumeMesh::Vec3d const & vtx = x_prop[vtxh];
		p.chip(node_idx++, 0) = to_tensor(vtx);
	}
	return p;
}

Node::Positions Tetrahedron::X(Mesh const & mesh, Vtxhs const & vtxhs)
{
	Node::Positions p;
	std::size_t node_idx = 0;
	for (Vtxh const & vtxh : vtxhs)
	{
		OpenVolumeMesh::Vec3d const & vtx = mesh.vertex(vtxh);
		p.chip(node_idx++, 0) = to_tensor(vtx);
	}
	return p;
}

Tetrahedron::Vtxhs Tetrahedron::vtxhs(
	OpenVolumeMesh::GeometricTetrahedralMeshV3d const & mesh,
	OpenVolumeMesh::CellHandle const & cellh)
{
	Vtxhs avtxhs;
	auto const vvtxhs = mesh.get_cell_vertices(cellh);
	std::copy_n(vvtxhs.begin(), avtxhs.size(), avtxhs.begin());
	return avtxhs;
}
} // namespace FeltElements
