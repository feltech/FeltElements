//
// Created by dave on 07/11/2019.
//
#include "TetGenIO.hpp"
#include "Tetrahedron.hpp"

namespace
{
using namespace FeltElements;

template <typename Index>
auto delta(Index const i, Index const j)
{
	return i == j;
};

template <typename Matrix>
using is_matrix = typename std::enable_if<
    (Matrix::RowsAtCompileTime > 1 && Matrix::ColsAtCompileTime > 1), void*>::type;

template <typename Matrix>
using is_vector = typename std::enable_if<
	(Matrix::RowsAtCompileTime == 1 || Matrix::ColsAtCompileTime == 1), void*>::type;

template <typename Matrix>
using is_ovm_vector = typename std::enable_if<
	(std::is_same<Matrix, typename Matrix::vector_type>::value), void*>::type;

template <typename Matrix, is_matrix<Matrix> = nullptr>
auto to_tensor(Matrix const & matrix)
{
	static_assert(
		!(Matrix::Options & Eigen::RowMajor),
		"Cannot map a row-major matrix to a tensor since row-major tensors are not yet"
  		" supported by Eigen");

	return Eigen::TensorMap<
	    Eigen::Tensor<FeltElements::Tetrahedron::Scalar const, 2,
	    Matrix::Options & Eigen::RowMajor ? Eigen::RowMajor : Eigen::ColMajor>>{
		matrix.data(), {Matrix::RowsAtCompileTime, Matrix::ColsAtCompileTime}
	};
}

template <typename Matrix, is_vector<Matrix> = nullptr >
auto to_tensor(Matrix const & matrix)
{
	return Eigen::TensorMap<Eigen::Tensor<FeltElements::Tetrahedron::Scalar const, 1>>{
		matrix.data(), {std::max(Matrix::RowsAtCompileTime, Matrix::ColsAtCompileTime)}
	};
}

template <typename Matrix, is_ovm_vector<Matrix> = nullptr >
auto to_tensor(Matrix const & matrix)
{
	return Eigen::TensorMap<Eigen::Tensor<FeltElements::Tetrahedron::Scalar const, 1>>{
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

const Tetrahedron::ElasticityTensor c_lambda = ([]() {
  std::size_t constexpr N = 3;
  Tetrahedron::ElasticityTensor c{};
  for (std::size_t i = 0; i < N; i++)
	  for (std::size_t j = 0; j < N; j++)
		  for (std::size_t k = 0; k < N; k++)
			  for (std::size_t l = 0; l < N; l++)
				  c(i, j, k, l) = delta(i, j) * delta(k, l);

  return c;
}());

const Tetrahedron::ElasticityTensor c_mu = ([]() {
  std::size_t constexpr N = 3;
  Tetrahedron::ElasticityTensor c{};
  for (std::size_t i = 0; i < N; i++)
	  for (std::size_t j = 0; j < N; j++)
		  for (std::size_t k = 0; k < N; k++)
			  for (std::size_t l = 0; l < N; l++)
				  c(i, j, k, l) = delta(i, k) * delta(j, l) + delta(i, l) * delta(j, k);

  return c;
}());

template <std::size_t dim = Tetrahedron::Node::dim>
const Tetrahedron::MatrixTensor<dim, dim> I = to_tensor(
	Tetrahedron::Matrix<dim, dim>{Tetrahedron::Matrix<dim, dim>::Identity()});
} // End anon namespace

namespace FeltElements
{
// clang-format off
Tetrahedron::IsoCoordDerivativeTensor const Tetrahedron::dL_by_dN = // NOLINT(cert-err58-cpp)
	to_tensor(IsoCoordDerivativeMatrix{(Eigen::Matrix4d{} <<
		// (1, L) = A * N
		1, 1, 1, 1,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1).finished().block<3, 4>(1, 0)});

Tetrahedron::ShapeDerivativeTensor const Tetrahedron::dN_by_dL = // NOLINT(cert-err58-cpp)
	to_tensor(ShapeDerivativeMatrix{(Eigen::Matrix4d{} <<
		1, 1, 1, 1,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1).finished().inverse().block<4, 3>(0, 1)});
// clang-format on
}

namespace ex
{
using namespace FeltElements;


auto const dX_by_dL = [](auto const & X) {
	constexpr Tetrahedron::IndexPairs<1> X_L{{{0, 0}}};
	return X.contract(Tetrahedron::dN_by_dL, X_L);
};

auto const dL_by_dX = [](Tetrahedron::GradientTensor const & dX_by_dL)
	-> Tetrahedron::GradientTensor {
	// Note:
	// * Can't accept an expression because may not have `Dimensions` for `to_matrix()` to use.
	// * Can't return expression because mapping a temporary.
	return to_tensor(Tetrahedron::GradientMatrix{to_matrix(dX_by_dL).inverse()});
};

auto const dN_by_dX = [](auto const & dL_by_dX) {
  // dN/dX^T = dX/dL^(-T) * dN/dL^T => dN/dX = dN/dL * dX/dL^(-1) = dN/dL * dL/dX
  constexpr Tetrahedron::IndexPairs<1> dN_dL{{{1, 0}}};
  return Tetrahedron::dN_by_dL.contract(dL_by_dX, dN_dL);
};

auto const dx_by_dX = [](auto const & x, auto const & dN_by_dX) {
	constexpr Tetrahedron::IndexPairs<1> x_dN{{{0, 0}}};
	return x.contract(dN_by_dX, x_dN);
};

auto const b = [](auto const & F) {
	constexpr Tetrahedron::IndexPairs<1> F_F{{{1, 1}}};
	return F.contract(F, F_F);
};

auto const sigma = [](
	Tetrahedron::Scalar const J, auto const & b,
	Tetrahedron::Scalar const lambda, Tetrahedron::Scalar const mu) {
	return (mu / J) * (b - I<>) + (lambda / J) * log(J) * I<>;
};

auto const T = [](auto const & dN_by_dx, Tetrahedron::Scalar const v, auto const & sigma) {
	// T = v * sigma * dN/dx^T
	constexpr Tetrahedron::IndexPairs<1> sigma_dN{{{1, 1}}};
	return v * sigma.contract(dN_by_dx, sigma_dN);
};

auto const c = [](Tetrahedron::Scalar J, Tetrahedron::Scalar lambda, Tetrahedron::Scalar mu) {
	Tetrahedron::Scalar const lambda_prime = lambda / J;
	Tetrahedron::Scalar const mu_prime = (mu - lambda * std::log(J)) / J;

	return lambda_prime * c_lambda + mu_prime * c_mu;
};

auto const Kc = [](auto const & dN_by_dx, Tetrahedron::Scalar const v, auto const & c) {
	constexpr Tetrahedron::IndexPairs<1> c_a{{{1, 1}}};
	constexpr Tetrahedron::IndexPairs<1> c_b{{{3, 1}}};
	using Shuffle = Eigen::array<Eigen::Index, 4>;
	constexpr Shuffle shuffle{1, 0, 2, 3};
	// Kc = v * dN_a/dx_k * c_ikjl * dN_b/dx_l
	return v * dN_by_dx.contract(c, c_a).contract(dN_by_dx, c_b).shuffle(shuffle);
};

auto const Ks = [](auto const & dN_by_dx, Tetrahedron::Scalar const v, auto const & s)
{
	constexpr Tetrahedron::IndexPairs<1> s_a{{{1, 0}}};
	constexpr Tetrahedron::IndexPairs<1> s_b{{{1, 1}}};
	constexpr Tetrahedron::IndexPairs<0> s_I{};
	using Shuffle = Eigen::array<Eigen::Index, 4>;
	constexpr Shuffle shuffle{2, 0, 3, 1};
	// Ks = v * dN_a/dx_k * sigma_kl * dN_b/dx_l * delta_ij
	return v * dN_by_dx.contract(s, s_a).contract(dN_by_dx, s_b).contract(I<>, s_I).shuffle(shuffle);
};
}


namespace FeltElements
{
Tetrahedron::StiffnessAndForces Tetrahedron::KT(
	Node::Positions const & x, ShapeDerivativeTensor const & dN_by_dX,
	Scalar const lambda, Scalar const mu)
{
	Scalar const v = Tetrahedron::V(x);

	GradientTensor const & F = ex::dx_by_dX(x, dN_by_dX);  // non-expression for calc'ing Jacobian
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

	return StiffnessAndForces(K, T);
}

Tetrahedron::StiffnessTensor Tetrahedron::Kc(
	ShapeDerivativeTensor const & dN_by_dx, Scalar const v, ElasticityTensor const & c)
{
	return ex::Kc(dN_by_dx, v, c);
}

Tetrahedron::StiffnessTensor Tetrahedron::Ks(
	ShapeDerivativeTensor const & dN_by_dx, Scalar const v, StressTensor const & s)
{
	return ex::Ks(dN_by_dx, v, s);
}

Tetrahedron::ElasticityTensor Tetrahedron::c(Scalar J, Scalar lambda, Scalar mu)
{
	return ex::c(J, lambda, mu);
}

Tetrahedron::Node::Forces Tetrahedron::T(
	ShapeDerivativeTensor const & dN_by_dx, Scalar const v, StressTensor const & sigma)
{
	return ex::T(dN_by_dx, v, sigma);
}

Tetrahedron::StressTensor Tetrahedron::sigma(
	Tetrahedron::Scalar const J, Tetrahedron::GradientTensor const & b,
	Tetrahedron::Scalar const lambda, Tetrahedron::Scalar const mu)
{
	return ex::sigma(J, b, lambda, mu);
}

Tetrahedron::Scalar Tetrahedron::J(Tetrahedron::GradientTensor const & dx_by_dX)
{
	return to_matrix(dx_by_dX).determinant();
}

Tetrahedron::GradientTensor Tetrahedron::b(Tetrahedron::GradientTensor const & F)
{
	return ex::b(F);
}

Tetrahedron::GradientTensor Tetrahedron::dx_by_dX(
	Node::Positions const & x, ShapeDerivativeTensor const & dN_by_dX)
{
	return ex::dx_by_dX(x, dN_by_dX);
}

Tetrahedron::GradientTensor Tetrahedron::dx_by_dX(
	GradientTensor const & dx_by_dL, GradientTensor const & dL_by_dX)
{
	constexpr IndexPairs<1> x_X{{{0, 1}}};
	return dx_by_dL.contract(dL_by_dX, x_X);
}

Tetrahedron::ShapeDerivativeTensor Tetrahedron::dN_by_dX(Node::Positions const & X)
{
	auto const & dX_by_dL = ex::dX_by_dL(X);
	auto const & dL_by_dX = ex::dL_by_dX(dX_by_dL);
	return ex::dN_by_dX(dL_by_dX);
}

Tetrahedron::CartesianDerivativeTensor Tetrahedron::dx_by_dN(
	ShapeCartesianTransform const & N_to_x)
{
	Matrix<3, 4> dx_by_dN = to_matrix(N_to_x).block<3, 4>(1, 0);
	return to_tensor(dx_by_dN);
}

Tetrahedron::GradientTensor Tetrahedron::dL_by_dX(GradientTensor const & dX_by_dL)
{
	return ex::dL_by_dX(dX_by_dL);
}

Tetrahedron::GradientTensor Tetrahedron::dX_by_dL(Node::Positions const & X)
{
	return ex::dX_by_dL(X);
}

Tetrahedron::ShapeDerivativeTensor Tetrahedron::dN_by_dX(ShapeCartesianTransform const & N_to_x)
{
	// Interpolation: (1, x, y, z)^T = N_to_x * N, where N is 4x natural coordinates (corners).
	// Invert then strip constant terms, leaving just coefficients, i.e. the derivative.
	Matrix<4, 3> dN_by_dX = to_matrix(N_to_x).inverse().block<4, 3>(0, 1);
	return to_tensor(dN_by_dX);
}

Tetrahedron::ShapeCartesianTransform Tetrahedron::N_to_x(Node::Positions const & X)
{
	using Shuffle = Eigen::array<Eigen::Index, 2>;
	using Padding = Eigen::array<std::pair<Eigen::Index, Eigen::Index>, 2>;
	constexpr Shuffle transpose{{1, 0}};
	constexpr Padding padding{{{1, 0}, {0, 0}}};
	return X.shuffle(transpose).pad(padding, 1.0);
}

Tetrahedron::Scalar Tetrahedron::V(Node::Positions const & x)
{
	using Indices = Eigen::array<Eigen::Index, 2>;
	auto const & start_3x3 = x.slice(Indices{0, 0}, Indices{3, 3});
	auto const & end_1x3 = x.slice(Indices{3, 0}, Indices{1, 3});
	auto const & end_3x3 = end_1x3.broadcast(Indices{3, 1});

	MatrixTensor<3, 3> const delta = start_3x3 - end_3x3;
	auto const & mat_delta = to_matrix(delta);
	return std::abs(mat_delta.determinant() / 6.0);
}

Tetrahedron::Node::SpatialCoordProp Tetrahedron::x(Mesh & mesh)
{
	Node::SpatialCoordProp x_prop =
		mesh.request_vertex_property<OpenVolumeMesh::Vec3d>("x");
	x_prop->set_persistent(true);

	for (
		OpenVolumeMesh::VertexIter it_vtx = mesh.vertices_begin(); it_vtx != mesh.vertices_end();
		it_vtx++)
		x_prop[*it_vtx] = mesh.vertex(*it_vtx);

	return x_prop;
}

Tetrahedron::Node::Positions Tetrahedron::x(
	Mesh const & mesh, Node::Vtxhs const & vtxhs, Node::SpatialCoordProp const & x_prop)
{
	Node::Positions p;
	std::size_t node_idx = 0;
	for (Node::Vtxh const & vtxh : vtxhs)
	{
		OpenVolumeMesh::Vec3d const & vtx = x_prop[vtxh];
		p.chip(node_idx++, 0) = to_tensor(vtx);
	}
	return p;
}

Tetrahedron::Node::Positions Tetrahedron::X(Mesh const & mesh, Node::Vtxhs const & vtxhs)
{
	Node::Positions p;
	std::size_t node_idx = 0;
	for (Node::Vtxh const & vtxh : vtxhs)
	{
		OpenVolumeMesh::Vec3d const & vtx = mesh.vertex(vtxh);
		p.chip(node_idx++, 0) = to_tensor(vtx);
	}
	return p;
}

Tetrahedron::Node::Vtxhs Tetrahedron::vtxhs(
	OpenVolumeMesh::GeometricTetrahedralMeshV3d const & mesh,
	OpenVolumeMesh::CellHandle const & cellh)
{
	return mesh.get_cell_vertices(cellh);
}
} // namespace FeltElements
