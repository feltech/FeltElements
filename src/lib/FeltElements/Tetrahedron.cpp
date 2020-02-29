//
// Created by dave on 07/11/2019.
//
#include "TetGenIO.hpp"
#include "Tetrahedron.hpp"

namespace
{

template <typename Index>
auto delta(Index const i, Index const j)
{
	return i == j;
};

template <typename Matrix, typename ...Dim>
auto to_tensor(Matrix const & matrix, Dim const... dims)
{
	constexpr int rank = sizeof... (dims);
	return Eigen::TensorMap<Eigen::Tensor<FeltElements::Tetrahedron::Scalar const, rank>>{
		matrix.data(), {dims...}
	};
}

template <typename Tensor, Eigen::Index dim = 3>
auto to_matrix(Tensor const & tensor)
{
	return Eigen::Map<Eigen::Matrix3d const>{tensor.data(), dim, dim};
}
}

namespace FeltElements
{

Tetrahedron::Tetrahedron(FeltElements::TetGenIO const & mesh, std::size_t const tet_idx)
	:
	m_vertices{
		Node::PosMap{mesh.tet_corner_point(tet_idx, 0)},
		Node::PosMap{mesh.tet_corner_point(tet_idx, 1)},
		Node::PosMap{mesh.tet_corner_point(tet_idx, 2)},
		Node::PosMap{mesh.tet_corner_point(tet_idx, 3)},
	},
	m_displacements{
		Node::PosMap{mesh.tet_corner_displacement(tet_idx, 0)},
		Node::PosMap{mesh.tet_corner_displacement(tet_idx, 1)},
		Node::PosMap{mesh.tet_corner_displacement(tet_idx, 2)},
		Node::PosMap{mesh.tet_corner_displacement(tet_idx, 3)},
	},
	m_material_volume{
		std::abs(
			(m_vertices[1] - m_vertices[3])
			.cross(m_vertices[2] - m_vertices[3])
			.dot(m_vertices[0] - m_vertices[3])
		) / 6.0
	}
{}


void init(Tetrahedron::Mesh const & mesh, Tetrahedron::CellHandle const & tet_h)
{

}


Tetrahedron::IsoCoordDerivativeMatrix const Tetrahedron::dL_by_dN = // NOLINT(cert-err58-cpp)
	(Eigen::Matrix4d{} <<
		// clang-format off
		1, 1, 1, 1,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1).finished().block<3, 4>(1, 0);
		// clang-format on

Tetrahedron::ShapeDerivativeMatrix const Tetrahedron::dN_by_dL = // NOLINT(cert-err58-cpp)
	(Eigen::Matrix4d{} <<
		// clang-format off
		1, 1, 1, 1,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1).finished().inverse().block<4, 3>(0, 1);
		// clang-format on

Tetrahedron::ShapeDerivativeMatrix Tetrahedron::dN_by_dX() const
{
	Eigen::Matrix4d dX_by_dN;
	// Interpolation: (1, x, y, z)^T = dX_by_dN * N, where N is 4x natural coordinates (corners).
	// clang-format off
	dX_by_dN <<
		1,	1,	1,	1,
		X(0), X(1), X(2), X(3);
	// clang-format on

	// dN_i/dX_j, where i is row, j is column
	return dX_by_dN.inverse().block<4, 3>(0, 1);
}

Tetrahedron::IsoCoordDerivativeMatrix Tetrahedron::dx_by_dN() const
{
	return (IsoCoordDerivativeMatrix() << x(0), x(1), x(2), x(3)).finished();
}

Tetrahedron::Node::Pos Tetrahedron::x(Tetrahedron::Node::Index const idx) const
{
	return X(idx) + u(idx);
}

Tetrahedron::Node::Positions Tetrahedron::X(Mesh const & mesh, CellHandle const cellh)
{
	Node::Positions p;
	std::size_t node_idx = 0;
	for (OpenVolumeMesh::CellVertexIter iter = mesh.cv_iter(cellh); iter.valid(); iter++)
	{
		OpenVolumeMesh::Vec3d const & vtx = mesh.vertex(*iter);
		p.chip(node_idx++, 0) = to_tensor(vtx, 3);
	}
	return p;
}

Tetrahedron::Node::PosMap const & Tetrahedron::X(Node::Index const idx) const
{
	return m_vertices[idx];
}

Tetrahedron::Node::Positions Tetrahedron::u(
	Mesh const & mesh, Node::PosProperty const & displacements, CellHandle const cellh)
{
	Node::Positions p;
	std::size_t node_idx = 0;
	for (OpenVolumeMesh::CellVertexIter iter = mesh.cv_iter(cellh); iter.valid(); iter++)
	{
		Node::Pos const & vtx = displacements[*iter];
		p.chip(node_idx++, 0) = to_tensor(vtx, 3);
	}
	return p;
}

Tetrahedron::Node::PosMap & Tetrahedron::u(Tetrahedron::Node::Index const idx)
{
	return m_displacements[idx];
}

Tetrahedron::Node::PosMap const & Tetrahedron::u(Tetrahedron::Node::Index const idx) const
{
	return m_displacements[idx];
}

Tetrahedron::GradientMatrix Tetrahedron::dx_by_dX(
	IsoCoordDerivativeMatrix const & dx_by_dN, ShapeDerivativeMatrix const & dN_by_dX)
{
	return dx_by_dN * dN_by_dX;
}

Tetrahedron::Scalar Tetrahedron::J(Tetrahedron::GradientMatrix const & F)
{
	return F.determinant();
}

Tetrahedron::GradientMatrix Tetrahedron::b(Tetrahedron::GradientMatrix const & F)
{
	return F * F.transpose();
}

Tetrahedron::ShapeDerivativeMatrix
Tetrahedron::dN_by_dx(ShapeDerivativeMatrix const & dN_by_dX, GradientMatrix const & F)
{
	return dN_by_dX * F.inverse();
}

Tetrahedron::Scalar Tetrahedron::V() const
{
	return m_material_volume;
}

Tetrahedron::Scalar Tetrahedron::v() const
{
	Node::Positions xs;
	xs.chip(0, 0) = to_tensor(Node::Pos{x(0)}, 3);
	xs.chip(1, 0) = to_tensor(Node::Pos{x(1)}, 3);
	xs.chip(2, 0) = to_tensor(Node::Pos{x(2)}, 3);
	xs.chip(3, 0) = to_tensor(Node::Pos{x(3)}, 3);

	using Indices = Eigen::array<Eigen::Index, 2>;
	auto const & start_3x3 = xs.slice(Indices{0, 0}, Indices{3, 3});
	auto const & end_1x3 = xs.slice(Indices{3, 0}, Indices{1, 3});
	auto const & end_3x3 = end_1x3.broadcast(Indices{3, 1});

	MatrixTensor<> const delta = start_3x3 - end_3x3;
	auto const & mat_delta = to_matrix(delta);
	return std::abs(mat_delta.determinant()) / 6.0;

//	return std::abs(
//			(x(1) - x(3))
//				.cross(x(2) - x(3))
//				.dot(x(0) - x(3))
//		) / 6.0;
}

Tetrahedron::ElasticityTensor Tetrahedron::neo_hookian_elasticity(
	Tetrahedron::Scalar const J, Tetrahedron::Scalar const lambda, Tetrahedron::Scalar const mu)
{
	Tetrahedron::Scalar const lambda_prime = lambda / J;
	Tetrahedron::Scalar const mu_prime = (mu - lambda * std::log(J));

	std::size_t constexpr N = 3;

	static const ElasticityTensor c_lambda = ([]() {
		std::size_t constexpr N = 3;
		ElasticityTensor c;
		for (std::size_t i = 0; i < N; i++)
			for (std::size_t j = 0; j < N; j++)
				for (std::size_t k = 0; k < N; k++)
					for (std::size_t l = 0; l < N; l++)
						c(i, j, k, l) = delta(i, j) * delta(k, l);

		return c;
	}());

	static const ElasticityTensor c_mu = ([]() {
		ElasticityTensor c;
		for (std::size_t i = 0; i < N; i++)
			for (std::size_t j = 0; j < N; j++)
				for (std::size_t k = 0; k < N; k++)
					for (std::size_t l = 0; l < N; l++)
						c(i, j, k, l) = delta(i, k) * delta(j, l) + delta(i, l) * delta(j, k);

		return c;
	}());

	return lambda_prime * c_lambda + mu_prime * c_mu;
}

Tetrahedron::StressMatrix Tetrahedron::neo_hookian_stress(
	Tetrahedron::Scalar const J, Tetrahedron::GradientMatrix const & b,
	Tetrahedron::Scalar const lambda, Tetrahedron::Scalar const mu)
{
	static auto const I = StressMatrix::Identity();
	return (mu / J) * (b - I) + (lambda / J) * log(J) * I;
}

Tetrahedron::StiffnessMatrix Tetrahedron::Kcab(
	Tetrahedron::ShapeDerivativeMatrix const & dN_by_dx,
	Tetrahedron::Scalar const v,
	Tetrahedron::ElasticityTensor const & c,
	Tetrahedron::Node::Index const a,
	Tetrahedron::Node::Index const b)
{
	constexpr Eigen::array<Eigen::IndexPair<int>, 1> c_a = {
			Eigen::IndexPair<int>(1, 0)
	};
	constexpr Eigen::array<Eigen::IndexPair<int>, 1> c_b = {
		Eigen::IndexPair<int>(3, 0)
	};

	auto const & dNa_by_dx = to_tensor(dN_by_dx.row(a), 3);
	auto const & dNb_by_dx = to_tensor(dN_by_dx.row(b), 3);

	using Tensor2 = Eigen::TensorFixedSize<Tetrahedron::Scalar, Eigen::Sizes<3, 3>>;

	Tensor2 K = v * c.contract(dNb_by_dx, c_b).contract(dNa_by_dx, c_a);

	return to_matrix(K);
}

Tetrahedron::StiffnessMatrix Tetrahedron::Ksab(
	Tetrahedron::ShapeDerivativeMatrix const & dN_by_dx,
	Tetrahedron::Scalar const v,
	Tetrahedron::StressMatrix const & sigma,
	Tetrahedron::Node::Index const a,
	Tetrahedron::Node::Index const b)
{
	static auto const I = GradientMatrix::Identity();
	auto const & dNa_by_dxT = dN_by_dx.row(a);
	auto const & dNb_by_dx = dN_by_dx.row(b).transpose();

	return v * dNa_by_dxT * sigma * dNb_by_dx * I;
}
Tetrahedron::Node::Pos Tetrahedron::null_pos{0,0,0};
} // namespace FeltElements
