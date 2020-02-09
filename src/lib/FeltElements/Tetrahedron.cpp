//
// Created by dave on 07/11/2019.
//
#include <boost/range/irange.hpp>

#include "Tetrahedron.hpp"
#include "MeshFile.hpp"

namespace
{
	auto to_tensor(auto const & matrix, auto const... dims)
	{
		using MatrixType = decltype(matrix);
		constexpr int rank = sizeof... (dims);
		return Eigen::TensorMap<Eigen::Tensor<double const, rank>>{
			matrix.data(), {dims...}
		};
	}

	auto from_tensor(auto const & tensor)
	{
		return Eigen::Map<Eigen::Matrix3d const>{tensor.data(), 3, 3};
	}
}

namespace FeltElements
{
Tetrahedron::Tetrahedron(FeltElements::MeshFile const & mesh, std::size_t const tet_idx)
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
	m_volume{
		std::abs(
			(m_vertices[1] - m_vertices[3])
			.cross(m_vertices[2] - m_vertices[3])
			.dot(m_vertices[0] - m_vertices[3])
		) / 6.0
	}
{}

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

Tetrahedron::ShapeDerivativeMatrix Tetrahedron::dN_by_dX()
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

Tetrahedron::Node::PosMap const & Tetrahedron::X(Node::Index const idx) const
{
	return m_vertices[idx];
}

Tetrahedron::Node::PosMap & Tetrahedron::u(Tetrahedron::Node::Index const idx)
{
	return m_displacements[idx];
}

Tetrahedron::Node::PosMap const & Tetrahedron::u(Tetrahedron::Node::Index const idx) const
{
	return m_displacements[idx];
}

Tetrahedron::GradientMatrix Tetrahedron::dx_by_dX(ShapeDerivativeMatrix const & dN_by_dX) const
{
	return dx_by_dN() * dN_by_dX;
}

Tetrahedron::ShapeDerivativeMatrix
Tetrahedron::dN_by_dx(ShapeDerivativeMatrix const & dN_by_dX)
{
	return dN_by_dX * dx_by_dX(dN_by_dX).inverse();
}

double Tetrahedron::V() const
{
	return m_volume;
}

Tetrahedron::ElasticityTensor Tetrahedron::neo_hookian_elasticity(
	Tetrahedron::GradientMatrix const & F, double const lambda, double const mu)
{
	constexpr std::size_t N = 3;
	Tetrahedron::ElasticityTensor c;
	double const J = F.determinant();
	double const lambda_prime = lambda / J;
	double const mu_prime = (mu - lambda * std::log(J));

	for (std::size_t i : boost::irange(N))
		for (std::size_t j : boost::irange(N))
			for (std::size_t k : boost::irange(N))
				for (std::size_t l : boost::irange(N))
				{
					c(i, j, k, l) = lambda_prime * delta(i, j) * delta(k, l) +
						mu_prime * (delta(i, k) * delta(j, l) + delta(i, l) * delta(j, k));
				}
	return c;
}

Tetrahedron::GradientMatrix Tetrahedron::Kcab(
	Tetrahedron::ShapeDerivativeMatrix const & dN_by_dx,
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

	Eigen::Tensor<double, 1> const & dNa_by_dx = to_tensor(dN_by_dx.row(a), 3);
	Eigen::Tensor<double, 1> const & dNb_by_dx = to_tensor(dN_by_dx.row(b), 3);

	Eigen::Tensor<double, 2> K = c.contract(dNb_by_dx, c_b).contract(dNa_by_dx, c_a);

	return from_tensor(K);
}

} // namespace FeltElements
