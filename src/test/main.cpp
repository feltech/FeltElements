#define EIGEN_DEFAULT_IO_FORMAT Eigen::IOFormat(3, 0, ", ", ",\n", "", "", "", "")
#define CATCH_CONFIG_MAIN
#include <unsupported/Eigen/CXX11/Tensor>

#include <FeltElements/MeshFile.hpp>
#include <FeltElements/Tetrahedron.hpp>
#include "util.hpp"
#include <catch2/catch.hpp>

std::string const file_name_single = "resources/single/tet.1";
std::string const file_name_double = "resources/double/tet.1";

using namespace FeltElements;


SCENARIO("Loading a tetrahedralisation")
{
//	char cwd[500];
//	getcwd(cwd, 500);
//	std::cerr << "Executing tests in " << cwd << std::endl;

	auto expected_counts = [](const FeltElements::MeshFile & mesh,
							  auto num_simplexes,
							  auto num_corners,
							  auto num_trifaces) {
		CHECK(mesh.num_simplexes() == num_simplexes);
		CHECK(mesh.num_corners() == num_corners);
		CHECK(mesh.num_trifaces() == num_trifaces);
	};

	GIVEN("single tetrahedron mesh is loaded")
	{
		auto const & mesh = FeltElements::MeshFile{file_name_single};
		THEN("expected counts are reported")
		{
			expected_counts(mesh, 1, 4, 4);
		}

		THEN("expected node positions are reported")
		{
			CHECK(mesh.points()[3 * mesh.corners()[0]] == 0);
			CHECK(mesh.points()[3 * mesh.corners()[0] + 1] == 0);
			CHECK(mesh.points()[3 * mesh.corners()[0] + 2] == 0);

			CHECK(mesh.points()[3 * mesh.corners()[1]] == 0);
			CHECK(mesh.points()[3 * mesh.corners()[1] + 1] == 1);
			CHECK(mesh.points()[3 * mesh.corners()[1] + 2] == 0);

			CHECK(mesh.points()[3 * mesh.corners()[2]] == 0);
			CHECK(mesh.points()[3 * mesh.corners()[2] + 1] == 0);
			CHECK(mesh.points()[3 * mesh.corners()[2] + 2] == 1);

			CHECK(mesh.points()[3 * mesh.corners()[3]] == 1);
			CHECK(mesh.points()[3 * mesh.corners()[3] + 1] == 0);
			CHECK(mesh.points()[3 * mesh.corners()[3] + 2] == 0);
		}

		THEN("displacements are initially zero")
		{
			CHECK(mesh.displacements()[3 * mesh.corners()[0]] == 0);
			CHECK(mesh.displacements()[3 * mesh.corners()[0] + 1] == 0);
			CHECK(mesh.displacements()[3 * mesh.corners()[0] + 2] == 0);

			CHECK(mesh.displacements()[3 * mesh.corners()[1]] == 0);
			CHECK(mesh.displacements()[3 * mesh.corners()[1] + 1] == 0);
			CHECK(mesh.displacements()[3 * mesh.corners()[1] + 2] == 0);

			CHECK(mesh.displacements()[3 * mesh.corners()[2]] == 0);
			CHECK(mesh.displacements()[3 * mesh.corners()[2] + 1] == 0);
			CHECK(mesh.displacements()[3 * mesh.corners()[2] + 2] == 0);

			CHECK(mesh.displacements()[3 * mesh.corners()[3]] == 0);
			CHECK(mesh.displacements()[3 * mesh.corners()[3] + 1] == 0);
			CHECK(mesh.displacements()[3 * mesh.corners()[3] + 2] == 0);
		}

		AND_WHEN("tetrahedron is loaded from mesh")
		{
			auto tet = FeltElements::Tetrahedron{mesh, 0};
			THEN("expected node positions are reported")
			{
				CHECK(tet.X(0) == Eigen::Vector3d{0, 0, 0});
				CHECK(tet.X(1) == Eigen::Vector3d{0, 1, 0});
				CHECK(tet.X(2) == Eigen::Vector3d{0, 0, 1});
				CHECK(tet.X(3) == Eigen::Vector3d{1, 0, 0});
			}

			THEN("expected initial volume is reported")
			{
				CHECK(tet.V() == 1.0/6);
			}

			THEN("deformed volume matches initial volume")
			{
				CHECK(tet.v() == tet.V());
			}
		}
	}

	GIVEN("double quadratic tetrahedron mesh is loaded")
	{
		auto const & mesh = FeltElements::MeshFile{file_name_double};
		THEN("expected counts are reported")
		{
			expected_counts(mesh, 2, 10, 6);
		}

		WHEN("tetrahedrons are loaded from mesh")
		{
			auto tet1 = FeltElements::Tetrahedron{mesh, 0};
			auto tet2 = FeltElements::Tetrahedron{mesh, 1};

			THEN("expected node positions are reported")
			{
				CHECK(tet1.X(0) == Eigen::Vector3d{0, 0, 0});
				CHECK(tet1.X(1) == Eigen::Vector3d{0, 0, 1});
				CHECK(tet1.X(2) == Eigen::Vector3d{1, 0, 0});
				CHECK(tet2.X(3) == Eigen::Vector3d{0, 0.5, 0.5});

				CHECK(tet2.X(0) == Eigen::Vector3d{0, 1, 0});
				CHECK(tet2.X(1) == Eigen::Vector3d{0, 0, 0});
				CHECK(tet2.X(2) == Eigen::Vector3d{1, 0, 0});
				CHECK(tet2.X(3) == Eigen::Vector3d{0, 0.5, 0.5});
			}

			THEN("expected volumes are reported")
			{
				CHECK(tet1.V() == 1.0 / 12);
				CHECK(tet2.V() == 1.0 / 12);
			}
		}
	}
}


SCENARIO("Shape function derivatives")
{
	GIVEN("simple tetrahedron")
	{
//		std::string const file_name = "resources/single/tet.1";
		std::string const file_name = "resources/double/tet.1";
		auto const & mesh = FeltElements::MeshFile{file_name};
		auto tet = FeltElements::Tetrahedron{mesh, 0};

		std::stringstream ss;
		ss << "Tetrahedron vertices: ";
		for (std::size_t i = 0; i < 4; i++)
			ss << tet.X(i).transpose() << " ; ";
		ss << std::endl;
		INFO(ss.str());

		THEN("derivative of shape function wrt isoparametric coords is correct")
		{
			Eigen::Matrix<double, 4, 3> expected;
			// clang-format off
			expected << // NOLINT
				-1, -1, -1,
				1, 0, 0,
				0, 1, 0,
				0, 0, 1;
			// clang-format on
			CHECK(tet.dN_by_dL == expected);
		}

		THEN("derivative of isoparametric coords wrt shape function is correct")
		{
			Eigen::Matrix<double, 3, 4> expected;
			// clang-format off
			expected << // NOLINT
				0, 1, 0, 0,
				0, 0, 1, 0,
				0, 0, 0, 1;
			// clang-format on
			CHECK(tet.dL_by_dN == expected);
		}

		WHEN("derivative of shape function wrt material coords is calculated")
		{
			auto const & dN_by_dX = tet.dN_by_dX();

			THEN("derivative is correct")
			{
				Eigen::Matrix<double, 4, 3> expected;
				// clang-format off
				expected << // NOLINT
					-1, -1, -1,
					0, -1, 1,
					1, 0, 0,
					0, 2, 0;
				// clang-format on
				CHECK(dN_by_dX == expected);
			}

			THEN("derivative of isoparametric wrt material coords is correct")
			{
				Eigen::Matrix3d expected;
				// clang-format off
				expected << // NOLINT
					0, -1, 1,
					1, 0, 0,
					0, 2, 0;
				// clang-format on
				CHECK(tet.dL_by_dN * dN_by_dX == expected);
			}

			AND_WHEN("deformation gradient is calculated")
			{
				auto const & dx_by_dN = tet.dx_by_dN();
				auto const & F = Tetrahedron::dx_by_dX(dx_by_dN, dN_by_dX);

				THEN("deformation gradient is the identity matrix")
				{
					CHECK(F == Eigen::Matrix3d::Identity());
				}

				AND_WHEN("derivative of shape function wrt spatial coords is calculated")
				{
					auto const & dN_by_dx = Tetrahedron::dN_by_dx(dN_by_dX, F);

					THEN("derivative is the same as in material coords")
					{
						CHECK(dN_by_dx == dN_by_dX);
					}
				}

			}

			AND_WHEN("a node is deformed")
			{
				tet.u(0) += Eigen::Vector3d{0.5, 0, 0};

				THEN("material coords are unchanged")
				{
					CHECK(tet.X(0) == Eigen::Vector3d{0, 0, 0});
					CHECK(tet.X(1) == Eigen::Vector3d{0, 0, 1});
					CHECK(tet.X(2) == Eigen::Vector3d{1, 0, 0});
					CHECK(tet.X(3) == Eigen::Vector3d{0, 0.5, 0.5});
				}

				THEN("spatial coords equal material with deformation")
				{
					CHECK(tet.x(0) == Eigen::Vector3d{0.5, 0, 0});
					CHECK(tet.x(1) == Eigen::Vector3d{0, 0, 1});
					CHECK(tet.x(2) == Eigen::Vector3d{1, 0, 0});
					CHECK(tet.x(3) == Eigen::Vector3d{0, 0.5, 0.5});
				}

				THEN("spatial volume has changed")
				{
					CHECK(tet.v() == tet.V() / 2);
				}

				AND_WHEN("derivative of shape function wrt material coords is calculated")
				{
					auto const & dN_by_dX_deformed = tet.dN_by_dX();

					THEN("derivative is unchanged from before it was deformed")
					{
						CHECK(dN_by_dX_deformed == dN_by_dX);
					}
				}
				AND_WHEN("deformation gradient is calculated")
				{
					auto const & dx_by_dN = tet.dx_by_dN();
					auto const & F = tet.dx_by_dX(dx_by_dN, dN_by_dX);

					THEN("deformation gradient is correct")
					{
						Eigen::Matrix3d expected;
						// clang-format off
						expected << // NOLINT
							0.5, -0.5, -0.5,
							0, 1, 0,
							0, 0, 1;
						// clang-format on
						CHECK(F == expected);
					}

					AND_WHEN("Jacobian of deformation gradient is calculated")
					{
						double const J = Tetrahedron::J(F);

						THEN("value is correct")
						{
							CHECK(J == 0.5);
						}
					}

					AND_WHEN("Left Cauchy-Green / Finger tensor is calculated")
					{
						Tetrahedron::GradientMatrix const & b = Tetrahedron::b(F);
						INFO(b)

						THEN("value is correct")
						{
							Tetrahedron::GradientMatrix expected;
							// clang-format off
							expected << // NOLINT
								 0.75, -0.5, -0.5,
								-0.5,    1,    0,
								-0.5,    0,    1;
							// clang-format on
							CHECK(b.isApprox(expected));
						}
					}

					AND_WHEN("derivative of shape function wrt spatial coords is calculated")
					{
						auto const & dN_by_dx = tet.dN_by_dx(dN_by_dX, F);

						std::stringstream ss;
						ss << "dN_by_dX = "  << std::endl << dN_by_dX << std::endl;
						ss << "F^-1 = "  << std::endl << F.inverse() << std::endl;
						INFO(ss.str());

						THEN("derivative is correct")
						{
							Eigen::Matrix<double, 4, 3> expected;
							// clang-format off
							expected << // NOLINT
								-2, -2, -2,
								0, -1, 1,
								2, 1, 1,
								0, 2, 0;
							// clang-format on
							CHECK(dN_by_dx == expected);
						}
					}
				}
			}
		} // WHEN("derivative of shape function wrt material coords is calculated")

	} // End GIVEN("simple tetrahedron")
}

SCENARIO("Neo-hookian tangent stiffness matrix")
{
	GIVEN("material properties")
	{
		using namespace FeltElements;

		// Material properties: https://www.azom.com/properties.aspx?ArticleID=920
		double const mu = 0.4; // Shear modulus: 0.0003 - 0.02
		double const E = 1;	   // Young's modulus: 0.001 - 0.05
		// Lame's first parameter: https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
		double lambda = (mu * (E - 2 * mu)) / (3 * mu - E);

		std::stringstream s;
		s << "Lambda = " << lambda << "; mu = " << mu;
		INFO(s.str());

		AND_GIVEN("an undeformed tetrahedron")
		{
			auto const mesh = FeltElements::MeshFile{file_name_single};
			auto const tet = FeltElements::Tetrahedron{mesh, 0};
			auto const v = tet.v();
			auto const & dN_by_dX = tet.dN_by_dX();
			auto const & dx_by_dN = tet.dx_by_dN();
			auto const & F = Tetrahedron::dx_by_dX(dx_by_dN, dN_by_dX);
			auto const & dN_by_dx = Tetrahedron::dN_by_dx(dN_by_dX, F);
			auto const J = Tetrahedron::J(F);
			auto const & b = Tetrahedron::b(F);

			WHEN("neo-hookian elasticity tensor is calculated")
			{
				auto const & c = Tetrahedron::neo_hookian_elasticity(J, lambda, mu);

				THEN("it has expected values")
				{
					FeltElements::Tetrahedron::ElasticityTensor expected;
					expected.setValues({{{{1.2, 0, 0}, {0, 0.4, 0}, {0, 0, 0.4}},
											{{0, 0.4, 0}, {0.4, 0, 0}, {0, 0, 0}},
											{{0, 0, 0.4}, {0, 0, 0}, {0.4, 0, 0}}},
										{{{0, 0.4, 0}, {0.4, 0, 0}, {0, 0, 0}},
											{{0.4, 0, 0}, {0, 1.2, 0}, {0, 0, 0.4}},
											{{0, 0, 0}, {0, 0, 0.4}, {0, 0.4, 0}}},
										{{{0, 0, 0.4}, {0, 0, 0}, {0.4, 0, 0}},
											{{0, 0, 0}, {0, 0, 0.4}, {0, 0.4, 0}},
											{{0.4, 0, 0}, {0, 0.4, 0}, {0, 0, 1.2}}}});
					Eigen::Tensor<bool, 0> comparison = ((c - expected).abs() < 0.00001).all();

					INFO("Expected:")
					INFO(expected);
					INFO("Actual:")
					INFO(c);
					CHECK(comparison(0));
				}

				AND_WHEN("constitutive component of tangent matrix is calculated")
				{
					auto const & Kcab = Tetrahedron::Kcab(dN_by_dx, c, v, 0, 1);

					INFO("Kcab:")
					INFO(Kcab)

					THEN("matrix of values are as expected")
					{
						Eigen::Matrix3d expected;
						expected <<
						// clang-format off
							-0.0667, -0.0667,    0,
							-0.0667, -0.2, -0.0667,
							0, -0.0667, -0.0667;
						// clang-format on
						CHECK(Kcab.isApprox(expected, 0.0005));
					}
				}

			} // End WHEN("neo-hookian elasticity tensor is calculated")

			WHEN("neo-hookian Cauchy stress tensor is calculated")
			{
				Tetrahedron::GradientMatrix const & sigma = Tetrahedron::neo_hookian_stress(
					J, b, lambda, mu);

				INFO("sigma:")
				INFO(sigma)

				THEN("Cauchy stress tensor is zero")
				{
					Eigen::Matrix3d expected;
					expected.setZero();
					CHECK(sigma == expected);
				}

				AND_WHEN("initial stress component of tangent stiffness matrix is calculated")
				{
					Tetrahedron::GradientMatrix const & Ksab = Tetrahedron::Ksab(
						dN_by_dx, sigma, v, 0, 1);

					THEN("stress component is zero")
					{
						Eigen::Matrix3d expected;
						expected.setZero();
						CHECK(Ksab == expected);
					}
				}
			}
		} // End AND_GIVEN("an undeformed tetrahedron")

		AND_GIVEN("a deformed tetrahedron")
		{
			auto const mesh = FeltElements::MeshFile{file_name_single};
			auto tet = FeltElements::Tetrahedron{mesh, 0};
			// Do the deformation.
			tet.u(0) += Eigen::Vector3d{0.5, 0, 0};

			auto const v = tet.v();
			auto const & dN_by_dX = tet.dN_by_dX();
			auto const & dx_by_dN = tet.dx_by_dN();
			auto const & F = Tetrahedron::dx_by_dX(dx_by_dN, dN_by_dX);
			auto const & dN_by_dx = Tetrahedron::dN_by_dx(dN_by_dX, F);
			auto const J = Tetrahedron::J(F);
			auto const & b = Tetrahedron::b(F);

			WHEN("neo-hookian elasticity tensor is calculated")
			{
				auto const & c = Tetrahedron::neo_hookian_elasticity(J, lambda, mu);

				THEN("it has expected values")
				{
					FeltElements::Tetrahedron::ElasticityTensor expected;
					expected.setValues({{{{1.2, 0, 0}, {0, 0.4, 0}, {0, 0, 0.4}},
											{{0, 0.4, 0}, {0.4, 0, 0}, {0, 0, 0}},
											{{0, 0, 0.4}, {0, 0, 0}, {0.4, 0, 0}}},
										{{{0, 0.4, 0}, {0.4, 0, 0}, {0, 0, 0}},
											{{0.4, 0, 0}, {0, 1.2, 0}, {0, 0, 0.4}},
											{{0, 0, 0}, {0, 0, 0.4}, {0, 0.4, 0}}},
										{{{0, 0, 0.4}, {0, 0, 0}, {0.4, 0, 0}},
											{{0, 0, 0}, {0, 0, 0.4}, {0, 0.4, 0}},
											{{0.4, 0, 0}, {0, 0.4, 0}, {0, 0, 1.2}}}});
					Eigen::Tensor<bool, 0> comparison = ((c - expected).abs() < 0.00001).all();

					INFO("Expected:")
					INFO(expected);
					INFO("Actual:")
					INFO(c);
					CHECK(comparison(0));
				}

				AND_WHEN("constitutive component of tangent matrix is calculated")
				{
					auto const & Kcab = Tetrahedron::Kcab(dN_by_dx, c, v, 0, 1);

					INFO("Kcab:")
					INFO(Kcab)

					THEN("matrix of values are as expected")
					{
						Eigen::Matrix3d expected;
						expected <<
								 // clang-format off
								 -0.0667, -0.0667,    0,
							-0.0667, -0.2, -0.0667,
							0, -0.0667, -0.0667;
						// clang-format on
						CHECK(Kcab.isApprox(expected, 0.0005));
					}
				}

			} // End WHEN("neo-hookian elasticity tensor is calculated")

			WHEN("neo-hookian Cauchy stress tensor is calculated")
			{
				Tetrahedron::GradientMatrix const & sigma = Tetrahedron::neo_hookian_stress(
					J, b, lambda, mu);

				INFO("sigma:")
				INFO(sigma)

				THEN("Cauchy stress tensor is zero")
				{
					Eigen::Matrix3d expected;
					expected.setZero();
					CHECK(sigma == expected);
				}

				AND_WHEN("initial stress component of tangent stiffness matrix is calculated")
				{
					Tetrahedron::GradientMatrix const & Ksab = Tetrahedron::Ksab(
						dN_by_dx, sigma, v, 0, 1);

					THEN("stress component is zero")
					{
						Eigen::Matrix3d expected;
						expected.setZero();
						CHECK(Ksab == expected);
					}
				}
			}
		} // End AND_GIVEN("an undeformed tetrahedron")
	}
}
