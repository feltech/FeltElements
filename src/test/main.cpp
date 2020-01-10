#define EIGEN_DEFAULT_IO_FORMAT Eigen::IOFormat(3, 0, "  ", "\n", "(", ")", "", "")
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <FeltElements/MeshFile.hpp>
#include <FeltElements/Tetrahedron.hpp>

SCENARIO("Loading a tetrahedralisation")
{
	auto expected_counts = [](const FeltElements::MeshFile & mesh,
							  auto num_simplexes,
							  auto num_corners,
							  auto num_trifaces) {
		CHECK(mesh.num_simplexes() == num_simplexes);
		CHECK(mesh.num_corners() == num_corners);
		CHECK(mesh.num_trifaces() == num_trifaces);
	};

	GIVEN("single tetrahedron mesh")
	{
		std::string const file_name = "resources/single/tet.1";
		WHEN("mesh is loaded")
		{
			auto const & mesh = FeltElements::MeshFile{file_name};
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

				THEN("expected volume is reported")
				{
					CHECK(tet.V() == 1.0/6);
				}
			}
		}
	}

	GIVEN("double quadratic tetrahedron mesh")
	{
		std::string const file_name = "resources/double/tet.1";
		WHEN("mesh is loaded")
		{
			auto const & mesh = FeltElements::MeshFile{file_name};
			THEN("expected counts are reported")
			{
				expected_counts(mesh, 2, 10, 6);
			}

			AND_WHEN("tetrahedrons are loaded from mesh")
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
}

SCENARIO("Simple deformation gradient")
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
				auto const & F = tet.dx_by_dX(dN_by_dX);

				THEN("deformation gradient is the identity matrix")
				{
					CHECK(F == Eigen::Matrix3d::Identity());
				}

				AND_WHEN("derivative of shape function wrt spatial coords is calculated")
				{
					auto const & dN_by_dx = tet.dN_by_dx(dN_by_dX);

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
					auto const & F = tet.dx_by_dX(dN_by_dX);

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

					AND_WHEN("derivative of shape function wrt spatial coords is calculated")
					{
						auto const & dN_by_dx = tet.dN_by_dx(dN_by_dX);

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
		}
	}
}