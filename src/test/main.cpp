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
					CHECK(tet.vertex(0) == Eigen::Vector3d{0, 0, 0});
					CHECK(tet.vertex(1) == Eigen::Vector3d{0, 1, 0});
					CHECK(tet.vertex(2) == Eigen::Vector3d{0, 0, 1});
					CHECK(tet.vertex(3) == Eigen::Vector3d{1, 0, 0});
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
					CHECK(tet1.vertex(0) == Eigen::Vector3d{0, 0, 0});
					CHECK(tet1.vertex(1) == Eigen::Vector3d{0, 0, 1});
					CHECK(tet1.vertex(2) == Eigen::Vector3d{1, 0, 0});
					CHECK(tet2.vertex(3) == Eigen::Vector3d{0, 0.5, 0.5});

					CHECK(tet2.vertex(0) == Eigen::Vector3d{0, 1, 0});
					CHECK(tet2.vertex(1) == Eigen::Vector3d{0, 0, 0});
					CHECK(tet2.vertex(2) == Eigen::Vector3d{1, 0, 0});
					CHECK(tet2.vertex(3) == Eigen::Vector3d{0, 0.5, 0.5});
				}
			}
		}
	}
}

SCENARIO("Simple deformation gradient")
{
	GIVEN("simple tetrahedron")
	{
		std::string const file_name = "resources/single/tet.1";
		auto const & mesh = FeltElements::MeshFile{file_name};
		auto tet = FeltElements::Tetrahedron{mesh, 0};

		THEN("deformation gradient is initially zero")
		{
			CHECK(tet.deformation_gradient() == Eigen::Matrix3d::Identity());
		}

		WHEN("a node is deformed")
		{
			tet.displacement(0) += Eigen::Vector3d{0.5, 0, 0};

			THEN("deformation gradient is correct")
			{
				Eigen::Matrix3d expected;
				// clang-format off
				expected << //NOLINT
					0.5, -0.5, -0.5,
					0, 1, 0,
					0, 0, 1;
				// clang-format on
				CHECK(tet.deformation_gradient() == expected);
			}
		}
	}
}