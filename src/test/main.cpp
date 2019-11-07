#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <FeltElements/MeshFile.hpp>

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
		std::string file_name = "resources/single/tet.1";
		WHEN("mesh is loaded")
		{
			auto const & mesh = FeltElements::MeshFile{file_name};
			THEN("expected counts are reported")
			{
				expected_counts(mesh, 1, 4, 4);
			}
		}
	}

	GIVEN("double quadratic tetrahedron mesh")
	{
		std::string file_name = "resources/double/tet.1";
		WHEN("mesh is loaded")
		{
			auto const & mesh = FeltElements::MeshFile{file_name};
			THEN("expected counts are reported")
			{
				expected_counts(mesh, 2, 10, 6);
			}
		}
	}
}