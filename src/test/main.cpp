#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <FeltElements/MeshFile.hpp>

SCENARIO("Loading a tetrahedralisation")
{
	GIVEN("single tetrahedron mesh")
	{
		std::string file_name = "resources/single/tet.1";
		WHEN("mesh is loaded")
		{
			auto const& mesh = FeltElements::MeshFile{file_name};
			THEN("expected number of simplexes are reported")
			{
				CHECK(mesh.num_simplexes() == 1);
			}
		}
	}

	GIVEN("double tetrahedron mesh")
	{
		std::string file_name = "resources/double/tet.1";
		WHEN("mesh is loaded")
		{
			auto const& mesh = FeltElements::MeshFile{file_name};
			THEN("expected number of simplexes are reported")
			{
				CHECK(mesh.num_simplexes() == 2);
			}
		}
	}
}