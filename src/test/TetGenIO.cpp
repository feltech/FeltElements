#include <FeltElements/TetGenIO.hpp>
#include <FeltElements/Typedefs.hpp>
#include <catch2/catch.hpp>

SCENARIO("Loading a tetrahedralisation")
{
	using namespace FeltElements;
	auto expected_counts = [](const TetGenIO & mesh,
							  std::size_t num_points,
							  std::size_t num_simplexes,
							  std::size_t num_corners) {
		CHECK(mesh.num_points() == num_points);
		CHECK(mesh.num_simplexes() == num_simplexes);
		CHECK(mesh.num_corners() == num_corners);
	};

	GIVEN("single tetrahedron mesh is loaded")
	{
		auto const & mesh = TetGenIO{"resources/tetgen/one/tet.1"};
		THEN("expected counts are reported")
		{
			expected_counts(mesh, 4, 1, 4);
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
	}

	GIVEN("two-element mesh is loaded")
	{
		auto const & io = TetGenIO{"resources/tetgen/two/tet.1"};
		FeltElements::Mesh ovm = io.to_mesh();

		WHEN("an OpenVolumeMesh is constructed")
		{
			THEN("mesh has correct number of elements")
			{
				CHECK(ovm.n_cells() == 2);
				CHECK(ovm.n_faces() == 7);
				CHECK(ovm.n_vertices() == 5);
			}
			THEN("expected node positions are reported")
			{
				auto const & vtxs1 = ovm.get_cell_vertices(OpenVolumeMesh::CellHandle{0});
				CHECK(ovm.vertex(vtxs1[0]) == OpenVolumeMesh::Vec3d{0, 0, 0});
				CHECK(ovm.vertex(vtxs1[1]) == OpenVolumeMesh::Vec3d{0, 0, 1});
				CHECK(ovm.vertex(vtxs1[2]) == OpenVolumeMesh::Vec3d{1, 0, 0});
				CHECK(ovm.vertex(vtxs1[3]) == OpenVolumeMesh::Vec3d{0, 0.5, 0.5});

				auto const & vtxs2 = ovm.get_cell_vertices(OpenVolumeMesh::CellHandle{1});
				CHECK(ovm.vertex(vtxs2[0]) == OpenVolumeMesh::Vec3d{0, 1, 0});
				CHECK(ovm.vertex(vtxs2[1]) == OpenVolumeMesh::Vec3d{0, 0, 0});
				CHECK(ovm.vertex(vtxs2[2]) == OpenVolumeMesh::Vec3d{1, 0, 0});
				CHECK(ovm.vertex(vtxs2[3]) == OpenVolumeMesh::Vec3d{0, 0.5, 0.5});
			}
		}
	}
}
