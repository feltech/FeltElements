#define EIGEN_DEFAULT_IO_FORMAT Eigen::IOFormat(3, 0, ", ", ",\n", "", "", "", "")
#include <boost/range/irange.hpp>

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_CONSOLE_WIDTH 200
#include "common.hpp"

char const * const file_name_one = "resources/single/tet.1";
char const * const file_name_two = "resources/two/tet.1";
char const * const file_name_two_quadratic = "resources/double/tet.1";

using namespace FeltElements;

SCENARIO("Loading a tetrahedralisation")
{
	char cwd[500];
	getcwd(cwd, 500);
	std::cerr << "Executing tests in " << cwd << std::endl;

	auto expected_counts =
		[](const TetGenIO & mesh, std::size_t num_points, std::size_t num_simplexes,
				std::size_t num_corners) {
			CHECK(mesh.num_points() == num_points);
			CHECK(mesh.num_simplexes() == num_simplexes);
			CHECK(mesh.num_corners() == num_corners);
		};

	GIVEN("single tetrahedron mesh is loaded")
	{
		auto const & mesh = TetGenIO{file_name_one};
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
}

SCENARIO("OpenVolumeMesh construction")
{
	GIVEN("two-element mesh is loaded")
	{
		auto const & io = TetGenIO{file_name_two};

		WHEN("an OpenVolumeMesh is constructed")
		{
			auto ovm = io.to_mesh();

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
} // End SCENARIO("OpenVolumeMesh construction")

SCENARIO("Metrics of undeformed mesh")
{
	GIVEN("one-element mesh")
	{
		auto const & io = TetGenIO{file_name_one};
		auto mesh = io.to_mesh();


		WHEN("vertex index mapping is fetched")
		{
			auto const & vtxhs = Tetrahedron::vtxhs(mesh, 0);

			THEN("mapping is expected")
			{
				CHECK(vtxhs == Tetrahedron::Node::Vtxhs{0, 2, 3, 1});
			}

			AND_WHEN("material node position tensor is constructed")
			{
				Tetrahedron::Node::Positions X = Tetrahedron::X(mesh, vtxhs);

				THEN("expected positions are reported")
				{
					// clang-format off
					check_equal(X, "X", {
						{0, 0, 0},
						{0, 1, 0},
						{0, 0, 1},
						{1, 0, 0}
					});
					// clang-format on
				}

				AND_WHEN("material volume is calculated")
				{
					Tetrahedron::Scalar V = Tetrahedron::V(X);

					THEN("volume is correct")
					{
						CHECK(V == 1.0 / 6.0);
					}
				}

				AND_WHEN("spatial coordinate property is created")
				{
					using namespace OpenVolumeMesh;
					VertexPropertyT<Vec3d> const & x_prop = Tetrahedron::x(mesh);

					THEN("spatial coordinates equal material coordinates")
					{
						std::vector<Vec3d> all_X{};
						std::vector<Vec3d> all_x{};

						for (VertexIter it_vtx = mesh.vertices_begin();
							 it_vtx != mesh.vertices_end();
							 it_vtx++)
						{
							all_X.push_back(mesh.vertex(*it_vtx));
							all_x.push_back(x_prop[*it_vtx]);
						}

						CHECK(all_X == all_x);
					}

					AND_WHEN("spatial node position tensor is constructed")
					{
						Tetrahedron::Node::Positions x = Tetrahedron::x(vtxhs, x_prop);

						THEN("spatial node positions equal material positions")
						{
							// clang-format off
							check_equal(x, "x", X, "X");
							// clang-format on
						}
					}
				}
			}
		}
	}
}

SCENARIO("Metrics of deformed mesh")
{
	GIVEN("one-element mesh")
	{
		auto [X, x] = load_tet(file_name_one);

		INFO("Tetrahedron vertices:")
		INFO(X)

		WHEN("a node is deformed")
		{
			x(0, 0) += 0.5;

			THEN("spatial position is updated")
			{
				// clang-format off
				check_equal(x, "x", {
					{0.5, 0, 0},
					{0, 1, 0},
					{0, 0, 1},
					{1, 0, 0}
				});
				// clang-format on
			}

			AND_WHEN("spatial volume is calculated")
			{
				Tetrahedron::Scalar v = Tetrahedron::V(x);

				THEN("volume is correct")
				{
					CHECK(v == 1.0 / 12.0);
				}
			}
		}
	}
}

SCENARIO("Coordinate derivatives in undeformed mesh")
{
	THEN("derivative of shape wrt local coords is correct")
	{
		check_equal(
			Attribute::Element::MaterialShapeDerivative::dN_by_dL,
			"dN_by_dL",
			{
				// clang-format off
				{-1, -1, -1},
				{1, 0, 0},
				{0, 1, 0},
				{0, 0, 1}
				// clang-format on
			});
	}

	THEN("derivative of local wrt shape coords is correct")
	{
		check_equal(
			Attribute::Element::MaterialShapeDerivative::dL_by_dN,
			"dL_by_dN",
			{
				// clang-format off
				{0, 1, 0, 0},
				{0, 0, 1, 0},
				{0, 0, 0, 1}
				// clang-format on
			});
	}

	THEN("derivative of local/shape wrt shape/local coords are inverses")
	{
		constexpr Tetrahedron::IndexPairs<1> LN_NL{{{1, 0}}};
		Tetrahedron::MatrixTensor<3> const identity =
			Attribute::Element::MaterialShapeDerivative::dL_by_dN.contract(
				Attribute::Element::MaterialShapeDerivative::dN_by_dL, LN_NL);

		check_equal(
			identity,
			"dL_by_dN * dN_by_dL",
			{
				// clang-format off
			{1, 0, 0},
			{0, 1, 0},
			{0, 0, 1}
				// clang-format on
			});
	}

	GIVEN("a one-element mesh")
	{
		auto [X, x] = load_tet(file_name_one);
		(void)x;

		AND_WHEN("derivative of material wrt local coords is calculated")
		{
			auto const & dX_by_dL = Tetrahedron::dX_by_dL(X);

			THEN("derivative is correct")
			{
				// clang-format off
				check_equal(dX_by_dL, "dX_by_dL", {
					{0, 0, 1},
					{1, 0, 0},
					{0, 1, 0}
				});
				// clang-format on
			}

			AND_WHEN("derivative of local wrt material coords is calculated")
			{
				auto const & dL_by_dX = Tetrahedron::dL_by_dX(dX_by_dL);

				THEN("derivative is correct")
				{
					// clang-format off
					check_equal(dL_by_dX, "dL_by_dX", {
						{0, 1, 0},
						{0, 0, 1},
						{1, 0, 0}
					});
					// clang-format on
				}

				AND_WHEN("derivative of material wrt shape coords is calculated")
				{
					auto const & dN_by_dX = Tetrahedron::dN_by_dX(X);

					THEN("derivative is correct")
					{
						// clang-format off
						check_equal(dN_by_dX, "dN_by_dX", {
							{-1, -1, -1},
							{0, 1, 0},
							{0, 0, 1},
							{1, 0, 0}
						});
						// clang-format on
					}
				}
			}
		}

		AND_WHEN("transformation from natural to cartesian coordinates is calculated")
		{
			auto const & N_to_x = Tetrahedron::N_to_x(X);

			THEN("matrix is correct")
			{
				// clang-format off
				check_equal(N_to_x, "N_to_x", {
					{1, 1, 1, 1},
					{0, 0, 0, 1},
					{0, 1, 0, 0},
					{0, 0, 1, 0}
				});
				// clang-format on
			}

			AND_WHEN("derivative of natural wrt cartesian coords is calculated")
			{
				auto const & dN_by_dX = Tetrahedron::dN_by_dX(N_to_x);

				THEN("derivative is correct")
				{
					// clang-format off
					check_equal(dN_by_dX, "dN_by_dX", {
						{-1, -1, -1},
						{0, 1, 0},
						{-0, 0, 1},
						{1, -0, 0}
					 });
					// clang-format on
				}
			}

			AND_WHEN("derivative of cartesian wrt natural coords is calculated")
			{
				auto const & dx_by_dN = Tetrahedron::dx_by_dN(N_to_x);

				THEN("derivative is correct")
				{
					// clang-format off
					check_equal(dx_by_dN, "dx_by_dN", {
						{0, 0, 0, 1},
						{0, 1, 0, 0},
						{0, 0, 1, 0}

					});
					// clang-format on
				}
			}
		}
	}

	GIVEN("First element of a two-element mesh")
	{
		auto [X, x] = load_tet(file_name_two);
		(void)x;

		THEN("expected positions are reported")
		{
			Tetrahedron::Node::Positions expected;
			// clang-format off
			expected.setValues({
			   {0, 0, 0},
			   {0, 0, 1},
			   {1, 0, 0},
			   {0, 0.5, 0.5}
			});
			// clang-format on

			INFO("X:")
			INFO(X)
			CHECK(equal(X, expected));
		}

		WHEN("derivative of material wrt local coords is calculated")
		{
			auto const & dX_by_dL = Tetrahedron::dX_by_dL(X);

			THEN("derivative is correct")
			{
				// clang-format off
				check_equal(dX_by_dL, "dX_by_dL", {
					{0, 1, 0},
					{0, 0, 0.5},
					{1, 0, 0.5}
				});
				// clang-format on
			}

			AND_WHEN("derivative of local wrt material coords is calculated")
			{
				auto const & dL_by_dX = Tetrahedron::dL_by_dX(dX_by_dL);

				THEN("derivative is correct")
				{
					// clang-format off
					check_equal(dL_by_dX, "dL_by_dX", {
						{0, -1, 1},
						{1, 0, 0},
						{0, 2, 0}
					});
					// clang-format on
				}

				THEN("derivative is inverse of material wrt local")
				{
					constexpr Tetrahedron::IndexPairs<1> XL_LX{{{1, 0}}};
					Tetrahedron::MatrixTensor<3> const identity =
						dX_by_dL.contract(dL_by_dX, XL_LX);
					// clang-format off
					check_equal(identity, "dX_by_dL * dL_by_dX", {
						{1, 0, 0},
						{0, 1, 0},
						{0, 0, 1}
					});
					// clang-format on
				}
			}
		}

		AND_WHEN("derivative of natural wrt material coords is calculated")
		{
			auto const & dN_by_dX = Tetrahedron::dN_by_dX(X);

			THEN("derivative is correct")
			{
				// clang-format off
				check_equal(dN_by_dX, "dN_by_dX", {
					{-1, -1, -1},
					{0, -1, 1},
					{1, 0, 0},
					{0, 2, 0}
				});
				// clang-format on
			}
		}

		AND_WHEN("transformation from natural to cartesian coordinates is calculated")
		{
			auto const & N_to_x = Tetrahedron::N_to_x(X);

			THEN("matrix is correct")
			{
				// clang-format off
				check_equal(N_to_x, "N_to_x", {
					{1, 1, 1, 1},
					{0, 0, 1, 0},
					{0, 0, 0, 0.5},
					{0, 1, 0, 0.5}

				});
				// clang-format on
			}

			AND_WHEN("derivative of natural wrt material coords is calculated")
			{
				auto const & dN_by_dX = Tetrahedron::dN_by_dX(N_to_x);

				THEN("derivative is correct")
				{
					// clang-format off
					check_equal(dN_by_dX, "dN_by_dX", {
						{-1, -1, -1},
						{0, -1, 1},
						{1, 0, -0},
						{0, 2, 0}
					});
					// clang-format on
				}
			}
		}
	}
}


SCENARIO("Coordinate derivatives in deformed mesh")
{
	GIVEN("a one-element mesh")
	{
		auto [X, x] = load_tet(file_name_one);

		INFO("Material vertices:")
		INFO(X)

		WHEN("a node is deformed")
		{
			x(0, 0) += 0.5;
			INFO("Spatial vertices:")
			INFO(x)

			AND_WHEN("derivative of spatial wrt local coords is calculated")
			{
				auto const & dx_by_dL = Tetrahedron::dX_by_dL(x);

				THEN("derivative is correct")
				{
					// clang-format off
					check_equal(dx_by_dL, "dx_by_dL", {
						{-0.5, -0.5, 0.5},
						{1, 0, 0},
						{0, 1, 0}
					});
					// clang-format on
				}

				AND_WHEN("derivative of local wrt spatial coords is calculated")
				{
					auto const & dL_by_dx = Tetrahedron::dL_by_dX(dx_by_dL);

					THEN("derivative is correct")
					{
						// clang-format off
						check_equal(dL_by_dx, "dL_by_dx", {
							{0, 1, -0},
							{0, -0, 1},
							{2, 1, 1}
						});
						// clang-format on
					}
				}
			}

			AND_WHEN("derivative of shape function wrt spatial coords is calculated")
			{
				auto const & dN_by_dx = Tetrahedron::dN_by_dX(x);

				THEN("derivative is correct")
				{
					// clang-format off
					check_equal(dN_by_dx, "dN_by_dx", {
						{-2, -2, -2},
						{0, 1, 0},
						{0, 0, 1},
						{2, 1, 1}
					});
					// clang-format on
				}
			}
		}
	}
}


SCENARIO("Deformation gradient of undeformed element")
{
	GIVEN("first element of two-element mesh")
	{
		auto [X, x] = load_tet(file_name_two);

		INFO("Material vertices:")
		INFO(X)

		WHEN("deformation gradient is calculated from natural coordinate derivative")
		{
			auto const & dN_by_dX = Tetrahedron::dN_by_dX(X);
			auto const & dx_by_dX = Tetrahedron::dx_by_dX(X, dN_by_dX);

			THEN("gradient is identity")
			{
				// clang-format off
				check_equal(dx_by_dX, "dx_by_dX", {
					{1, 0, 0},
					{0, 1, 0},
					{0, 0, 1}
				});
				// clang-format on
			}
		}

		WHEN("deformation gradient is calculated from local coordinate derivative")
		{
			auto const & dX_by_dL = Tetrahedron::dX_by_dL(X);
			auto const & dL_by_dX = Tetrahedron::dL_by_dX(dX_by_dL);
			auto const & dx_by_dX = Tetrahedron::dx_by_dX(dX_by_dL, dL_by_dX);

			THEN("gradient is identity")
			{
				// clang-format off
				check_equal(dx_by_dX, "dx_by_dX", {
					{1, 0, 0},
					{0, 1, 0},
					{0, 0, 1}
				});
				// clang-format on
			}
		}

		WHEN("element is translated")
		{
			x(0, 0) += 0.5;
			x(1, 0) += 0.5;
			x(2, 0) += 0.5;
			x(3, 0) += 0.5;

			INFO("Spatial vertices:")
			INFO(x)

			AND_WHEN("deformation gradient is calculated")
			{
				 auto const & dN_by_dX = Tetrahedron::dN_by_dX(X);
				 auto const & dx_by_dX = Tetrahedron::dx_by_dX(x, dN_by_dX);

				 THEN("gradient is identity")
				 {
					 // clang-format off
					 check_equal(dx_by_dX, "dx_by_dX", {
						 {1, 0, 0},
						 {0, 1, 0},
						 {0, 0, 1}
					 });
					 // clang-format on
				 }
			}
		 }
	}
}

SCENARIO("Deformation gradient of deformed element")
{
	GIVEN("a one-element mesh")
	{
		auto [X, x] = load_tet(file_name_one);

		auto const & dN_by_dX = Tetrahedron::dN_by_dX(X);
		INFO("Material vertices:")
		INFO(X)

		WHEN("a node is deformed")
		{
			x(0, 0) += 0.5;
			INFO("Spatial vertices:")
			INFO(x)

			AND_WHEN("derivative of shape function wrt spatial coords is calculated")
			{
				auto const & dN_by_dx = Tetrahedron::dN_by_dX(x);

				THEN("derivative is correct")
				{
					// clang-format off
					check_equal(dN_by_dx, "dN_by_dx", {
						{-2, -2, -2},
						{0, 1, 0},
						{0, 0, 1},
						{2, 1, 1}
					});
					// clang-format on
				}
			}
			AND_WHEN("deformation gradient is calculated")
			{
				auto const & F = Tetrahedron::dx_by_dX(x, dN_by_dX);

				THEN("deformation gradient is correct")
				{
					// clang-format off
					check_equal(F, "F", {
						{0.5, -0.5, -0.5},
						{0, 1, 0},
						{0, 0, 1}
					});
					// clang-format on
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
					auto const & b = Tetrahedron::b(F);

					THEN("tensor is correct")
					{
						// clang-format off
						check_equal(b, "b", {
							{0.75, -0.5, -0.5},
							{-0.5,    1,    0},
							{-0.5,    0,    1}
						});
						// clang-format on
					}
				}
			}
		}
	}

	GIVEN("first element of two-element mesh")
	{
		auto [X, x] = load_tet(file_name_two);

		auto const & dN_by_dX = Tetrahedron::dN_by_dX(X);

		INFO("Material vertices:")
		INFO(X)

		WHEN("a node is deformed")
		{
			x(0, 0) += 0.5;

			AND_WHEN("derivative of shape function wrt spatial coords is calculated")
			{
				auto const & dN_by_dx = Tetrahedron::dN_by_dX(x);

				THEN("derivative is correct")
				{
					// clang-format off
					check_equal(dN_by_dx, "dN_by_dx", {
						{-2, -2, -2},
						{0, -1, 1},
						{2, 1, 1},
						{0, 2, 0}

					});
					// clang-format on
				}
			}
			AND_WHEN("deformation gradient is calculated")
			{
				auto const & F = Tetrahedron::dx_by_dX(x, dN_by_dX);

				THEN("deformation gradient is correct")
				{
					// clang-format off
					check_equal(F, "F", {
						{0.5, -0.5, -0.5},
						{0, 1, 0},
						{0, 0, 1}
					});
					// clang-format on
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
					auto const & b = Tetrahedron::b(F);

					THEN("tensor is correct")
					{
						// clang-format off
						check_equal(b, "b", {
							{0.75, -0.5, -0.5},
							{-0.5,    1,    0},
							{-0.5,    0,    1}
						});
						// clang-format on
					}
				}
			}
		}
	}
}

SCENARIO("Internal equivalent nodal force")
{
	GIVEN("an undeformed one-element mesh")
	{
		auto [X, x] = load_tet(file_name_one);

		auto const & dN_by_dX = Tetrahedron::dN_by_dX(X);
		//		Tetrahedron::Scalar const lambda = 1;
		//		Tetrahedron::Scalar const mu = 1;
		// Material properties: https://www.azom.com/properties.aspx?ArticleID=920
		double const mu = 0.4; // Shear modulus: 0.0003 - 0.02
		double const E = 1;	   // Young's modulus: 0.001 - 0.05
		// Lame's first parameter: https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
		double lambda = (mu * (E - 2 * mu)) / (3 * mu - E);

		INFO("Material vertices:")
		INFO(X)

		WHEN("neo-hookean stress is calculated")
		{
			auto const & F = Tetrahedron::dx_by_dX(x, dN_by_dX);
			auto const J = Tetrahedron::J(F);
			auto const & b = Tetrahedron::b(F);

			auto const & sigma = Tetrahedron::sigma(J, b, lambda, mu);

			THEN("stress is zero")
			{
				// clang-format off
				check_equal(sigma, "sigma", {
					{0, 0, 0},
					{0, 0, 0},
					{0, 0, 0}
				});
				// clang-format on
			}

			AND_WHEN("internal equivalent nodal forces are calculated")
			{
				auto const v = Tetrahedron::V(x);
				auto const & dN_by_dx = Tetrahedron::dN_by_dX(x);
				auto const & T = Tetrahedron::T(dN_by_dx, v, sigma);

				THEN("nodal forces are zero")
				{
					// clang-format off
					check_equal(T, "T", {
						{0, 0, 0},
						{0, 0, 0},
						{0, 0, 0}
					});
					// clang-format on
				}
			}
		}

		WHEN("a node is deformed")
		{
			x(0, 0) += 0.5;
			INFO("Spatial vertices:")
			INFO(x)

			AND_WHEN("neo-hookean stress is calculated")
			{
				auto const & F = Tetrahedron::dx_by_dX(x, dN_by_dX);
				auto const J = Tetrahedron::J(F);
				auto const & b = Tetrahedron::b(F);

				auto const & sigma = Tetrahedron::sigma(J, b, lambda, mu);

				THEN("stress is correct")
				{
					// clang-format off
					check_equal(sigma, "sigma", {
						{-0.754518, -0.4, -0.4},
						{-0.4, -0.554518, 0},
						{-0.4, 0, -0.554518}
					 });
					// clang-format on
				}

				AND_WHEN("internal equivalent nodal forces are calculated")
				{
					auto const v = Tetrahedron::V(x);
					auto const & dN_by_dx = Tetrahedron::dN_by_dX(x);
					auto const & T = Tetrahedron::T(dN_by_dx, v, sigma);

					THEN("nodal forces are correct")
					{
						// clang-format off
						check_equal(T, "T", {
							{0.259086, -0.0333333, -0.0333333, -0.19242},
							{0.159086, -0.0462098, 0, -0.112876},
							{0.159086, 0, -0.0462098, -0.112876}
						});
						// clang-format on
					}
				}
			}
		}
	}
}

SCENARIO("Neo-hookian tangent stiffness tensor")
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
		INFO(s.str())

		AND_GIVEN("an undeformed tetrahedron")
		{
			auto [X, x] = load_tet(file_name_one);

			auto const & dN_by_dX = Tetrahedron::dN_by_dX(X);
			auto const & F = Tetrahedron::dx_by_dX(x, dN_by_dX);
			auto const J = Tetrahedron::J(F);
			auto const & dN_by_dx = Tetrahedron::dN_by_dX(x);
			auto const v = Tetrahedron::V(x);

			INFO("Material vertices:")
			INFO(X)

			WHEN("neo-hookian elasticity tensor is calculated")
			{
				auto const & c = Tetrahedron::c(J, lambda, mu);

				THEN("it has expected values")
				{
					// clang-format off
					check_equal(c, "c",   {
						{
							{{1.2, 0, 0}, {0, 0.4, 0}, {0, 0, 0.4}},
							{{0, 0.4, 0}, {0.4, 0, 0}, {0, 0, 0}},
							{{0, 0, 0.4}, {0, 0, 0}, {0.4, 0, 0}}
						}, {
							{{0, 0.4, 0}, {0.4, 0, 0}, {0, 0, 0}},
							{{0.4, 0, 0}, {0, 1.2, 0}, {0, 0, 0.4}},
							{{0, 0, 0}, {0, 0, 0.4}, {0, 0.4, 0}}
						}, {
							{{0, 0, 0.4}, {0, 0, 0}, {0.4, 0, 0}},
							{{0, 0, 0}, {0, 0, 0.4}, {0, 0.4, 0}},
							{{0.4, 0, 0}, {0, 0.4, 0}, {0, 0, 1.2}}
						}
					});
					// clang-format on
				}
				AND_WHEN("constitutive component of tangent tensor is calculated")
				{
					auto const & Kc = Tetrahedron::Kc(dN_by_dx, v, c);

					THEN("it has expected values")
					{
						// clang-format off
						check_equal(Kc, "Kc", {
							{
								{{0.333333, -0.0666667, -0.0666667, -0.2}, {0.133333, -0.0666667, 0, -0.0666667}, {0.133333, 0, -0.0666667, -0.0666667}},
								{{-0.0666667, 0.0666667, 0, 0}, {-0.0666667, 0, 0, 0.0666667}, {0, 0, 0, 0}},
								{{-0.0666667, 0, 0.0666667, 0}, {0, 0, 0, 0}, {-0.0666667, 0, 0, 0.0666667}},
								{{-0.2, 0, 0, 0.2}, {-0.0666667, 0.0666667, 0, 0}, {-0.0666667, 0, 0.0666667, 0}}
							}, {
								{{0.133333, -0.0666667, 0, -0.0666667}, {0.333333, -0.2, -0.0666667, -0.0666667}, {0.133333, -0.0666667, -0.0666667, 0}},
								{{-0.0666667, 0, 0, 0.0666667}, {-0.2, 0.2, 0, 0}, {-0.0666667, 0, 0.0666667, 0}},
								{{0, 0, 0, 0}, {-0.0666667, 0, 0.0666667, 0}, {-0.0666667, 0.0666667, 0, 0}},
								{{-0.0666667, 0.0666667, 0, 0}, {-0.0666667, 0, 0, 0.0666667}, {0, 0, 0, 0}}
							}, {
								{{0.133333, 0, -0.0666667, -0.0666667}, {0.133333, -0.0666667, -0.0666667, 0}, {0.333333, -0.0666667, -0.2, -0.0666667}},
								{{0, 0, 0, 0}, {-0.0666667, 0, 0.0666667, 0}, {-0.0666667, 0.0666667, 0, 0}},
								{{-0.0666667, 0, 0, 0.0666667}, {-0.0666667, 0.0666667, 0, 0}, {-0.2, 0, 0.2, 0}},
								{{-0.0666667, 0, 0.0666667, 0}, {0, 0, 0, 0}, {-0.0666667, 0, 0, 0.0666667}}
							}

						});
						// clang-format on
					}
				}
			}

			WHEN("neo-hookian Cauchy stress tensor is calculated")
			{
				auto const & b = Tetrahedron::b(F);
				auto const & sigma = Tetrahedron::sigma(J, b, lambda, mu);

				AND_WHEN("initial stress component of tangent stiffness matrix is calculated")
				{
					Tetrahedron::StiffnessTensor const & Ks = Tetrahedron::Ks(dN_by_dx, v, sigma);

					THEN("stress component is zero")
					{
						Tetrahedron::StiffnessTensor expected;
						expected.setZero();
						check_equal(Ks, "Ks", expected, "zero");
					}
				}
			}
		} // End AND_GIVEN("an undeformed tetrahedron")

		AND_GIVEN("a deformed tetrahedron")
		{
			auto [X, x] = load_tet(file_name_one);
			x(0, 0) += 0.5;

			auto const & dN_by_dX = Tetrahedron::dN_by_dX(X);
			auto const & F = Tetrahedron::dx_by_dX(x, dN_by_dX);
			auto const J = Tetrahedron::J(F);
			auto const & dN_by_dx = Tetrahedron::dN_by_dX(x);
			auto const v = Tetrahedron::V(x);

			INFO("Material vertices:")
			INFO(X)
			INFO("Material volume:")
			INFO(Tetrahedron::V(X))
			INFO("Spatial vertices:")
			INFO(x)

			WHEN("neo-hookian elasticity tensor is calculated")
			{
				auto const & c = Tetrahedron::c(J, lambda, mu);

				THEN("it has expected values")
				{
					// clang-format off
					check_equal(c, "c",   {
						{
							{{3.50904, 0, 0}, {0, 0.8, 0}, {0, 0, 0.8}},
							{{0, 1.35452, 0}, {1.35452, 0, 0}, {0, 0, 0}},
							{{0, 0, 1.35452}, {0, 0, 0}, {1.35452, 0, 0}}
						}, {
							{{0, 1.35452, 0}, {1.35452, 0, 0}, {0, 0, 0}},
							{{0.8, 0, 0}, {0, 3.50904, 0}, {0, 0, 0.8}},
							{{0, 0, 0}, {0, 0, 1.35452}, {0, 1.35452, 0}}
						}, {
							{{0, 0, 1.35452}, {0, 0, 0}, {1.35452, 0, 0}},
							{{0, 0, 0}, {0, 0, 1.35452}, {0, 1.35452, 0}},
							{{0.8, 0, 0}, {0, 0.8, 0}, {0, 0, 3.50904}}
						}

					});
					// clang-format on
				}
				AND_WHEN("constitutive component of tangent tensor is calculated")
				{
					auto const & Kc = Tetrahedron::Kc(dN_by_dx, v, c);

					THEN("it has expected values")
					{
						// clang-format off
						check_equal(Kc, "Kc", {
							{
								{{2.07269, -0.225753, -0.225753, -1.62118}, {0.718173, -0.133333, 0, -0.584839}, {0.718173, 0, -0.133333, -0.584839}},
								{{-0.225753, 0.112876, 0, 0.112876}, {-0.225753, 0, 0, 0.225753}, {0, 0, 0, 0}},
								{{-0.225753, 0, 0.112876, 0.112876}, {0, 0, 0, 0}, {-0.225753, 0, 0, 0.225753}},
								{{-1.62118, 0.112876, 0.112876, 1.39543}, {-0.49242, 0.133333, 0, 0.359086}, {-0.49242, 0, 0.133333, 0.359086}}
							}, {
								{{0.718173, -0.225753, 0, -0.49242}, {2.07269, -0.584839, -0.225753, -1.2621}, {0.718173, -0.225753, -0.133333, -0.359086}},
								{{-0.133333, 0, 0, 0.133333}, {-0.584839, 0.29242, 0, 0.29242}, {-0.133333, 0, 0.0666667, 0.0666667}},
								{{0, 0, 0, 0}, {-0.225753, 0, 0.112876, 0.112876}, {-0.225753, 0.112876, 0, 0.112876}},
								{{-0.584839, 0.225753, 0, 0.359086}, {-1.2621, 0.29242, 0.112876, 0.856802}, {-0.359086, 0.112876, 0.0666667, 0.179543}}
							}, {
								{{0.718173, 0, -0.225753, -0.49242}, {0.718173, -0.133333, -0.225753, -0.359086}, {2.07269, -0.225753, -0.584839, -1.2621}},
								{{0, 0, 0, 0}, {-0.225753, 0, 0.112876, 0.112876}, {-0.225753, 0.112876, 0, 0.112876}},
								{{-0.133333, 0, 0, 0.133333}, {-0.133333, 0.0666667, 0, 0.0666667}, {-0.584839, 0, 0.29242, 0.29242}},
								{{-0.584839, 0, 0.225753, 0.359086}, {-0.359086, 0.0666667, 0.112876, 0.179543}, {-1.2621, 0.112876, 0.29242, 0.856802}}
							}

						});
						// clang-format on
					}
				}
			}

			WHEN("neo-hookian Cauchy stress tensor is calculated")
			{
				auto const & b = Tetrahedron::b(F);
				auto const & sigma = Tetrahedron::sigma(J, b, lambda, mu);

				THEN("stress is correct")
				{
					// clang-format off
					check_equal(sigma, "sigma", {
						{-0.754518, -0.4, -0.4},
						{-0.4, -0.554518, 0},
						{-0.4, 0, -0.554518}
					});
					// clang-format on
				}

				AND_WHEN("initial stress component of tangent stiffness matrix is calculated")
				{
					Tetrahedron::StiffnessTensor const & Ks = Tetrahedron::Ks(dN_by_dx, v, sigma);

					THEN("stress component is correct")
					{
						// clang-format off
						check_equal(Ks, "Ks", {
							{
								{{-1.15452, 0.159086, 0.159086, 0.836345}, {0, 0, 0, 0}, {0, 0, 0, 0}},
								{{0.159086, -0.0462098, 0, -0.112876}, {0, 0, 0, 0}, {0, 0, 0, 0}},
								{{0.159086, 0, -0.0462098, -0.112876}, {0, 0, 0, 0}, {0, 0, 0, 0}},
								{{0.836345, -0.112876, -0.112876, -0.610592}, {0, 0, 0, 0}, {0, 0, 0, 0}}
							}, {
								{{0, 0, 0, 0}, {-1.15452, 0.159086, 0.159086, 0.836345}, {0, 0, 0, 0}},
								{{0, 0, 0, 0}, {0.159086, -0.0462098, 0, -0.112876}, {0, 0, 0, 0}},
								{{0, 0, 0, 0}, {0.159086, 0, -0.0462098, -0.112876}, {0, 0, 0, 0}},
								{{0, 0, 0, 0}, {0.836345, -0.112876, -0.112876, -0.610592}, {0, 0, 0, 0}}
							}, {
								{{0, 0, 0, 0}, {0, 0, 0, 0}, {-1.15452, 0.159086, 0.159086, 0.836345}},
								{{0, 0, 0, 0}, {0, 0, 0, 0}, {0.159086, -0.0462098, 0, -0.112876}},
								{{0, 0, 0, 0}, {0, 0, 0, 0}, {0.159086, 0, -0.0462098, -0.112876}},
								{{0, 0, 0, 0}, {0, 0, 0, 0}, {0.836345, -0.112876, -0.112876, -0.610592}}
							}
						});
						// clang-format on

					}
				}
			}

		} // End AND_GIVEN("an undeformed tetrahedron")
	}
}

SCENARIO("Solution of a single element")
{
	GIVEN("tangent stiffness and equivalent node force tensors")
	{
		using namespace FeltElements;
		// Material properties: https://www.azom.com/properties.aspx?ArticleID=920
		double const mu = 0.4; // Shear modulus: 0.0003 - 0.02
		double const E = 1;	   // Young's modulus: 0.001 - 0.05
		// Lame's first parameter: https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
		double lambda = (mu * (E - 2 * mu)) / (3 * mu - E);

		auto [X, x] = load_tet(file_name_one);
		x(0, 0) += 0.5;
		auto const & dN_by_dX = Tetrahedron::dN_by_dX(X);

		std::stringstream s;
		s << "Lambda = " << lambda << "; mu = " << mu;
		INFO(s.str())
		INFO("Material vertices:")
		INFO(X)

		WHEN("displacement is solved")
		{
			std::stringstream ss;

			for (int i = 0; i < 7; i++)
			{
				auto [K, T] = Tetrahedron::KT(x, dN_by_dX, lambda, mu);

				using Displacements = Eigen::Matrix<double, 12, 1>;
				Eigen::Map<Displacements> B{T.data(), 12, 1};
				B = -B;

				Eigen::Map<
					Eigen::Matrix<double, 12, 12, Tetrahedron::StiffnessTensor::Layout>> const
					A{K.data(), 12, 12};

				Displacements u = A.colPivHouseholderQr().solve(B);

				FeltElements::Tetrahedron::Node::Positions delta =
					Eigen::TensorMap<Eigen::Tensor<double, 2>>{u.data(), {3, 4}}.shuffle(
						Eigen::array<Eigen::Index, 2>{1, 0});

				ss << ">>>>>>>>>>>>>>>>>>>>>>> Iteration " << i << "\n\n";
				ss << "K (tensor)" << "\n";
				ss << K << "\n";
				ss << "v" << "\n";
				ss << Tetrahedron::V(x) << "\n";
				ss << "-T" << "\n";
				ss << B << "\n";
				ss << "K" << "\n";
				ss << A << "\n";
				ss << "u" << "\n";
				ss << u << "\n";
				ss << "x" << "\n";
				ss << x << "\n";
				ss << "delta" << "\n";
				ss << delta << "\n";

				x += delta;
			}
			INFO(ss.str())

			THEN("volume returns and strain is zero")
			{
				CHECK(Tetrahedron::V(x) == Approx(1.0 / 6));
				// clang-format off
				check_equal(Tetrahedron::b(Tetrahedron::dx_by_dX(x, dN_by_dX)), "b", {
					{1, 0, 0},
					{0, 1, 0},
					{0, 0, 1}
				});
				check_equal(x, "x", {
					{0.287126, 0.0989334, 0.0478705},
					{-0.0143856, 1.0512, 0.000137554},
					{-0.0143856, 0.0512004, 1.00014},
					{1.19166, 0.400445, 0.349382}
				});
				// clang-format on
			}
		}
	}
}


SCENARIO("Mesh attributes")
{
	GIVEN("a two-element mesh")
	{
		auto const & io = TetGenIO{file_name_two};
		auto mesh = io.to_mesh();

		WHEN("nodal force attributes are constructed")
		{
			FeltElements::Attribute::Node::Force const attribs{mesh};

			THEN("properties are default initialised")
			{
				for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
					check_equal(
						attribs[*itvtxh], "force",
						Eigen::TensorFixedSize<double, Eigen::Sizes<3>>{}.setZero(), "zero");
			}
		}

		WHEN("spatial position attributes are constructed")
		{
			FeltElements::Attribute::Node::SpatialPosition const attribs{mesh};

			THEN("attributes are initialised to material position")
			{
				for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
				{
					check_equal(
						attribs[*itvtxh],
						"x",
						reinterpret_cast<Tetrahedron::Node::Pos const &>(mesh.vertex(*itvtxh)),
						"X");
				}
			}
		}

		WHEN("element vertex mapping attributes are constructed")
		{
			FeltElements::Attribute::Element::VertexHandles const vtxh_attribs{mesh};

			THEN("properties are populated")
			{
				auto itvtxh = mesh.cells_begin();
				CHECK(
					vtxh_attribs[*itvtxh] ==
					std::array<OpenVolumeMesh::VertexHandle, 4>{0, 3, 1, 4});
				itvtxh++;
				CHECK(
					vtxh_attribs[*itvtxh] ==
					std::array<OpenVolumeMesh::VertexHandle, 4>{2, 0, 1, 4});
			}

			AND_WHEN("element natural wrt material coords attributes are constructed")
			{
				FeltElements::Attribute::Element::MaterialShapeDerivative const dn_by_dX_attribs{
					mesh, vtxh_attribs};

				THEN("properties are populated")
				{
					auto itvtxh = mesh.cells_begin();
					check_equal(
						dn_by_dX_attribs[*itvtxh],
						std::string{"dN_by_dX "} + std::to_string(itvtxh),
						Tetrahedron::dN_by_dX(
							  Tetrahedron::X(mesh, Tetrahedron::vtxhs(mesh, *itvtxh))),
						"dN_by_dX");
					itvtxh++;
					check_equal(
						dn_by_dX_attribs[*itvtxh],
						std::string{"dN_by_dX "} + std::to_string(itvtxh),
						Tetrahedron::dN_by_dX(
							Tetrahedron::X(mesh, Tetrahedron::vtxhs(mesh, *itvtxh))),
						"dN_by_dX");
				}
			}
		}
	}
}
