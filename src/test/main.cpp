#define EIGEN_DEFAULT_IO_FORMAT Eigen::IOFormat(3, 0, ", ", ",\n", "", "", "", "")
#include <FeltElements/Solver.hpp>
#include <boost/range/irange.hpp>

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_CONSOLE_WIDTH 200
#include "common.hpp"

char const * const file_name_one = "resources/one.ovm";
char const * const file_name_two = "resources/two.ovm";

using namespace FeltElements;

//	char cwd[500];
//	getcwd(cwd, 500);
//	std::cerr << "Executing tests in " << cwd << std::endl;

SCENARIO("Mesh attributes")
{
	GIVEN("a two-element mesh")
	{
		auto mesh = load_ovm_mesh(file_name_two);

		WHEN("element nodal force attributes are constructed")
		{
			Element::Attribute::NodalForces const attrib_forces{mesh};

			THEN("properties are default initialised")
			{
				auto itcellh = mesh.cells_begin();
				Node::Forces T = attrib_forces[*itcellh];
				CHECK(Tensor::Func::all_of(T == 0));
				itcellh++;
				T = attrib_forces[*itcellh];
				CHECK(Tensor::Func::all_of(T == 0));
			}
		}

		WHEN("element stiffness attributes are constructed")
		{
			Element::Attribute::Stiffness const attrib_stiffness{mesh};

			THEN("properties are default initialised")
			{
				auto itcellh = mesh.cells_begin();
				Element::Stiffness K = attrib_stiffness[*itcellh];
				CHECK(Tensor::Func::all_of(K == 0));
				itcellh++;
				K = attrib_stiffness[*itcellh];
				CHECK(Tensor::Func::all_of(K == 0));
			}
		}

		WHEN("element vertex mapping attributes are constructed")
		{
			FeltElements::Element::Attribute::VertexHandles const attrib_vtxhs{mesh};

			THEN("properties are populated")
			{
				auto itcellh = mesh.cells_begin();
				CHECK(
					attrib_vtxhs[*itcellh] ==
					std::array<OpenVolumeMesh::VertexHandle, 4>{0, 3, 1, 4});
				itcellh++;
				CHECK(
					attrib_vtxhs[*itcellh] ==
					std::array<OpenVolumeMesh::VertexHandle, 4>{2, 0, 1, 4});
			}

			AND_WHEN("element natural wrt material coords attributes are constructed")
			{
				Node::Attribute::MaterialPosition const attrib_X{mesh};
				FeltElements::Element::Attribute::MaterialShapeDerivative const attrib_dN_by_dX{
					mesh, attrib_vtxhs, attrib_X};

				THEN("properties are populated")
				{
					auto itcellh = mesh.cells_begin();
					check_equal(
						attrib_dN_by_dX[*itcellh],
						std::string{"dN_by_dX "} + std::to_string(itcellh),
						Derivatives::dN_by_dX(attrib_X.for_element(attrib_vtxhs[*itcellh])));
					itcellh++;
					check_equal(
						attrib_dN_by_dX[*itcellh],
						std::string{"dN_by_dX "} + std::to_string(itcellh),
						Derivatives::dN_by_dX(attrib_X.for_element(attrib_vtxhs[*itcellh])));
				}
			}
		}

		WHEN("spatial position attributes are constructed")
		{
			Node::Attribute::SpatialPosition const attrib_x{mesh};

			THEN("attributes are initialised to material position")
			{
				for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
				{
					check_equal(
						attrib_x[*itvtxh],
						"x",
						reinterpret_cast<Node::Pos const &>(mesh.vertex(*itvtxh)),
						"X");
				}
			}

			AND_WHEN("element nodal spatial position tensor is constructed")
			{
				Element::Attribute::VertexHandles const attrib_vtxhs{mesh};

				auto itcellh = mesh.cells_begin();
				Node::Positions x1 = attrib_x.for_element(attrib_vtxhs[*itcellh]);
				itcellh++;
				Node::Positions x2 = attrib_x.for_element(attrib_vtxhs[*itcellh]);

				THEN("positions are expected")
				{
					// clang-format off
					check_equal(x1, "x1", {
						{0.000000, 0.000000, 0.000000},
						{0.000000, 0.000000, 1.000000},
						{1.000000, 0.000000, 0.000000},
						{0.000000, 0.500000, 0.500000}
					});
					check_equal(x2, "x2", {
						{0.000000, 1.000000, 0.000000},
						{0.000000, 0.000000, 0.000000},
						{1.000000, 0.000000, 0.000000},
						{0.000000, 0.500000, 0.500000}
					});
					// clang-format on
				}
			}
		}
	}
}

SCENARIO("Metrics of undeformed mesh")
{
	GIVEN("one-element mesh")
	{
		FeltElements::Mesh mesh;
		OpenVolumeMesh::IO::FileManager{}.readFile(file_name_one, mesh);

		WHEN("vertex index mapping is fetched")
		{
			Element::Attribute::VertexHandles const attrib_vtxhs{mesh};
			auto const & vtxhs = attrib_vtxhs[0];

			THEN("mapping is expected")
			{
				CHECK(vtxhs == Element::Vtxhs{ 0, 1, 2, 3 });
			}

			AND_WHEN("material node position tensor is constructed")
			{
				Node::Attribute::MaterialPosition const attrib_X{mesh};
				Node::Positions const X = attrib_X.for_element(vtxhs);

				THEN("expected positions are reported")
				{
					// clang-format off
					check_equal(X, "X", {
						{0, 0, 0},
						{0, 0, 1},
						{1, 0, 0},
						{0, 1, 0}
					});
					// clang-format on
				}

				AND_WHEN("material volume is calculated")
				{
					Scalar V = Derivatives::V(X);

					THEN("volume is correct")
					{
						CHECK(V == 1.0 / 6.0);
					}
				}

				AND_WHEN("spatial coordinate property is created")
				{
					using namespace OpenVolumeMesh;
					Node::Attribute::SpatialPosition const attrib_x{mesh};

					THEN("spatial coordinates equal material coordinates")
					{
						std::vector<Vec3d> all_X{};
						std::vector<Vec3d> all_x{};

						for (VertexIter it_vtx = mesh.vertices_begin();
							 it_vtx != mesh.vertices_end();
							 it_vtx++)
						{
							auto const & x = attrib_x[*it_vtx];
							all_X.push_back(mesh.vertex(*it_vtx));
							all_x.emplace_back(x(0), x(1), x(2));
						}

						CHECK(all_X == all_x);
					}

					AND_WHEN("spatial node position tensor is constructed")
					{
						Node::Attribute::SpatialPosition attrib_x{mesh};
						auto const x = attrib_x.for_element(vtxhs);

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
					{0.5, 0.0, 0.0},
					{0.0, 0.0, 1.0},
					{1.0, 0.0, 0.0},
					{0.0, 1.0, 0.0}
				});
				// clang-format on
			}

			AND_WHEN("spatial volume is calculated")
			{
				Scalar v = Derivatives::V(x);

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
			Element::Attribute::MaterialShapeDerivative::dN_by_dL,
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
			Element::Attribute::MaterialShapeDerivative::dL_by_dN,
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
		enum {i, j, k};
		using namespace Tensor;
		using namespace Tensor::Func;
		Matrix<3> const identity = einsum<Indices<i, k>, Indices<k, j>>(
			Element::Attribute::MaterialShapeDerivative::dL_by_dN,
			Element::Attribute::MaterialShapeDerivative::dN_by_dL);

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
			auto const & dX_by_dL = Derivatives::dX_by_dL(X);

			THEN("derivative is correct")
			{
				// clang-format off
				check_equal(dX_by_dL, "dX_by_dL", {
					{0, 1, 0},
					{0, 0, 1},
					{1, 0, 0}
				});
				// clang-format on
			}

			AND_WHEN("derivative of local wrt material coords is calculated")
			{
				auto const & dL_by_dX = Derivatives::dL_by_dX(dX_by_dL);

				THEN("derivative is correct")
				{
					// clang-format off
					check_equal(dL_by_dX, "dL_by_dX", {
						{0, 0, 1},
						{1, 0, 0},
						{0, 1, 0}
					});
					// clang-format on
				}

				AND_WHEN("derivative of material wrt shape coords is calculated")
				{
					auto const & dN_by_dX = Derivatives::dN_by_dX(X);

					THEN("derivative is correct")
					{
						// clang-format off
						check_equal(dN_by_dX, "dN_by_dX", {
							{-1, -1, -1},
							{0, 0, 1},
							{1, 0, 0},
							{0, 1, 0}
						});
						// clang-format on
					}
				}
			}
		}

		AND_WHEN("transformation from natural to cartesian coordinates is calculated")
		{
			auto const & N_to_x = Derivatives::N_to_x(X);

			THEN("matrix is correct")
			{
				// clang-format off
				check_equal(N_to_x, "N_to_x", {
					{1, 1, 1, 1},
					{0, 0, 1, 0},
					{0, 0, 0, 1},
					{0, 1, 0, 0}
				});
				// clang-format on
			}

			AND_WHEN("derivative of natural wrt cartesian coords is calculated")
			{
				auto const & dN_by_dX = Derivatives::dN_by_dX(N_to_x);

				THEN("derivative is correct")
				{
					// clang-format off
					check_equal(dN_by_dX, "dN_by_dX", {
						{-1, -1, -1},
						{0, 0, 1},
						{1, 0, 0},
						{0, 1, 0}
					 });
					// clang-format on
				}
			}

			AND_WHEN("derivative of cartesian wrt natural coords is calculated")
			{
				auto const & dx_by_dN = Derivatives::dx_by_dN(N_to_x);

				THEN("derivative is correct")
				{
					// clang-format off
					check_equal(dx_by_dN, "dx_by_dN", {
						{0, 0, 1, 0},
						{0, 0, 0, 1},
						{0, 1, 0, 0}
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
			Node::Positions expected{
			   {0.0, 0.0, 0.0},
			   {0.0, 0.0, 1.0},
			   {1.0, 0.0, 0.0},
			   {0.0, 0.5, 0.5}
			};
			// clang-format on

			INFO("X:")
			INFO(X)
			CHECK(equal(X, expected));
		}

		WHEN("derivative of material wrt local coords is calculated")
		{
			auto const & dX_by_dL = Derivatives::dX_by_dL(X);

			THEN("derivative is correct")
			{
				// clang-format off
				check_equal(dX_by_dL, "dX_by_dL", {
					{0.0, 1.0, 0.0},
					{0.0, 0.0, 0.5},
					{1.0, 0.0, 0.5}
				});
				// clang-format on
			}

			AND_WHEN("derivative of local wrt material coords is calculated")
			{
				auto const & dL_by_dX = Derivatives::dL_by_dX(dX_by_dL);

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
					using namespace Tensor;
					using namespace Tensor::Func;
					Matrix<3> const identity = einsum<Idxs<i, k>, Idxs<k, j>>(
						dX_by_dL, dL_by_dX);
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
			auto const & dN_by_dX = Derivatives::dN_by_dX(X);

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
			auto const & N_to_x = Derivatives::N_to_x(X);

			THEN("matrix is correct")
			{
				// clang-format off
				check_equal(N_to_x, "N_to_x", {
					{1.0, 1.0, 1.0, 1.0},
					{0.0, 0.0, 1.0, 0.0},
					{0.0, 0.0, 0.0, 0.5},
					{0.0, 1.0, 0.0, 0.5}
				});
				// clang-format on
			}

			AND_WHEN("derivative of natural wrt material coords is calculated")
			{
				auto const & dN_by_dX = Derivatives::dN_by_dX(N_to_x);

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
				auto const & dx_by_dL = Derivatives::dX_by_dL(x);

				THEN("derivative is correct")
				{
					// clang-format off
					check_equal(dx_by_dL, "dx_by_dL", {
						{-0.5, 0.5, -0.5},
						{0.0, 0.0, 1.0},
						{1.0, 0.0, 0.0}
					});
					// clang-format on
				}

				AND_WHEN("derivative of local wrt spatial coords is calculated")
				{
					auto const & dL_by_dx = Derivatives::dL_by_dX(dx_by_dL);

					THEN("derivative is correct")
					{
						// clang-format off
						check_equal(dL_by_dx, "dL_by_dx", {
							{0, 0, 1},
							{2, 1, 1},
							{0, 1, 0}
						});
						// clang-format on
					}
				}
			}

			AND_WHEN("derivative of shape function wrt spatial coords is calculated")
			{
				auto const & dN_by_dx = Derivatives::dN_by_dX(x);

				THEN("derivative is correct")
				{
					// clang-format off
					check_equal(dN_by_dx, "dN_by_dx", {
						{-2, -2, -2},
						{0, 0, 1},
						{2, 1, 1},
						{0, 1, 0}
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
			auto const & dN_by_dX = Derivatives::dN_by_dX(X);
			auto const & dx_by_dX = Derivatives::dx_by_dX(X, dN_by_dX);

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
			auto const & dX_by_dL = Derivatives::dX_by_dL(X);
			auto const & dL_by_dX = Derivatives::dL_by_dX(dX_by_dL);
			auto const & dx_by_dX = Derivatives::dx_by_dX(dX_by_dL, dL_by_dX);

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
				 auto const & dN_by_dX = Derivatives::dN_by_dX(X);
				 auto const & dx_by_dX = Derivatives::dx_by_dX(x, dN_by_dX);

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

		auto const & dN_by_dX = Derivatives::dN_by_dX(X);
		INFO("Material vertices:")
		INFO(X)

		WHEN("a node is deformed")
		{
			x(0, 0) += 0.5;
			INFO("Spatial vertices:")
			INFO(x)

			AND_WHEN("derivative of shape function wrt spatial coords is calculated")
			{
				auto const & dN_by_dx = Derivatives::dN_by_dX(x);

				THEN("derivative is correct")
				{
					// clang-format off
					check_equal(dN_by_dx, "dN_by_dx", {
						{-2, -2, -2},
						{0, 0, 1},
						{2, 1, 1},
						{0, 1, 0}
					});
					// clang-format on
				}
			}
			AND_WHEN("deformation gradient is calculated")
			{
				auto const & F = Derivatives::dx_by_dX(x, dN_by_dX);

				THEN("deformation gradient is correct")
				{
					// clang-format off
					check_equal(F, "F", {
						{0.5, -0.5, -0.5},
						{0.0, 1.0, 0.0},
						{0.0, 0.0, 1.0}
					});
					// clang-format on
				}

				AND_WHEN("Jacobian of deformation gradient is calculated")
				{
					double const J = Derivatives::J(F);

					THEN("value is correct")
					{
						CHECK(J == 0.5);
					}
				}

				AND_WHEN("Left Cauchy-Green / Finger tensor is calculated")
				{
					auto const & b = Derivatives::b(F);

					THEN("tensor is correct")
					{
						// clang-format off
						check_equal(b, "b", {
							{0.75, -0.5, -0.5},
							{-0.5,    1.0,    0.0},
							{-0.5,    0.0,    1.0}
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

		auto const & dN_by_dX = Derivatives::dN_by_dX(X);

		INFO("Material vertices:")
		INFO(X)

		WHEN("a node is deformed")
		{
			x(0, 0) += 0.5;

			AND_WHEN("derivative of shape function wrt spatial coords is calculated")
			{
				auto const & dN_by_dx = Derivatives::dN_by_dX(x);

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
				auto const & F = Derivatives::dx_by_dX(x, dN_by_dX);

				THEN("deformation gradient is correct")
				{
					// clang-format off
					check_equal(F, "F", {
						{0.5, -0.5, -0.5},
						{0.0, 1.0, 0.0},
						{0.0, 0.0, 1.0}
					});
					// clang-format on
				}

				AND_WHEN("Jacobian of deformation gradient is calculated")
				{
					double const J = Derivatives::J(F);

					THEN("value is correct")
					{
						CHECK(J == 0.5);
					}
				}

				AND_WHEN("Left Cauchy-Green / Finger tensor is calculated")
				{
					auto const & b = Derivatives::b(F);

					THEN("tensor is correct")
					{
						// clang-format off
						check_equal(b, "b", {
							{0.75, -0.5, -0.5},
							{-0.5,    1.0,    0.0},
							{-0.5,    0.0,    1.0}
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

		INFO("Material vertices:")
		INFO(X)

		auto const & dN_by_dX = Derivatives::dN_by_dX(X);
		//		Tetrahedron::Scalar const lambda = 1;
		//		Tetrahedron::Scalar const mu = 1;
		// Material properties: https://www.azom.com/properties.aspx?ArticleID=920
		double const mu = 0.4; // Shear modulus: 0.0003 - 0.02
		double const E = 1;	   // Young's modulus: 0.001 - 0.05
		// Lame's first parameter: https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
		double lambda = (mu * (E - 2 * mu)) / (3 * mu - E);

		WHEN("neo-hookean stress is calculated")
		{
			auto const & F = Derivatives::dx_by_dX(x, dN_by_dX);
			auto const J = Derivatives::J(F);
			auto const & b = Derivatives::b(F);

			auto const & sigma = Derivatives::sigma(J, b, lambda, mu);

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
				auto const v = Derivatives::V(x);
				auto const & dN_by_dx = Derivatives::dN_by_dX(x);
				auto const & T = Derivatives::T(dN_by_dx, v, sigma);

				THEN("nodal forces are zero")
				{
					// clang-format off
					check_equal(T, "T", {
						{0, 0, 0},
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
				auto const & F = Derivatives::dx_by_dX(x, dN_by_dX);
				auto const J = Derivatives::J(F);
				auto const & b = Derivatives::b(F);

				auto const & sigma = Derivatives::sigma(J, b, lambda, mu);

				THEN("stress is correct")
				{
					// clang-format off
					check_equal(sigma, "sigma", {
						{-0.754518, -0.4, -0.4},
						{-0.4, -0.554518, 0.0},
						{-0.4, 0.0, -0.554518}
					 });
					// clang-format on
				}

				AND_WHEN("internal equivalent nodal forces are calculated")
				{
					auto const v = Derivatives::V(x);
					auto const & dN_by_dx = Derivatives::dN_by_dX(x);
					auto const & T = Derivatives::T(dN_by_dx, v, sigma);

					INFO("Derivative of natural wrt spatial position")
					INFO(dN_by_dx)

					THEN("nodal forces are correct")
					{
//						Node::Forces const & zero = 0;
//						check_equal(T, "T", zero, "zero");
						// clang-format off
						check_equal(T, "T", {
							{0.259086, 0.159086, 0.159086},
							{-0.0333333, 0.0, -0.0462098},
							{-0.19242, -0.112876, -0.112876},
							{-0.0333333, -0.0462098, 0.0}
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

			auto const & dN_by_dX = Derivatives::dN_by_dX(X);
			auto const & F = Derivatives::dx_by_dX(x, dN_by_dX);
			auto const J = Derivatives::J(F);
			auto const & dN_by_dx = Derivatives::dN_by_dX(x);
			auto const v = Derivatives::V(x);

			INFO("Material vertices:")
			INFO(X)

			WHEN("neo-hookian elasticity tensor is calculated")
			{
				auto const & c = Derivatives::c(J, lambda, mu);

				THEN("it has expected values")
				{
					// clang-format off
					check_equal(c, "c",   {
						{
							{{1.2, 0.0, 0.0}, {0.0, 0.4, 0.0}, {0.0, 0.0, 0.4}},
							{{0.0, 0.4, 0.0}, {0.4, 0.0, 0.0}, {0.0, 0.0, 0.0}},
							{{0.0, 0.0, 0.4}, {0.0, 0.0, 0.0}, {0.4, 0.0, 0.0}}
						}, {
							{{0.0, 0.4, 0.0}, {0.4, 0.0, 0.0}, {0.0, 0.0, 0.0}},
							{{0.4, 0.0, 0.0}, {0.0, 1.2, 0.0}, {0.0, 0.0, 0.4}},
							{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.4}, {0.0, 0.4, 0.0}}
						}, {
							{{0.0, 0.0, 0.4}, {0.0, 0.0, 0.0}, {0.4, 0.0, 0.0}},
							{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.4}, {0.0, 0.4, 0.0}},
							{{0.4, 0.0, 0.0}, {0.0, 0.4, 0.0}, {0.0, 0.0, 1.2}}
						}
					});
					// clang-format on
				}
				AND_WHEN("constitutive component of tangent tensor is calculated")
				{
					auto const & Kc = Derivatives::Kc(dN_by_dx, v, c);

					THEN("it has expected values")
					{
//						Element::Stiffness const & zero = 0;
//						check_equal(Kc, "Kc", zero, "zero");
						// clang-format off
						check_equal(Kc, "Kc", {
							{
								{{0.333333, 0.133333, 0.133333}, {-0.066667, 0.000000, -0.066667}, {-0.200000, -0.066667, -0.066667}, {-0.066667, -0.066667, 0.000000}},
								{{0.133333, 0.333333, 0.133333}, {0.000000, -0.066667, -0.066667}, {-0.066667, -0.066667, 0.000000}, {-0.066667, -0.200000, -0.066667}},
								{{0.133333, 0.133333, 0.333333}, {-0.066667, -0.066667, -0.200000}, {-0.066667, 0.000000, -0.066667}, {0.000000, -0.066667, -0.066667}}
							}, {
								{{-0.066667, 0.000000, -0.066667}, {0.066667, 0.000000, 0.000000}, {0.000000, 0.000000, 0.066667}, {0.000000, 0.000000, 0.000000}},
								{{0.000000, -0.066667, -0.066667}, {0.000000, 0.066667, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.066667}},
								{{-0.066667, -0.066667, -0.200000}, {0.000000, 0.000000, 0.200000}, {0.066667, 0.000000, 0.000000}, {0.000000, 0.066667, 0.000000}}
							}, {
								{{-0.200000, -0.066667, -0.066667}, {0.000000, 0.000000, 0.066667}, {0.200000, 0.000000, 0.000000}, {0.000000, 0.066667, 0.000000}},
								{{-0.066667, -0.066667, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.066667, 0.000000}, {0.066667, 0.000000, 0.000000}},
								{{-0.066667, 0.000000, -0.066667}, {0.066667, 0.000000, 0.000000}, {0.000000, 0.000000, 0.066667}, {0.000000, 0.000000, 0.000000}}
							}, {
								{{-0.066667, -0.066667, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.066667, 0.000000}, {0.066667, 0.000000, 0.000000}},
								{{-0.066667, -0.200000, -0.066667}, {0.000000, 0.000000, 0.066667}, {0.066667, 0.000000, 0.000000}, {0.000000, 0.200000, 0.000000}},
								{{0.000000, -0.066667, -0.066667}, {0.000000, 0.066667, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.066667}}
							}
						});
						// clang-format on
					}
				}
			}

			WHEN("neo-hookian Cauchy stress tensor is calculated")
			{
				auto const & b = Derivatives::b(F);
				auto const & sigma = Derivatives::sigma(J, b, lambda, mu);

				AND_WHEN("initial stress component of tangent stiffness matrix is calculated")
				{
					Element::Stiffness const & Ks = Derivatives::Ks(dN_by_dx, v, sigma);

					THEN("stress component is zero")
					{
						Element::Stiffness const expected = 0;
						check_equal(Ks, "Ks", expected, "zero");
					}
				}
			}
		} // End AND_GIVEN("an undeformed tetrahedron")

		AND_GIVEN("a deformed tetrahedron")
		{
			auto [X, x] = load_tet(file_name_one);
			x(0, 0) += 0.5;

			auto const & dN_by_dX = Derivatives::dN_by_dX(X);
			auto const & F = Derivatives::dx_by_dX(x, dN_by_dX);
			auto const J = Derivatives::J(F);
			auto const & dN_by_dx = Derivatives::dN_by_dX(x);
			auto const v = Derivatives::V(x);

			INFO("Material vertices:")
			INFO(X)
			INFO("Material volume:")
			INFO(Derivatives::V(X))
			INFO("Spatial vertices:")
			INFO(x)

			WHEN("neo-hookian elasticity tensor is calculated")
			{
				auto const & c = Derivatives::c(J, lambda, mu);

				THEN("it has expected values")
				{
					// clang-format off
					check_equal(c, "c",   {
						{
							{{3.50904, 0.0, 0.0}, {0.0, 0.8, 0.0}, {0.0, 0.0, 0.8}},
							{{0.0, 1.35452, 0.0}, {1.35452, 0.0, 0.0}, {0.0, 0.0, 0.0}},
							{{0.0, 0.0, 1.35452}, {0.0, 0.0, 0.0}, {1.35452, 0.0, 0.0}}
						}, {
							{{0.0, 1.35452, 0.0}, {1.35452, 0.0, 0.0}, {0.0, 0.0, 0.0}},
							{{0.8, 0.0, 0.0}, {0.0, 3.50904, 0.0}, {0.0, 0.0, 0.8}},
							{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.35452}, {0.0, 1.35452, 0.0}}
						}, {
							{{0.0, 0.0, 1.35452}, {0.0, 0.0, 0.0}, {1.35452, 0.0, 0.0}},
							{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.35452}, {0.0, 1.35452, 0.0}},
							{{0.8, 0.0, 0.0}, {0.0, 0.8, 0.0}, {0.0, 0.0, 3.50904}}
						}

					});
					// clang-format on
				}
				AND_WHEN("constitutive component of tangent tensor is calculated")
				{
					auto const & Kc = Derivatives::Kc(dN_by_dx, v, c);

					THEN("it has expected values")
					{
//						Element::Stiffness const & check = 0;
//						check_equal(Kc, "Kc", check, "zero");
						// clang-format off
						check_equal(Kc, "Kc", {
							{
								{{2.072690, 0.718173, 0.718173}, {-0.225753, 0.000000, -0.133333}, {-1.621184, -0.584839, -0.584839}, {-0.225753, -0.133333, 0.000000}},
								{{0.718173, 2.072690, 0.718173}, {0.000000, -0.225753, -0.133333}, {-0.492420, -1.262098, -0.359086}, {-0.225753, -0.584839, -0.225753}},
								{{0.718173, 0.718173, 2.072690}, {-0.225753, -0.225753, -0.584839}, {-0.492420, -0.359086, -1.262098}, {0.000000, -0.133333, -0.225753}}
							}, {
								{{-0.225753, 0.000000, -0.225753}, {0.112876, 0.000000, 0.000000}, {0.112876, 0.000000, 0.225753}, {0.000000, 0.000000, 0.000000}},
								{{0.000000, -0.225753, -0.225753}, {0.000000, 0.112876, 0.000000}, {0.000000, 0.112876, 0.112876}, {0.000000, 0.000000, 0.112876}},
								{{-0.133333, -0.133333, -0.584839}, {0.000000, 0.000000, 0.292420}, {0.133333, 0.066667, 0.292420}, {0.000000, 0.066667, 0.000000}}
							}, {
								{{-1.621184, -0.492420, -0.492420}, {0.112876, 0.000000, 0.133333}, {1.395431, 0.359086, 0.359086}, {0.112876, 0.133333, 0.000000}},
								{{-0.584839, -1.262098, -0.359086}, {0.000000, 0.112876, 0.066667}, {0.359086, 0.856802, 0.179543}, {0.225753, 0.292420, 0.112876}},
								{{-0.584839, -0.359086, -1.262098}, {0.225753, 0.112876, 0.292420}, {0.359086, 0.179543, 0.856802}, {0.000000, 0.066667, 0.112876}}
							}, {
								{{-0.225753, -0.225753, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.112876, 0.225753, 0.000000}, {0.112876, 0.000000, 0.000000}},
								{{-0.133333, -0.584839, -0.133333}, {0.000000, 0.000000, 0.066667}, {0.133333, 0.292420, 0.066667}, {0.000000, 0.292420, 0.000000}},
								{{0.000000, -0.225753, -0.225753}, {0.000000, 0.112876, 0.000000}, {0.000000, 0.112876, 0.112876}, {0.000000, 0.000000, 0.112876}}
							}
						});
						// clang-format on
					}
				}
			}

			WHEN("neo-hookian Cauchy stress tensor is calculated")
			{
				auto const & b = Derivatives::b(F);
				auto const & sigma = Derivatives::sigma(J, b, lambda, mu);

				THEN("stress is correct")
				{
					// clang-format off
					check_equal(sigma, "sigma", {
						{-0.754518, -0.4, -0.4},
						{-0.4, -0.554518, 0.0},
						{-0.4, 0.0, -0.554518}
					});
					// clang-format on
				}

				AND_WHEN("initial stress component of tangent stiffness matrix is calculated")
				{
					Element::Stiffness const & Ks = Derivatives::Ks(dN_by_dx, v, sigma);

					THEN("stress component is correct")
					{
//						Element::Stiffness const & check = 0;
//						check_equal(Ks, "Ks", check, "zero");
						// clang-format off
						check_equal(Ks, "Ks", {
							{
								{{-1.154518, -0.000000, -0.000000}, {0.159086, 0.000000, 0.000000}, {0.836345, 0.000000, 0.000000}, {0.159086, 0.000000, 0.000000}},
								{{-0.000000, -1.154518, -0.000000}, {0.000000, 0.159086, 0.000000}, {0.000000, 0.836345, 0.000000}, {0.000000, 0.159086, 0.000000}},
								{{-0.000000, -0.000000, -1.154518}, {0.000000, 0.000000, 0.159086}, {0.000000, 0.000000, 0.836345}, {0.000000, 0.000000, 0.159086}}
							}, {
								{{0.159086, 0.000000, 0.000000}, {-0.046210, -0.000000, -0.000000}, {-0.112876, -0.000000, -0.000000}, {0.000000, 0.000000, 0.000000}},
								{{0.000000, 0.159086, 0.000000}, {-0.000000, -0.046210, -0.000000}, {-0.000000, -0.112876, -0.000000}, {0.000000, 0.000000, 0.000000}},
								{{0.000000, 0.000000, 0.159086}, {-0.000000, -0.000000, -0.046210}, {-0.000000, -0.000000, -0.112876}, {0.000000, 0.000000, 0.000000}}
							}, {
								{{0.836345, 0.000000, 0.000000}, {-0.112876, -0.000000, -0.000000}, {-0.610592, -0.000000, -0.000000}, {-0.112876, -0.000000, -0.000000}},
								{{0.000000, 0.836345, 0.000000}, {-0.000000, -0.112876, -0.000000}, {-0.000000, -0.610592, -0.000000}, {-0.000000, -0.112876, -0.000000}},
								{{0.000000, 0.000000, 0.836345}, {-0.000000, -0.000000, -0.112876}, {-0.000000, -0.000000, -0.610592}, {-0.000000, -0.000000, -0.112876}}
							}, {
								{{0.159086, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {-0.112876, -0.000000, -0.000000}, {-0.046210, -0.000000, -0.000000}},
								{{0.000000, 0.159086, 0.000000}, {0.000000, 0.000000, 0.000000}, {-0.000000, -0.112876, -0.000000}, {-0.000000, -0.046210, -0.000000}},
								{{0.000000, 0.000000, 0.159086}, {0.000000, 0.000000, 0.000000}, {-0.000000, -0.000000, -0.112876}, {-0.000000, -0.000000, -0.046210}}
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
		using Tensor::Func::all;
		using Tensor::Func::last;
		using Tensor::Func::seq;

		// Material properties: https://www.azom.com/properties.aspx?ArticleID=920
		double  mu = 0.4; // Shear modulus: 0.0003 - 0.02
		double const E = 1;	   // Young's modulus: 0.001 - 0.05
		// Lame's first parameter: https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
		double lambda = (mu * (E - 2 * mu)) / (3 * mu - E);
		mu = lambda = 4;

		auto mesh = load_ovm_mesh(file_name_one);
		Attributes attrib{mesh};
		auto const & vtxhs = attrib.vtxh[0];
		auto const & X = attrib.X.for_element(vtxhs);
		// Push top-most node to the right slightly.
		attrib.x[3](0) += 0.5;
		auto x = attrib.x.for_element(vtxhs);
		auto const & dN_by_dX = attrib.dN_by_dX[0];

		std::stringstream s;
		s << "Lambda = " << lambda << "; mu = " << mu;
		INFO(s.str())
		INFO("Material vertices:")
		INFO(X)

		// Construct gravity force tensor.
//		Scalar constexpr m = 0.1;
//		Node::Forces G = 0;
//		G(all, 1) = -9.81 * m;
//		INFO("G")
//		INFO(G)

		WHEN("element stiffness and internal force tensor attributes are constructed")
		{
			Solver::update_elements_stiffness_and_internal_forces(mesh, attrib, lambda, mu);

			THEN("attributes hold correct solutions")
			{
				Cellh const cellh{0};
				auto const &x = attrib.x.for_element(attrib.vtxh[cellh]);
				Scalar const v = Derivatives::V(x);
				auto const &dN_by_dX = attrib.dN_by_dX[cellh];

				auto const &F = Derivatives::dx_by_dX(x, dN_by_dX);
				auto const FFt = Derivatives::b(F);
				Scalar const J = Derivatives::J(F);

				auto const &dx_by_dL = Derivatives::dX_by_dL(x);
				auto const &dL_by_dx = Derivatives::dL_by_dX(dx_by_dL);
				auto dN_by_dx = Derivatives::dN_by_dX(dL_by_dx);

				auto const &sigma = Derivatives::sigma(J, FFt, lambda, mu);

				auto const &c = Derivatives::c(J, lambda, mu);
				Node::Forces const & T = Derivatives::T(dN_by_dx, v, sigma);

				auto const &Kc = Derivatives::Kc(dN_by_dx, v, c);
				auto const &Ks = Derivatives::Ks(dN_by_dx, v, sigma);
				Element::Stiffness const & K = Kc + Ks;

				check_equal(attrib.T[cellh], "T (attribute)", T, "T (check)");
				check_equal(attrib.K[cellh], "K (attribute)", K, "K (check)");
			}
		}

		WHEN("displacement is solved")
		{
			std::size_t step;
			std::size_t constexpr max_steps = 1000;
			std::string log;

			for (step = 0; step < max_steps; step++)
			{
				log += fmt::format("\n\n>>>>>>>>>>>> Iteration {} <<<<<<<<<<<<\n", step);

				for (auto node_idx : ranges::views::ints(0ul, X.dimension(0)))
				{
					using OvmVec = decltype(mesh)::PointT;
					OvmVec mesh_vtx{x(node_idx, 0), x(node_idx, 1), x(node_idx, 2)};
					mesh.set_vertex(vtxhs[node_idx], mesh_vtx);
				}
				OpenVolumeMesh::IO::FileManager{}.writeFile(
					fmt::format("artefacts/single_tet_solve.{}.ovm", step), mesh);

				log += fmt::format("\nx\n{}", x);

				Scalar const v = Derivatives::V(x);
				log += fmt::format("\nv = {}", v);

				auto const & F = Derivatives::dx_by_dX(x, dN_by_dX);
				auto const FFt = Derivatives::b(F);
				Scalar const J = Derivatives::J(F);

				auto const & dx_by_dL = Derivatives::dX_by_dL(x);
				auto const & dL_by_dx = Derivatives::dL_by_dX(dx_by_dL);
				auto dN_by_dx = Derivatives::dN_by_dX(dL_by_dx);

				auto const & sigma = Derivatives::sigma(J, FFt, lambda, mu);

				auto const & c = Derivatives::c(J, lambda, mu);
				Node::Forces T = Derivatives::T(dN_by_dx, v, sigma);
//				T -= G;
				T(seq(0, 3), seq(0, 3)) = 0;
				log += fmt::format("\nT (constrained)\n{}", T);

				auto const & Kc = Derivatives::Kc(dN_by_dx, v, c);
				auto const & Ks = Derivatives::Ks(dN_by_dx, v, sigma);
				Element::Stiffness K = Kc + Ks;
				// Boundary condition
				K(0, 0, 0, 0) = 10e4;
				K(0, 1, 0, 1) = 10e4;
				K(0, 2, 0, 2) = 10e4;
				K(1, 0, 1, 0) = 10e4;
				K(1, 1, 1, 1) = 10e4;
//				K(1, 2, 1, 2) = 10e4;
//				K(2, 0, 2, 0) = 10e4;
				K(2, 1, 2, 1) = 10e4;
//				K(2, 2, 2, 2) = 10e4;
				log += fmt::format("\nK (constrained)\n{}", K);

				using Displacements = Eigen::Matrix<double, 12, 1>;
				Eigen::Map<Displacements> T_vec{T.data(), 12, 1};

				using StiffnessMatrix = Eigen::Matrix<Scalar, 12, 12>;
				StiffnessMatrix K_mat = Eigen::Map<StiffnessMatrix>{K.data(), 12, 12};
				log += fmt::format("\nK (matrix)\n{}", K_mat);
				REQUIRE(K_mat.determinant() != Approx(0).margin(0.00001));

				Displacements u_vec = K_mat.ldlt().solve(-T_vec);
				log += fmt::format("\nu (vector)\n{}", u_vec);

				Tensor::Map<4, 3> u{u_vec.data()};
				// Boundary condition.
				u(0, all) = 0;
				u(1, 0) = 0;
				u(1, 1) = 0;
				u(2, 1) = 0;
				log += fmt::format("\nu\n{}", u);

				x += u;

				if (norm(u) < 0.00001)
					break;
			}

			INFO(log)

			THEN("volume returns and strain is zero")
			{
				CHECK(step < max_steps);
				CHECK(Derivatives::V(x) == Approx(1.0 / 6));
				// clang-format off
				check_equal(x, "x", {
					{0.000000, 0.000000, 0.000000},
					{0.000000, 0.000000, 1.000000},
					{1.000000, 0.000000, 0.000000},
					{0.000000, 1.000000, 0.000000}
				});
				// clang-format on
			}
		} // WHEN("displacement is solved")

		WHEN("displacement is solved Gauss-Seidel style")
		{
			std::size_t step;
			std::size_t constexpr max_steps = 10;
			std::string log;

			for (step = 0; step < max_steps; step++)
			{
				log += fmt::format("\n\n>>>>>>>>>>>> Iteration {} <<<<<<<<<<<<\n", step);

				for (auto node_idx : ranges::views::ints(0ul, X.dimension(0)))
				{
					using OvmVec = decltype(mesh)::PointT;
					OvmVec mesh_vtx{x(node_idx, 0), x(node_idx, 1), x(node_idx, 2)};
					mesh.set_vertex(vtxhs[node_idx], mesh_vtx);
				}
				OpenVolumeMesh::IO::FileManager{}.writeFile(
					fmt::format("artefacts/single_tet_solve.{}.ovm", step), mesh);

				log += fmt::format("\nx\n{}", x);

				Scalar const v = Derivatives::V(x);
				log += fmt::format("\nv = {}", v);

				auto const &F = Derivatives::dx_by_dX(x, dN_by_dX);
				auto const FFt = Derivatives::b(F);
				Scalar const J = Derivatives::J(F);

				auto const &dx_by_dL = Derivatives::dX_by_dL(x);
				auto const &dL_by_dx = Derivatives::dL_by_dX(dx_by_dL);
				auto dN_by_dx = Derivatives::dN_by_dX(dL_by_dx);

				auto const &sigma = Derivatives::sigma(J, FFt, lambda, mu);

				auto const &c = Derivatives::c(J, lambda, mu);
				Node::Forces T = Derivatives::T(dN_by_dx, v, sigma);
				//				T -= G;
				T(seq(0, 3), seq(0, 3)) = 0;
				log += fmt::format("\nT (constrained)\n{}", T);

				auto const &Kc = Derivatives::Kc(dN_by_dx, v, c);
				auto const &Ks = Derivatives::Ks(dN_by_dx, v, sigma);
				Element::Stiffness K = Kc + Ks;
				// Boundary condition
				K(0, 0, 0, 0) = 10e4;
				K(0, 1, 0, 1) = 10e4;
				K(0, 2, 0, 2) = 10e4;
				K(1, 0, 1, 0) = 10e4;
				K(1, 1, 1, 1) = 10e4;
				//				K(1, 2, 1, 2) = 10e4;
				//				K(2, 0, 2, 0) = 10e4;
				K(2, 1, 2, 1) = 10e4;
				//				K(2, 2, 2, 2) = 10e4;
				log += fmt::format("\nK (constrained)\n{}", K);

				Node::Positions u = 0;

				for (Tensor::Index a = 0; a < u.dimension(0); a++)
				{
					using namespace Tensor;
					using Func::einsum;
					auto u_a = u(a, all);
//					auto const & T_a = T(a, all);
//					auto const & K_a = K(a, all, all, all);
//					auto const & K_a_a = K(a, all, a, all);
					Tensor::Map<3> const & T_a{&T(a, 0)};
					Tensor::Map<3, 4, 3> const & K_a{&K(a, 0, 0, 0)};
					Tensor::Matrix<3> const & K_a_a = K(a, all, a, all);
//					Tensor::Vector<3> const & Ku = einsum<Idxs<i, b, j>, Idxs<b, j>>(K_a, u);
//					Tensor::Vector<3> const & TKu = T_a - Ku;
//					Tensor::Matrix<3> const & invKa = inv(K_a_a);
//					Tensor::Vector<3> delta = invKa % TKu;
//					u_a += delta;
					u_a = inv(K_a_a) % (-T_a - einsum<Idxs<i, b, j>, Idxs<b, j>>(K_a, u));
				}

				// Boundary condition.
				u(0, all) = 0;
				u(1, 0) = 0;
				u(1, 1) = 0;
				u(2, 1) = 0;
				log += fmt::format("\nu\n{}", u);

				x += u;

				if (norm(u) < 0.00001)
					break;
			}

			INFO(log)

			THEN("solution converges to material configuration")
			{
				CHECK(step < max_steps);
				CHECK(Derivatives::V(x) == Approx(1.0 / 6));
				// clang-format off
				check_equal(x, "x", {
					{0.000000, 0.000000, 0.000000},
					{0.000000, 0.000000, 1.000000},
					{1.000000, 0.000000, 0.000000},
					{0.000000, 1.000000, 0.000000}
				});
				// clang-format on
			}
		}
	}
}

#define EIGEN_FASTOR_ALIGN BOOST_PP_CAT(Eigen::Aligned, FASTOR_MEMORY_ALIGNMENT_VALUE)
template <int rows, int cols>
using EigenConstTensorMap = Eigen::Map<
	Eigen::Matrix<Scalar, rows, cols, Eigen::RowMajor> const, EIGEN_FASTOR_ALIGN>;

static auto constexpr const index_of = [](auto const & haystack, auto && needle) {
  auto const & it = std::find(
	  haystack.cbegin(), haystack.cend(), std::forward<decltype(needle)>(needle));
  return std::distance(haystack.cbegin(), it);
};

SCENARIO("Solution of two elements")
{
	GIVEN("tangent stiffness and equivalent node force tensors")
	{
		using namespace FeltElements;
		using Tensor::Func::all;
		using Tensor::Func::last;
		using Tensor::Func::seq;
		using OvmVtx = Mesh::PointT;
		using VerticesMatrix = Eigen::Matrix<
			Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
		using EigenMapOvmVertices = Eigen::Map<VerticesMatrix const>;
		using EigenMapTensorVertices = Eigen::Map<
			VerticesMatrix const, Eigen::Unaligned,
			Eigen::Stride<(FASTOR_MEMORY_ALIGNMENT_VALUE / sizeof(Scalar)), 1>>;
		using FixedDOF = Eigen::Vector3d;
		using FixedDOFs = Eigen::VectorXd;

		// Material properties: https://www.azom.com/properties.aspx?ArticleID=920
		double mu = 0.4;	 // Shear modulus: 0.0003 - 0.02
		double const E = 1;	 // Young's modulus: 0.001 - 0.05
		// Lame's first parameter: https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
		double lambda = (mu * (E - 2 * mu)) / (3 * mu - E);
		mu = lambda = 4;

		auto mesh = load_ovm_mesh(file_name_two);
		Attributes attrib{mesh};

		FixedDOFs fixedDofs{3 * mesh.n_vertices()};
		for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
		{
			OvmVtx const & vtx = mesh.vertex(*itvtxh);
			Eigen::Index vtx_idx = int{*itvtxh};

			if (vtx == OvmVtx{0, 0, 0} || vtx == OvmVtx{1, 0, 0} || vtx == OvmVtx{0, 0, 1})
				fixedDofs.block<3, 1>(3 * vtx_idx, 0) = FixedDOF{1.0, 1.0, 1.0};
			else
				fixedDofs.block<3, 1>(3 * vtx_idx, 0) = FixedDOF{0, 0, 0};
		}

		auto const material_volume =
			Derivatives::V(attrib.x.for_element(attrib.vtxh[0])) +
			Derivatives::V(attrib.x.for_element(attrib.vtxh[1]));

		// Deform vertex.
		attrib.x[0](0) += 0.5;

		auto const initial_deformed_volume =
			Derivatives::V(attrib.x.for_element(attrib.vtxh[0])) +
			Derivatives::V(attrib.x.for_element(attrib.vtxh[1]));

		CAPTURE(lambda, mu, material_volume, initial_deformed_volume);

		auto const rows = static_cast<Eigen::Index>(mesh.n_vertices());
		Eigen::Index constexpr cols = 3;

		INFO("Mesh material vertices")
		EigenMapOvmVertices mat_vtxs{mesh.vertex(0).data(), rows, cols};
		INFO(mat_vtxs)

		INFO("Mesh spatial vertices")
		EigenMapTensorVertices const & mat_x{attrib.x[0].data(), rows, cols};
		INFO(mat_x)

		WHEN("displacement is solved")
		{
			std::size_t step;
			std::size_t constexpr max_steps = 10;
			std::string log;

			for (step = 0; step < max_steps; step++)
			{
				log += fmt::format("\n\n>>>>>>>>>>>> Iteration {} <<<<<<<<<<<<", step);

				Solver::update_elements_stiffness_and_internal_forces(mesh, attrib, lambda, mu);

				Eigen::VectorXd mat_T{3 * mesh.n_vertices()};
				mat_T.setZero();
				for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
				{
					auto const vtx_idx = (*itvtxh).idx();
					for (auto itcellh = mesh.vc_iter(*itvtxh); itcellh.valid(); itcellh++)
					{
						auto const & cell_vtxhs = attrib.vtxh[*itcellh];
						auto const & cell_vtx_idx = index_of(cell_vtxhs, *itvtxh);
						auto const & cell_T = attrib.T[*itcellh];

						mat_T.block<3, 1>(3 * vtx_idx, 0) +=
						    EigenConstTensorMap<4, 3>{cell_T.data()}.block<1, 3>(cell_vtx_idx, 0);
					}
				}
				mat_T.array() *= (Eigen::VectorXd::Ones(fixedDofs.size()) - fixedDofs).array();

				Eigen::MatrixXd mat_K{3 * mesh.n_vertices(), 3 * mesh.n_vertices()};
				mat_K.setZero();

				auto const update_submatrix = [&log, &mat_K, &attrib](
					auto const vtxh_src, auto const vtxh_dst, auto const cellh)
				{
				  auto const & cell_K = attrib.K[cellh];
				  auto const & cell_vtxhs = attrib.vtxh[cellh];
				  auto const cell_a = index_of(cell_vtxhs, vtxh_src);
				  auto const cell_b = index_of(cell_vtxhs, vtxh_dst);
				  auto const a = vtxh_src.idx();
				  auto const b = vtxh_dst.idx();

				  Tensor::Matrix<3> const Kab = cell_K(cell_a, all, cell_b, all);

//				  log += fmt::format("\nK_{},{}\n{}", a, b, Kab);
//				  log += fmt::format("\nK_{},{}\n{} (mapped)", a, b, EigenMap{Kab.data()});
				  mat_K.block<3, 3>(3 * a, 3 * b) += EigenConstTensorMap<3, 3>{Kab.data()};
				};

				for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
				{
					for (auto itcellh = mesh.vc_iter(*itvtxh); itcellh.valid(); itcellh++)
					{
						update_submatrix(*itvtxh, *itvtxh, *itcellh);
					}
				}
				for (auto itheh = mesh.halfedges_begin(); itheh != mesh.halfedges_end(); itheh++)
				{
					auto const & halfedge = mesh.halfedge(*itheh);
					auto const & vtxh_src = halfedge.from_vertex();
					auto const & vtxh_dst = halfedge.to_vertex();
					for (auto itcellh = mesh.hec_iter(*itheh); itcellh.valid(); itcellh++)
					{
						update_submatrix(vtxh_src, vtxh_dst, *itcellh);
					}
				}

				log += fmt::format("\nT\n{}", mat_T);
//				log += fmt::format("\nK (unconstrained)\n{}", mat_K);
				mat_K += 10e5 * fixedDofs.asDiagonal();
				log += fmt::format("\nK (constrained)\n{}", mat_K);

				REQUIRE(mat_K.determinant() != Approx(0).margin(epsilon));

				Eigen::VectorXd mat_u{3 * mesh.n_vertices()};
				mat_u.setZero();
				mat_u = mat_K.ldlt().solve(-mat_T);
//				log += fmt::format("\nu (unconstrained)\n{}", mat_u);
				mat_u.array() *= (Eigen::VectorXd::Ones(fixedDofs.size()) - fixedDofs).array();
				log += fmt::format("\nu (constrained)\n{}", mat_u);

				for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
				{
					auto const vtx_idx = (*itvtxh).idx();
					attrib.x[*itvtxh] += Tensor::Map<3>{mat_u.block<3, 1>(3 * vtx_idx, 0).data()};
				}

				if (mat_u.lpNorm<Eigen::Infinity>() < epsilon)
					break;
			}

			INFO(log)

			THEN("solution converges to deformed mesh")
			{
				WARN(fmt::format("Converged in {} steps", step));
				CHECK(step < max_steps);
				// clang-format off
				check_equal(mat_x, "x", (VerticesMatrix{mat_x.rows(), mat_x.cols()} <<
					0.500000, 0.000000, 0.000000,
					1.000000, 0.000000, 0.000000,
					0.333333, 1.154972, -0.122244,
					0.000000, 0.000000, 1.000000,
					0.333333, 0.619084, 0.518097
				).finished(), "expected");
				// clang-format on

				auto const total_volume =
					Derivatives::V(attrib.x.for_element(attrib.vtxh[0])) +
					Derivatives::V(attrib.x.for_element(attrib.vtxh[1]));
				CHECK(total_volume == Approx(0.1077625528));

			}
		} // WHEN("displacement is solved")

		WHEN("displacement is solved Gauss-Seidel style")
		{
			std::size_t step;
			std::size_t constexpr max_steps = 20;
			std::string log;

			for (step = 0; step < max_steps; step++)
			{
				log += fmt::format("\n\n>>>>>>>>>>>> Iteration {} <<<<<<<<<<<<", step);

				Solver::update_elements_stiffness_and_internal_forces(mesh, attrib, lambda, mu);

				std::vector<Node::Pos> u{};
				u.resize(mesh.n_vertices());
				for (auto & u_a : u)
					u_a.zeros();

				Scalar max_norm = 0;

				for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
				{
					auto const & vtxh_src = *itvtxh;
					auto const a = vtxh_src.idx();
					auto const vtx = mesh.vertex(vtxh_src);
					Scalar fixed = 0.0;
					if (vtx == OvmVtx{0, 0, 0} || vtx == OvmVtx{1, 0, 0} || vtx == OvmVtx{0, 0, 1})
						fixed = 1.0;

					Node::Force Ta = 0;
					Tensor::Matrix<3> Kaa = 0;
					Node::Force Ka_u = 0;
					for (auto itcellh = mesh.vc_iter(*itvtxh); itcellh.valid(); itcellh++)
					{
						auto const & cell_vtxhs = attrib.vtxh[*itcellh];
						auto const & cell_a = index_of(cell_vtxhs, *itvtxh);
						auto const & cell_T = attrib.T[*itcellh];
						auto const & cell_K = attrib.K[*itcellh];

						Ta += cell_T(cell_a, all);
						Kaa += cell_K(cell_a, all, cell_a, all);
					}
					diag(Kaa) += 10e5 * fixed;

					for (auto itheh = mesh.voh_iter(vtxh_src); itheh.valid(); itheh++)
					{
						Tensor::Multi<Node::dim, Node::dim> Kab = 0;
						auto const & halfedge = mesh.halfedge(*itheh);
						auto const & vtxh_dst = halfedge.to_vertex();
						auto const b = vtxh_dst.idx();

						for (auto itcellh = mesh.hec_iter(*itheh); itcellh.valid(); itcellh++)
						{
							auto const & cell_K = attrib.K[*itcellh];
							auto const & cell_vtxhs = attrib.vtxh[*itcellh];
							auto const cell_a = index_of(cell_vtxhs, vtxh_src);
							auto const cell_b = index_of(cell_vtxhs, vtxh_dst);
//							Tensor::Multi<Node::dim, Node::dim> cell_Kab = cell_K(cell_a, all, cell_b, all);
//							log += fmt::format("\nK[{}, {}, {}]{}{} = {}", (*itcellh).idx(), cell_a, cell_b, a, b, cell_Kab);
							Kab += cell_K(cell_a, all, cell_b, all);
						}
//						log += fmt::format("\nK{}{} = {}", a, b, Kab);
						Ka_u += Kab % u[b];
					}
					using Tensor::Func::inv;
					u[a] = inv(Kaa) % (-Ta - Ka_u) * (1.0 - fixed);

//					log += fmt::format("\nT{} = {}", a, Ta);
//					log += fmt::format("\nK{}{} = {}", a, a, Kaa);
//					log += fmt::format("\nK{} * u = {}", a, Ka_u);
					log += fmt::format("\nu[{}] = {}", a, u[a]);
				}

				for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
				{
					auto const a = (*itvtxh).idx();
					attrib.x[a] += u[a];
					max_norm = std::max(norm(u[a]), max_norm);
				}

				if (max_norm < epsilon)
					break;
			}

			INFO(log)

			THEN("solution converges to deformed mesh")
			{
				WARN(fmt::format("Converged in {} steps", step));
				CHECK(step < max_steps);
				check_equal(attrib.x[0], "x0", {0.500000, 0.000000, 0.000000});
				check_equal(attrib.x[1], "x1", {1.000000, 0.000000, 0.000000});
				check_equal(attrib.x[2], "x2", {0.333333, 1.154972, -0.122244});
				check_equal(attrib.x[3], "x3", {0.000000, 0.000000, 1.000000});
				check_equal(attrib.x[4], "x4", {0.333333, 0.619084, 0.518097});

				auto const total_volume =
					Derivatives::V(attrib.x.for_element(attrib.vtxh[0])) +
					Derivatives::V(attrib.x.for_element(attrib.vtxh[1]));
				CHECK(total_volume == Approx(0.1077625528));
			}
		}
	}
}
