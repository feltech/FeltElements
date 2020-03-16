#define EIGEN_DEFAULT_IO_FORMAT Eigen::IOFormat(3, 0, ", ", ",\n", "", "", "", "")
#include <unsupported/Eigen/CXX11/Tensor>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <boost/range/irange.hpp>

#include <FeltElements/TetGenIO.hpp>
#include <FeltElements/Tetrahedron.hpp>

#define CATCH_CONFIG_MAIN
#include "common.hpp"

char const * file_name_one = "resources/single/tet.1";
char const * file_name_two = "resources/two/tet.1";
char const * file_name_two_quadratic = "resources/double/tet.1";

using namespace FeltElements;


SCENARIO("Loading a tetrahedralisation")
{
//	char cwd[500];
//	getcwd(cwd, 500);
//	std::cerr << "Executing tests in " << cwd << std::endl;

	auto expected_counts = [](const TetGenIO & mesh,
							  auto num_points,
							  auto num_simplexes,
							  auto num_corners) {
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

		AND_WHEN("tetrahedron is loaded from mesh")
		{
			auto tet = Tetrahedron{mesh, 0};
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
		auto const & mesh = TetGenIO{file_name_two_quadratic};
		THEN("expected counts are reported")
		{
			expected_counts(mesh, 14, 2, 10);
		}

		WHEN("tetrahedrons are loaded from mesh")
		{
			auto tet1 = Tetrahedron{mesh, 0};
			auto tet2 = Tetrahedron{mesh, 1};

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


SCENARIO("Coordinate derivatives in undeformed mesh")
{
	THEN("derivative of shape wrt local coords is correct")
	{
		check_equal(Tetrahedron::dN_by_dL, "dN_by_dL", {
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
		check_equal(Tetrahedron::dL_by_dN, "dL_by_dN", {
				// clang-format off
				{0, 1, 0, 0},
				{0, 0, 1, 0},
				{0, 0, 0, 1}
				// clang-format on
		});
	}

	GIVEN("a one-element mesh")
	{
		auto mesh = TetGenIO{file_name_one}.to_mesh();

		WHEN("material node position tensor is constructed")
		{
			Tetrahedron::Node::Positions X = Tetrahedron::X(mesh, 0);

			THEN("expected positions are reported")
			{
				Tetrahedron::Node::Positions expected;
				expected.setValues({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}});

				INFO("X:")
				INFO(X)
				CHECK(equal(X, expected));
			}

			AND_WHEN("derivative of material wrt local coords is calculated")
			{
				auto const & dX_by_dL = Tetrahedron::dX_by_dL(X);

				THEN("derivative is correct")
				{
					// clang-format off
					check_equal(dX_by_dL, "dX_by_dL", {
						{1, 0, 0},
						{0, 1, 0},
						{0, 0, 1}
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
							{1, 0, 0},
							{0, 1, 0},
							{0, 0, 1}
						});
						// clang-format on
					}

					AND_WHEN("derivative of material wrt shape coords is calculated")
					{
						auto const & dN_by_dX = Tetrahedron::dN_by_dX(dL_by_dX);

						THEN("derivative is correct")
						{
							// clang-format off
							check_equal(dN_by_dX, "dN_by_dX", {
								{-1, -1, -1},
								{1, 0, 0},
								{0, 1, 0},
								{0, 0, 1}
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
						{0, 1, 0, 0},
						{0, 0, 1, 0},
						{0, 0, 0, 1}
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
							 {1, 0, 0},
							 {0, 1, 0},
							 {0, 0, 1}
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
							{0, 1, 0, 0},
							{0, 0, 1, 0},
							{0, 0, 0, 1}
						});
						// clang-format on
					}
				}
			}
		}
	}

	GIVEN("First element of a two-element mesh")
	{
		auto mesh = TetGenIO{file_name_two}.to_mesh();

		WHEN("material node position tensor is constructed")
		{
			Tetrahedron::Node::Positions X = Tetrahedron::X(mesh, 0);

			THEN("expected positions are reported")
			{
				Tetrahedron::Node::Positions expected;
				expected.setValues({{0, 0, 0}, {1, 0, 0}, {0, 0, 1}, {0, 0.5, 0.5}});

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
							{1, 0, 0},
							{0, 0, 0.5},
							{0, 1, 0.5}
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
								{1, 0, 0},
								{0, -1, 1},
								{0, 2, 0}
						});
						// clang-format on
					}

					AND_WHEN("derivative of natural wrt material coords is calculated")
					{
						auto const & dN_by_dX = Tetrahedron::dN_by_dX(dL_by_dX);

						THEN("derivative is correct")
						{
							// clang-format off
							check_equal(dN_by_dX, "dN_by_dX", {
									{-1, -1, -1},
									{1, 0, 0},
									{0, -1, 1},
									{0, 2, 0}
							});
							// clang-format on
						}
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
							{1, 0, 0},
							{0, -1, 1},
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
							{0, 1, 0, 0},
							{0, 0, 0, 0.5},
							{0, 0, 1, 0.5}
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
								{1, 0, 0},
								{0, -1, 1},
								{0, 2, 0}
						});
						// clang-format on
					}
				}
			}
		}
	}
}

SCENARIO("Displacements")
{
	GIVEN("a one-element mesh")
	{
		auto mesh = TetGenIO{file_name_one}.to_mesh();
		Tetrahedron::Node::Positions X = Tetrahedron::X(mesh, 0);

		AND_GIVEN("uninitialised displacements")
		{
			Tetrahedron::Node::PosProperty displacements =
				mesh.request_vertex_property<Tetrahedron::Node::Pos>("displacement");
			displacements->set_persistent(true);

			WHEN("displacement tensor is constructed")
			{
				Tetrahedron::Node::Positions u = Tetrahedron::u(mesh, displacements, 0);

				THEN("displacements are zero")
				{
					Tetrahedron::Node::Positions expected;
					expected.setZero();
					check_equal(u, "u", expected, "zero");
				}

				AND_WHEN("spatial node position tensor is calculated")
				{
					Tetrahedron::Node::Positions x = Tetrahedron::x(X, u);

					THEN("spatial position is equal to material position")
					{
						check_equal(X, "X", x, "x");
					}
				}

				AND_WHEN("a node is displaced")
				{
					Tetrahedron::VectorTensor<> delta;
					delta.setValues({0.5, 0, 0});
					u.chip(0, 0) += delta;

					AND_WHEN("spatial node position tensor is calculated")
					{
						Tetrahedron::Node::Positions x = Tetrahedron::x(X, u);

						THEN("spatial position is updated")
						{
							// clang-format off
							check_equal(x, "x", {
								{0.5, 0, 0},
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
	}
}


SCENARIO("Deformation gradient of undeformed element")
{
	GIVEN("first element of two-element mesh")
	{
		auto const & io = TetGenIO{file_name_two};
		auto const & mesh = io.to_mesh();

		auto const & X = Tetrahedron::X(mesh, 0);

		std::stringstream ss;
		INFO("Tetrahedron vertices:")
		INFO(X);

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
	}
}

SCENARIO("Deformation gradient of deformed element")
{
	GIVEN("first element of two-element mesh")
	{
		auto const & io = TetGenIO{file_name_two};
		Tetrahedron tet{io, 0};
		auto mesh = io.to_mesh();

		auto const & X = Tetrahedron::X(mesh, 0);

		std::stringstream ss;
		INFO("Tetrahedron vertices:")
		INFO(X);

		//		WHEN("a node is deformed")
		//		{
		//			// Create displacement per vertex.
		//			Tetrahedron::Node::PosProperty displacements =
		//					mesh.request_vertex_property<Tetrahedron::Node::Pos>("displacement");
		//			displacements->set_persistent(true);
		//
		//			// Fetch displacements of element as a tensor.
		//			Tetrahedron::Node::Positions u = Tetrahedron::u(mesh, displacements, 0);
		//
		//			// Update element displacements.
		//			Tetrahedron::VectorTensor<> delta;
		//			delta.setValues({0.5, 0, 0});
		//			u.chip(0, 0) += delta;
		//
		//			Tetrahedron::Node::Positions x = Tetrahedron::x(X, u);
		//
		//		}
	}
}
//			AND_WHEN("a node is deformed")
//			{
//				tet.u(0) += Eigen::Vector3d{0.5, 0, 0};
//
//				THEN("material coords are unchanged")
//				{
//					CHECK(tet.X(0) == Eigen::Vector3d{0, 0, 0});
//					CHECK(tet.X(1) == Eigen::Vector3d{0, 0, 1});
//					CHECK(tet.X(2) == Eigen::Vector3d{1, 0, 0});
//					CHECK(tet.X(3) == Eigen::Vector3d{0, 0.5, 0.5});
//				}
//
//				THEN("spatial coords equal material with deformation")
//				{
//					CHECK(tet.x(0) == Eigen::Vector3d{0.5, 0, 0});
//					CHECK(tet.x(1) == Eigen::Vector3d{0, 0, 1});
//					CHECK(tet.x(2) == Eigen::Vector3d{1, 0, 0});
//					CHECK(tet.x(3) == Eigen::Vector3d{0, 0.5, 0.5});
//				}
//
//				THEN("spatial volume has changed")
//				{
//					CHECK(tet.v() == tet.V() / 2);
//				}
//
//				AND_WHEN("derivative of shape function wrt material coords is calculated")
//				{
//					auto const & dN_by_dX_deformed = tet.dN_by_dX();
//
//					THEN("derivative is unchanged from before it was deformed")
//					{
//						CHECK(dN_by_dX_deformed == dN_by_dX);
//					}
//				}
//				AND_WHEN("deformation gradient is calculated")
//				{
//					auto const & dx_by_dN = tet.dx_by_dN();
//					auto const & F = tet.dx_by_dX(dx_by_dN, dN_by_dX);
//
//					THEN("deformation gradient is correct")
//					{
//						Eigen::Matrix3d expected;
//						// clang-format off
//						expected << // NOLINT
//							0.5, -0.5, -0.5,
//							0, 1, 0,
//							0, 0, 1;
//						// clang-format on
//						CHECK(F == expected);
//					}
//
//					AND_WHEN("Jacobian of deformation gradient is calculated")
//					{
//						double const J = Tetrahedron::J(F);
//
//						THEN("value is correct")
//						{
//							CHECK(J == 0.5);
//						}
//					}
//
//					AND_WHEN("Left Cauchy-Green / Finger tensor is calculated")
//					{
//						Tetrahedron::GradientMatrix const & b = Tetrahedron::b(F);
//						INFO(b)
//
//						THEN("value is correct")
//						{
//							Tetrahedron::GradientMatrix expected;
//							// clang-format off
//							expected << // NOLINT
//								 0.75, -0.5, -0.5,
//								-0.5,    1,    0,
//								-0.5,    0,    1;
//							// clang-format on
//							CHECK(b.isApprox(expected));
//						}
//					}
//
//					AND_WHEN("derivative of shape function wrt spatial coords is calculated")
//					{
//						auto const & dN_by_dx = tet.dN_by_dx(dN_by_dX, F);
//
//						std::stringstream ss;
//						ss << "dN_by_dX = "  << std::endl << dN_by_dX << std::endl;
//						ss << "F^-1 = "  << std::endl << F.inverse() << std::endl;
//						INFO(ss.str());
//
//						THEN("derivative is correct")
//						{
//							Eigen::Matrix<double, 4, 3> expected;
//							// clang-format off
//							expected << // NOLINT
//								-2, -2, -2,
//								0, -1, 1,
//								2, 1, 1,
//								0, 2, 0;
//							// clang-format on
//							CHECK(dN_by_dx == expected);
//						}
//					}
//				}
//			}
//}

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
			auto const mesh = FeltElements::TetGenIO{file_name_one};
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
					auto const & Kcab = Tetrahedron::Kcab(dN_by_dx, v, c, 0, 1);

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
						dN_by_dx, v, sigma, 0, 1);

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
			auto const mesh = FeltElements::TetGenIO{file_name_one};
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
					expected.setValues({
						// clang-format off
						{
							{{2.15452, 0, 0}, {0, 0.8, 0}, {0, 0, 0.8}},
							{{0, 0.677259, 0}, {0.677259, 0, 0}, {0, 0, 0}},
							{{0, 0, 0.677259}, {0, 0, 0}, {0.677259, 0, 0}}
						}, {
							{{0, 0.677259, 0}, {0.677259, 0, 0}, {0, 0, 0}},
							{{0.8, 0, 0}, {0, 2.15452, 0}, {0, 0, 0.8}},
							{{0, 0, 0}, {0, 0, 0.677259}, {0, 0.677259, 0}}
						}, {
							{{0, 0, 0.677259}, {0, 0, 0}, {0.677259, 0, 0}},
							{{0, 0, 0}, {0, 0, 0.677259}, {0, 0.677259, 0}},
							{{0.8, 0, 0}, {0, 0.8, 0}, {0, 0, 2.15452}}
						}
						// clang-format on
					});
					Eigen::Tensor<bool, 0> comparison = ((c - expected).abs() < 0.00005).all();

					INFO("Expected:")
					INFO(expected);
					INFO("Actual:")
					INFO(c);
					CHECK(comparison(0));
				}

				AND_WHEN("constitutive component of tangent matrix is calculated")
				{
					auto const & Kcab = Tetrahedron::Kcab(dN_by_dx, v, c, 0, 1);

					INFO("Kcab:")
					INFO(Kcab)

					THEN("matrix of values are as expected")
					{
						Eigen::Matrix3d expected;
						expected <<
						// clang-format off
							-0.113, -0.133,      0,
							-0.113, -0.359, -0.113,
							0, -0.133, -0.113;
						// clang-format on
						CHECK(Kcab.isApprox(expected, 0.002));
					}
				}

			} // End WHEN("neo-hookian elasticity tensor is calculated")

			WHEN("neo-hookian Cauchy stress tensor is calculated")
			{
				auto const & sigma = Tetrahedron::neo_hookian_stress(J, b, lambda, mu);

				INFO("sigma:")
				INFO(sigma)

				THEN("Cauchy stress tensor is correct")
				{
					Eigen::Matrix3d expected;
					expected <<
					// clang-format off
						-0.755,   -0.4,   -0.4,
						-0.4, -0.555,      0,
						-0.4,      0, -0.555;
					// clang-format on
					CHECK(sigma.isApprox(expected, 0.002));
				}

				AND_WHEN("initial stress component of tangent stiffness matrix is calculated")
				{
					auto const & Ksab = Tetrahedron::Ksab(dN_by_dx, v, sigma, 0, 1);

					INFO("Ksab:")
					INFO(Ksab)

					THEN("stress component is correct")
					{
						Eigen::Matrix3d expected;
						expected <<
						// clang-format off
							0.159,  0,		0,
							0,		0.159,  0,
							0,		0,		0.159;
						// clang-format on
						CHECK(Ksab.isApprox(expected, 0.001));
					}
				}
			}
		} // End AND_GIVEN("a deformed tetrahedron")
	}
}
