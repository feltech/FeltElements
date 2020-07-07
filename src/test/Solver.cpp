#include <FeltElements/Attributes.hpp>
#include <FeltElements/Derivatives.hpp>
#include <FeltElements/Solver.hpp>
// clang-format off
#include "util/Format.hpp"
#include <catch2/catch.hpp>
// clang-format on
#include <range/v3/view/iota.hpp>

#include "util/Assert.hpp"
#include "util/IO.hpp"

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
		auto mesh = Test::load_ovm_mesh(file_name_two);

		WHEN("body force attribute is constructed")
		{
			Attribute::Body::Force const attrib_body_force{mesh};

			THEN("it is zero initialised")
			{
				Node::Force const & F = *attrib_body_force;
				CHECK(Tensor::Func::all_of(F == 0));
			}
		}

		WHEN("material properties attribute is constructed")
		{
			Attribute::Body::Properties const attrib_material_properties{mesh};

			THEN("it is zero initialised")
			{
				Attribute::MaterialProperties const & M = *attrib_material_properties;
				CHECK(M.rho == 0);
				CHECK(M.lambda == 0);
				CHECK(M.mu == 0);
			}
		}

		WHEN("element nodal force attributes are constructed")
		{
			Attribute::Cell::NodalForces const attrib_forces{mesh};

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
			Attribute::Cell::Stiffness const attrib_stiffness{mesh};

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
			FeltElements::Attribute::Cell::VertexHandles const attrib_vtxhs{mesh};

			THEN("properties are populated")
			{
				auto itcellh = mesh.cells_begin();
				CHECK(
					attrib_vtxhs[*itcellh] ==
					std::array<Vtxh, 4>{Vtxh{0}, Vtxh{3}, Vtxh{1}, Vtxh{4}});
				itcellh++;
				CHECK(
					attrib_vtxhs[*itcellh] ==
					std::array<Vtxh, 4>{Vtxh{2}, Vtxh{0}, Vtxh{1}, Vtxh{4}});
			}

			AND_WHEN("element natural wrt material coords attributes are constructed")
			{
				Attribute::Vertex::MaterialPosition const attrib_X{mesh};
				FeltElements::Attribute::Cell::MaterialShapeDerivative const attrib_dN_by_dX{
					mesh, attrib_vtxhs, attrib_X};

				THEN("properties are populated")
				{
					auto itcellh = mesh.cells_begin();
					Test::check_equal(
						attrib_dN_by_dX[*itcellh],
						std::string{"dN_by_dX "} + std::to_string(itcellh),
						Derivatives::dN_by_dX(attrib_X.for_element(attrib_vtxhs[*itcellh])));
					itcellh++;
					Test::check_equal(
						attrib_dN_by_dX[*itcellh],
						std::string{"dN_by_dX "} + std::to_string(itcellh),
						Derivatives::dN_by_dX(attrib_X.for_element(attrib_vtxhs[*itcellh])));
				}
			}
		}

		WHEN("material position attributes are constructed")
		{
			Attribute::Vertex::MaterialPosition const attrib_x{mesh};

			THEN("attributes are initialised to material position")
			{
				for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
				{
					Tensor::ConstMap<3> const vtx{mesh.vertex(*itvtxh).data()};
					Test::check_equal(attrib_x[*itvtxh], "x", vtx, "X");
				}
			}

			AND_WHEN("element nodal spatial position tensor is constructed")
			{
				Attribute::Cell::VertexHandles const attrib_vtxhs{mesh};

				auto itcellh = mesh.cells_begin();
				Node::Positions x1 = attrib_x.for_element(attrib_vtxhs[*itcellh]);
				itcellh++;
				Node::Positions x2 = attrib_x.for_element(attrib_vtxhs[*itcellh]);

				THEN("positions are expected")
				{
					// clang-format off
					Test::check_equal(x1, "x1", {
						{0.000000, 0.000000, 0.000000},
						{0.000000, 0.000000, 1.000000},
						{1.000000, 0.000000, 0.000000},
						{0.000000, 0.500000, 0.500000}
					});
					Test::check_equal(x2, "x2", {
						{0.000000, 1.000000, 0.000000},
						{0.000000, 0.000000, 0.000000},
						{1.000000, 0.000000, 0.000000},
						{0.000000, 0.500000, 0.500000}
					});
					// clang-format on
				}

				AND_WHEN("element material volume attribute is constructed")
				{
					Attribute::Cell::MaterialVolume const attrib_material_volume{
						mesh, attrib_vtxhs, attrib_x};

					THEN("it is initialised to volume")
					{
						itcellh = mesh.cells_begin();
						Scalar V = attrib_material_volume[*itcellh];
						CHECK(V == 1.0 / 12.0);
						itcellh++;
						V = attrib_material_volume[*itcellh];
						CHECK(V == 1.0 / 12.0);
					}
				}

				AND_WHEN("element spatial volume attribute is constructed")
				{
					Attribute::Cell::SpatialVolume const attrib_spatial_volume{mesh};

					THEN("it is initialised to zero")
					{
						itcellh = mesh.cells_begin();
						Scalar V = attrib_spatial_volume[*itcellh];
						CHECK(V == 0);
						itcellh++;
						V = attrib_spatial_volume[*itcellh];
						CHECK(V == 0);
					}
				}
			}
		}
		WHEN("spatial position attributes are constructed")
		{
			Attribute::Vertex::SpatialPosition const attrib_x{mesh};

			THEN("attributes are initialised to material position")
			{
				for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
				{
					Tensor::ConstMap<3> const vtx{mesh.vertex(*itvtxh).data()};
					Test::check_equal(attrib_x[*itvtxh], "x", vtx, "X");
				}
			}

			AND_WHEN("element nodal spatial position tensor is constructed")
			{
				Attribute::Cell::VertexHandles const attrib_vtxhs{mesh};

				auto itcellh = mesh.cells_begin();
				Node::Positions x1 = attrib_x.for_element(attrib_vtxhs[*itcellh]);
				itcellh++;
				Node::Positions x2 = attrib_x.for_element(attrib_vtxhs[*itcellh]);

				THEN("positions are expected")
				{
					// clang-format off
					Test::check_equal(x1, "x1", {
						{0.000000, 0.000000, 0.000000},
						{0.000000, 0.000000, 1.000000},
						{1.000000, 0.000000, 0.000000},
						{0.000000, 0.500000, 0.500000}
					});
					Test::check_equal(x2, "x2", {
						{0.000000, 1.000000, 0.000000},
						{0.000000, 0.000000, 0.000000},
						{1.000000, 0.000000, 0.000000},
						{0.000000, 0.500000, 0.500000}
					});
					// clang-format on
				}
			}
		}

		WHEN("fixed degree of freedom attributes are constructed")
		{
			Attribute::Vertex::FixedDOF const attrib{mesh};

			THEN("attributes are initialised to zero")
			{
				for (auto itvtxh = mesh.vertices_begin(); itvtxh != mesh.vertices_end(); itvtxh++)
				{
					Node::Pos zero;
					zero.zeros();
					Test::check_equal(attrib[*itvtxh], "fixed DOFs", zero, "zero");
				}
			}
		}
	}
}

SCENARIO("Metrics of undeformed mesh")
{
	GIVEN("one-element mesh")
	{
		FeltElements::Mesh mesh = Test::load_ovm_mesh(file_name_one);

		WHEN("vertex index mapping is fetched")
		{
			Attribute::Cell::VertexHandles const attrib_vtxhs{mesh};
			auto const & vtxhs = attrib_vtxhs[Cellh{0}];

			THEN("mapping is expected")
			{
				CHECK(vtxhs == Element::Vtxhs{Vtxh{0}, Vtxh{1}, Vtxh{2}, Vtxh{3}});
			}

			AND_WHEN("material node position tensor is constructed")
			{
				Attribute::Vertex::MaterialPosition const attrib_X{mesh};
				Node::Positions const X = attrib_X.for_element(vtxhs);

				THEN("expected positions are reported")
				{
					// clang-format off
					Test::check_equal(X, "X", {
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
					Attribute::Vertex::SpatialPosition const attrib_x{mesh};

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
						auto const x = attrib_x.for_element(vtxhs);

						THEN("spatial node positions equal material positions")
						{
							// clang-format off
							Test::check_equal(x, "x", X, "X");
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
		Mesh mesh = Test::load_ovm_mesh(file_name_one);
		Attributes attrib{mesh};

		Cellh cellh{0};
		auto const & vtxh0 = attrib.vtxh[cellh];
		auto const & X = attrib.X.for_element(vtxh0);
		auto x = attrib.x.for_element(vtxh0);

		INFO("Tetrahedron vertices:")
		INFO(X)

		WHEN("a node is deformed")
		{
			x(0, 0) += 0.5;

			THEN("spatial position is updated")
			{
				// clang-format off
				Test::check_equal(x, "x", {
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

			AND_WHEN("spatial volume is calculated from material volume")
			{
				auto const & V = attrib.V[cellh];

				Scalar v = Derivatives::v(V, x);

				THEN("volume is correct")
				{
					CHECK(v == 1.0 / 12.0);
				}

				AND_WHEN("another node is deformed")
				{
					x(1, 2) += 0.1;

					THEN("spatial volume calculated from material volume equals naive calculation")
					{
						CHECK(Derivatives::V(x) == Derivatives::v(V, x));
					}
				}
			}
		}
	}
}

SCENARIO("Coordinate derivatives in undeformed mesh")
{
	THEN("derivative of shape wrt local coords is correct")
	{
		Test::check_equal(
			Derivatives::dN_by_dL,
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

	THEN("derivative of surface shape wrt local coords is correct")
	{
		Test::check_equal(
			Derivatives::dN_by_dS,
			"dN_by_dS",
			{
				// clang-format off
				{-1, -1},
				{1, 0},
				{0, 1},
				// clang-format on
			});
	}

	THEN("derivative of local wrt shape coords is correct")
	{
		Test::check_equal(
			Derivatives::dL_by_dN,
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
		enum
		{
			i,
			j,
			k
		};
		using namespace Tensor;
		using namespace Tensor::Func;
		Matrix<3> const identity =
			einsum<Indices<i, k>, Indices<k, j>>(Derivatives::dL_by_dN, Derivatives::dN_by_dL);

		Test::check_equal(
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

	THEN("determinant tensor of derivative of shape w.r.t local coords is correct")
	{
		Element::ShapeDerivativeDeterminant zero;
		zero.zeros();
		// clang-format off
		Test::check_equal(
			Derivatives::det_dN_by_dL,
			"det_dN_by_dL", {
				{
					{0.000000, 0.000000, 0.000000, 0.000000},
					{0.000000, 0.000000, -1.000000, 1.000000},
					{0.000000, 1.000000, 0.000000, -1.000000},
					{0.000000, -1.000000, 1.000000, 0.000000}
				}, {
					{0.000000, 0.000000, 1.000000, -1.000000},
					{0.000000, 0.000000, 0.000000, 0.000000},
					{-1.000000, 0.000000, 0.000000, 1.000000},
					{1.000000, 0.000000, -1.000000, 0.000000}
				}, {
					{0.000000, -1.000000, 0.000000, 1.000000},
					{1.000000, 0.000000, 0.000000, -1.000000},
					{0.000000, 0.000000, 0.000000, 0.000000},
					{-1.000000, 1.000000, 0.000000, 0.000000}
				}, {
					{0.000000, 1.000000, -1.000000, 0.000000},
					{-1.000000, 0.000000, 1.000000, 0.000000},
					{1.000000, -1.000000, 0.000000, 0.000000},
					{0.000000, 0.000000, 0.000000, 0.000000}
				}
			});
		// clang-format on
	}

	GIVEN("a one-element mesh")
	{
		auto [X, x] = Test::load_tet(file_name_one);
		(void)x;

		WHEN("derivative of material wrt local coords is calculated")
		{
			auto const & dX_by_dL = Derivatives::dX_by_dL(X);

			THEN("derivative is correct")
			{
				// clang-format off
				Test::check_equal(dX_by_dL, "dX_by_dL", {
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
					Test::check_equal(dL_by_dX, "dL_by_dX", {
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
						Test::check_equal(dN_by_dX, "dN_by_dX", {
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

		WHEN("transformation from natural to cartesian coordinates is calculated")
		{
			auto const & N_to_x = Derivatives::N_to_x(X);

			THEN("matrix is correct")
			{
				// clang-format off
				Test::check_equal(N_to_x, "N_to_x", {
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
					Test::check_equal(dN_by_dX, "dN_by_dX", {
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
					Test::check_equal(dx_by_dN, "dx_by_dN", {
						{0, 0, 1, 0},
						{0, 0, 0, 1},
						{0, 1, 0, 0}
					});
					// clang-format on
				}
			}
		}  // WHEN("transformation from natural to cartesian coordinates is calculated")

		// TODO: check normals calculated from halffaces vs. winding order
		WHEN("derivative of surface material wrt local coords is calculated")
		{
			using Tensor::Func::all;
			using Tensor::Func::fseq;
			using Tensor::Func::last;

			INFO("X")
			INFO(X)

			Node::SurfacePositions const & Xs1 = X(fseq<1, last>(), all);

			Node::SurfacePositions Xs2;
			Xs2(0, all) = X(0, all);
			Xs2(fseq<1, last>(), all) = X(fseq<2, last>(), all);

			Node::SurfacePositions Xs3;
			Xs3(fseq<0, 2>(), all) = X(fseq<0, 2>(), all);
			Xs3(fseq<2, last>(), all) = X(fseq<3, last>(), all);

			Node::SurfacePositions Xs4;
			Xs4 = X(fseq<0, 3>(), all);

			auto const & dX1_by_dS = Derivatives::dX_by_dS(Xs1);
			auto const & dX2_by_dS = Derivatives::dX_by_dS(Xs2);
			auto const & dX3_by_dS = Derivatives::dX_by_dS(Xs3);
			auto const & dX4_by_dS = Derivatives::dX_by_dS(Xs4);

			THEN("derivative is correct")
			{
				// clang-format off
				Test::check_equal(Xs1, "Xs1", {
					{0.000000, 0.000000, 1.000000},
					{1.000000, 0.000000, 0.000000},
					{0.000000, 1.000000, 0.000000}
				});
				Test::check_equal(Xs2, "Xs2", {
					{0.000000, 0.000000, 0.000000},
					{1.000000, 0.000000, 0.000000},
					{0.000000, 1.000000, 0.000000}
				});
				Test::check_equal(Xs3, "Xs3", {
					{0.000000, 0.000000, 0.000000},
					{0.000000, 0.000000, 1.000000},
					{0.000000, 1.000000, 0.000000}
				});
				Test::check_equal(Xs4, "Xs4", {
					{0.000000, 0.000000, 0.000000},
					{0.000000, 0.000000, 1.000000},
					{1.000000, 0.000000, 0.000000}
				});
				Test::check_equal(dX1_by_dS, "dX1_by_dS", {
					{1.000000, 0.000000},
					{0.000000, 1.000000},
					{-1.000000, -1.000000}
				});
				Test::check_equal(dX2_by_dS, "dX2_by_dS", {
					{1.000000, 0.000000},
					{0.000000, 1.000000},
					{0.000000, 0.000000}
				});
				Test::check_equal(dX3_by_dS, "dX3_by_dS", {
					{0.000000, 0.000000},
					{0.000000, 1.000000},
					{1.000000, 0.000000}
				});
				Test::check_equal(dX4_by_dS, "dX4_by_dS", {
					{0.000000, 1.000000},
					{0.000000, 0.000000},
					{1.000000, 0.000000}
				});
				// clang-format on
			}

			AND_WHEN("normal is calculated from tangent vectors")
			{
				using Tensor::Func::cross;
				using Tensor::Func::norm;

				Tensor::Vector<3> dX1_by_dS1 = dX1_by_dS(all, 0);
				Tensor::Vector<3> dX1_by_dS2 = dX1_by_dS(all, 1);
				Tensor::Vector<3> na1 = cross(dX1_by_dS1, dX1_by_dS2);
				Tensor::Vector<3> n1 = na1 / norm(na1);

				Tensor::Vector<3> dX2_by_dS1 = dX2_by_dS(all, 0);
				Tensor::Vector<3> dX2_by_dS2 = dX2_by_dS(all, 1);
				Tensor::Vector<3> na2 = cross(dX2_by_dS1, dX2_by_dS2);
				Tensor::Vector<3> n2 = na2 / norm(na2);

				Tensor::Vector<3> dX3_by_dS1 = dX3_by_dS(all, 0);
				Tensor::Vector<3> dX3_by_dS2 = dX3_by_dS(all, 1);
				Tensor::Vector<3> na3 = cross(dX3_by_dS1, dX3_by_dS2);
				Tensor::Vector<3> n3 = na3 / norm(na3);

				Tensor::Vector<3> dX4_by_dS1 = dX4_by_dS(all, 0);
				Tensor::Vector<3> dX4_by_dS2 = dX4_by_dS(all, 1);
				Tensor::Vector<3> na4 = cross(dX4_by_dS1, dX4_by_dS2);
				Tensor::Vector<3> n4 = na4 / norm(na4);

				THEN("normal is pointing away from volume")
				{
					Test::check_equal(
						n1, "n1", {1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0)});
					Test::check_equal(n2, "n2", {0.0, 0.0, 1.0});
					Test::check_equal(n3, "n3", {-1.0, 0.0, 0.0});
					Test::check_equal(n4, "n4", {0.0, 1.0, 0.0});
				}
			}
		}

		WHEN("determinant of derivative of material wrt local coords (Jacobian) is calculated")
		{
			Scalar det_dX_by_dL = Derivatives::det_dx_by_dL(X);

			THEN("Jacobian is 1")
			{
				CHECK(det_dX_by_dL == 1.0);
			}
		}
	}

	GIVEN("First element of a two-element mesh")
	{
		auto [X, x] = Test::load_tet(file_name_two);
		(void)x;

		THEN("expected positions are reported")
		{
			// clang-format off
			Node::Positions expected{
			   {0.0, 0.0, 0.0},
			   {0.0, 0.0, 1.0},
			   {1.0, 0.0, 0.0},
			   {0.0, 0.5, 0.5}
			};
			// clang-format on

			INFO("X:")
			INFO(X)
			CHECK(Test::equal(X, expected));
		}

		WHEN("derivative of material wrt local coords is calculated")
		{
			auto const & dX_by_dL = Derivatives::dX_by_dL(X);

			THEN("derivative is correct")
			{
				// clang-format off
				Test::check_equal(dX_by_dL, "dX_by_dL", {
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
					Test::check_equal(dL_by_dX, "dL_by_dX", {
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
					Matrix<3> const identity = einsum<Idxs<i, k>, Idxs<k, j>>(dX_by_dL, dL_by_dX);
					// clang-format off
					Test::check_equal(identity, "dX_by_dL * dL_by_dX", {
						{1, 0, 0},
						{0, 1, 0},
						{0, 0, 1}
					});
					// clang-format on
				}
			}
		}  // WHEN("derivative of material wrt local coords is calculated")

		WHEN(
			"determinant of derivative of material wrt local coords (functional Jacobian)"
			" is calculated")
		{
			Scalar det_dX_by_dL = Derivatives::det_dx_by_dL(X);

			THEN("Jacobian is 1/2")
			{
				CHECK(det_dX_by_dL == 0.5);
			}
		}

		AND_WHEN("derivative of natural wrt material coords is calculated")
		{
			auto const & dN_by_dX = Derivatives::dN_by_dX(X);

			THEN("derivative is correct")
			{
				// clang-format off
				Test::check_equal(dN_by_dX, "dN_by_dX", {
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
				Test::check_equal(N_to_x, "N_to_x", {
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
					Test::check_equal(dN_by_dX, "dN_by_dX", {
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
		auto [X, x] = Test::load_tet(file_name_one);

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
					Test::check_equal(dx_by_dL, "dx_by_dL", {
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
						Test::check_equal(dL_by_dx, "dL_by_dx", {
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
					Test::check_equal(dN_by_dx, "dN_by_dx", {
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
		auto [X, x] = Test::load_tet(file_name_two);

		INFO("Material vertices:")
		INFO(X)

		WHEN("deformation gradient is calculated from natural coordinate derivative")
		{
			auto const & dN_by_dX = Derivatives::dN_by_dX(X);
			auto const & dx_by_dX = Derivatives::dx_by_dX(X, dN_by_dX);

			THEN("gradient is identity")
			{
				// clang-format off
				Test::check_equal(dx_by_dX, "dx_by_dX", {
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
				Test::check_equal(dx_by_dX, "dx_by_dX", {
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
					 Test::check_equal(dx_by_dX, "dx_by_dX", {
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
		auto [X, x] = Test::load_tet(file_name_one);

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
					Test::check_equal(dN_by_dx, "dN_by_dx", {
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
					Test::check_equal(F, "F", {
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
						Test::check_equal(b, "b", {
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
		auto [X, x] = Test::load_tet(file_name_two);

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
					Test::check_equal(dN_by_dx, "dN_by_dx", {
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
					Test::check_equal(F, "F", {
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
						Test::check_equal(b, "b", {
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
		auto [X, x] = Test::load_tet(file_name_one);

		INFO("Material vertices:")
		INFO(X)

		auto const & dN_by_dX = Derivatives::dN_by_dX(X);
		//		Tetrahedron::Scalar const lambda = 1;
		//		Tetrahedron::Scalar const mu = 1;
		// Material properties: https://www.azom.com/properties.aspx?ArticleID=920
		double const mu = 0.4;	// Shear modulus: 0.0003 - 0.02
		double const E = 1;		// Young's modulus: 0.001 - 0.05
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
				Test::check_equal(sigma, "sigma", {
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
					Test::check_equal(T, "T", {
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
					Test::check_equal(sigma, "sigma", {
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
						// clang-format off
						Test::check_equal(T, "T", {
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
		double const mu = 0.4;	// Shear modulus: 0.0003 - 0.02
		double const E = 1;		// Young's modulus: 0.001 - 0.05
		// Lame's first parameter: https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
		double lambda = (mu * (E - 2 * mu)) / (3 * mu - E);

		std::stringstream s;
		s << "Lambda = " << lambda << "; mu = " << mu;
		INFO(s.str())

		AND_GIVEN("an undeformed tetrahedron")
		{
			auto [X, x] = Test::load_tet(file_name_one);

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
					Test::check_equal(c, "c",   {
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
						// clang-format off
						Test::check_equal(Kc, "Kc", {
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
						Test::check_equal(Ks, "Ks", expected, "zero");
					}
				}
			}
		}  // End AND_GIVEN("an undeformed tetrahedron")

		AND_GIVEN("a deformed tetrahedron")
		{
			auto [X, x] = Test::load_tet(file_name_one);
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
					Test::check_equal(c, "c",   {
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
						// clang-format off
						Test::check_equal(Kc, "Kc", {
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
					Test::check_equal(sigma, "sigma", {
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
						// clang-format off
						Test::check_equal(Ks, "Ks", {
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

		}  // End AND_GIVEN("an undeformed tetrahedron")
	}
}

SCENARIO("Solution of a single element")
{
	GIVEN("single element mesh and material properties")
	{
		using namespace FeltElements;
		using Tensor::Func::all;
		using Tensor::Func::last;
		using Tensor::Func::seq;

		// Material properties: https://www.azom.com/properties.aspx?ArticleID=920
		double mu = 0.4;	 // Shear modulus: 0.0003 - 0.02
		double const E = 1;	 // Young's modulus: 0.001 - 0.05
		// Lame's first parameter: https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
		double lambda = (mu * (E - 2 * mu)) / (3 * mu - E);
		mu = lambda = 4;

		auto mesh = Test::load_ovm_mesh(file_name_one);
		Attributes attrib{mesh};
		(*attrib.material).lambda = lambda;
		(*attrib.material).mu = mu;
		auto const & vtxhs = attrib.vtxh[Cellh{0}];
		auto const & X = attrib.X.for_element(vtxhs);
		// Push top-most node to the right slightly.
		attrib.x[Vtxh{3}](0) += 0.5;

		std::stringstream s;
		s << "Lambda = " << lambda << "; mu = " << mu;
		INFO(s.str())
		INFO("Material vertices:")
		INFO(X)

		WHEN("element stiffness and internal force tensor attributes are constructed")
		{
			Solver::update_elements_stiffness_and_forces(mesh, attrib);

			THEN("attributes hold correct solutions")
			{
				Cellh const cellh{0};
				auto const & x = attrib.x.for_element(attrib.vtxh[cellh]);
				Scalar const v = Derivatives::V(x);
				auto const & dN_by_dX = attrib.dN_by_dX[cellh];

				auto const & F = Derivatives::dx_by_dX(x, dN_by_dX);
				auto const FFt = Derivatives::b(F);
				Scalar const J = Derivatives::J(F);

				auto const & dx_by_dL = Derivatives::dX_by_dL(x);
				auto const & dL_by_dx = Derivatives::dL_by_dX(dx_by_dL);
				auto dN_by_dx = Derivatives::dN_by_dX(dL_by_dx);

				auto const & sigma = Derivatives::sigma(J, FFt, lambda, mu);

				auto const & c = Derivatives::c(J, lambda, mu);
				Node::Forces const & T = Derivatives::T(dN_by_dx, v, sigma);

				auto const & Kc = Derivatives::Kc(dN_by_dx, v, c);
				auto const & Ks = Derivatives::Ks(dN_by_dx, v, sigma);
				Element::Stiffness const & K = Kc + Ks;

				CHECK(attrib.v[cellh] == v);
				Test::check_equal(attrib.T[cellh], "T (attribute)", T, "T (check)");
				Test::check_equal(attrib.K[cellh], "K (attribute)", K, "K (check)");
			}
		}

		WHEN("displacement is solved")
		{
			std::size_t step;
			std::size_t constexpr max_steps = 10;
			std::string log;

			auto x = attrib.x.for_element(vtxhs);
			auto const & dN_by_dX = attrib.dN_by_dX[Cellh{0}];

			for (step = 0; step < max_steps; step++)
			{
				log += fmt::format("\n\n>>>>>>>>>>>> Iteration {} <<<<<<<<<<<<\n", step);
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
				K(2, 1, 2, 1) = 10e4;
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

			THEN("solution converges to material configuration")
			{
				CHECK(step < max_steps);
				CHECK(Derivatives::V(x) == Approx(1.0 / 6));
				// clang-format off
				Test::check_equal(x, "x", {
					{0.000000, 0.000000, 0.000000},
					{0.000000, 0.000000, 1.000000},
					{1.000000, 0.000000, 0.000000},
					{0.000000, 1.000000, 0.000000}
				});
				// clang-format on
			}
		}  // WHEN("displacement is solved")

		WHEN("displacement is solved Gauss-Seidel style")
		{
			std::size_t step;
			std::size_t constexpr max_steps = 10;
			std::string log;

			auto x = attrib.x.for_element(vtxhs);
			auto const & dN_by_dX = attrib.dN_by_dX[Cellh{0}];

			for (step = 0; step < max_steps; step++)
			{
				log += fmt::format("\n\n>>>>>>>>>>>> Iteration {} <<<<<<<<<<<<\n", step);
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

				Node::Positions u = 0;

				for (Tensor::Index a = 0; a < u.dimension(0); a++)
				{
					using namespace Tensor;
					using Func::einsum;
					auto u_a = u(a, all);
					Tensor::Map<3> const & T_a{&T(a, 0)};
					Tensor::Map<3, 4, 3> const & K_a{&K(a, 0, 0, 0)};
					Tensor::Matrix<3> const & K_a_a = K(a, all, a, all);
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
				Test::check_equal(x, "x", {
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

SCENARIO("Solution of two elements")
{
	GIVEN("two element mesh and material properties")
	{
		using namespace FeltElements;
		//		// Material properties: https://www.azom.com/properties.aspx?ArticleID=920
		//		Scalar mu = 0.4;	 // Shear modulus: 0.0003 - 0.02
		//		Scalar const E = 1;	 // Young's modulus: 0.001 - 0.05
		//		// Lame's first parameter: https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
		//		Scalar lambda = (mu * (E - 2 * mu)) / (3 * mu - E);

		auto mesh = Test::load_ovm_mesh(file_name_two);
		Attributes attrib{mesh};
		attrib.material->rho = 4;
		attrib.material->lambda = 4;
		attrib.material->mu = 4;
		*attrib.f = Node::Force{0.0, -9.81, 0.0};
		Test::write_ovm_mesh(mesh, attrib.x, "two_elem_initial");

		auto const total_volume = [&attrib]() {
			return Derivatives::V(attrib.x.for_element(attrib.vtxh[Cellh{0}])) +
				Derivatives::V(attrib.x.for_element(attrib.vtxh[Cellh{1}]));
		};

		auto const material_volume = total_volume();

		// Set boundary condition.
		for (auto vtxh : boost::make_iterator_range(mesh.vertices()))
		{
			using OvmVtx = Mesh::PointT;
			OvmVtx const & vtx = mesh.vertex(vtxh);

			if (vtx == OvmVtx{0, 0, 0} || vtx == OvmVtx{1, 0, 0} || vtx == OvmVtx{0, 0, 1})
				attrib.fixed_dof[vtxh] = Node::Pos{1.0, 1.0, 1.0};
		}
		// Set initial condition.
		attrib.x[Vtxh{0}](0) += 0.5;
		auto const initial_deformed_volume = total_volume();

		CAPTURE(
			attrib.material->rho,
			attrib.material->lambda,
			attrib.material->mu,
			material_volume,
			initial_deformed_volume);

		auto const rows = static_cast<Eigen::Index>(mesh.n_vertices());
		Eigen::Index constexpr cols = 3;

		INFO("Mesh material vertices")
		Solver::LDLT::EigenMapOvmVertices mat_vtxs{mesh.vertex(Vtxh{0}).data(), rows, cols};
		INFO(mat_vtxs)

		INFO("Mesh spatial vertices")
		Solver::LDLT::EigenMapTensorVertices const & mat_x{attrib.x[Vtxh{0}].data(), rows, cols};
		INFO(mat_x)

		auto const check_converges = [&mat_x, &total_volume](
										 auto const max_step, auto const final_step) {
			WARN(fmt::format("Converged in {} steps", final_step));
			CHECK(final_step < max_step);
			// clang-format off
			Test::check_equal(mat_x, "x", (
				Solver::LDLT::VerticesMatrix{mat_x.rows(), mat_x.cols()} <<
					0.5,             0,               0,
					1,               0,               0,
					0.333333333333,  0.627935241331, -0.460080182332,
					0,               0,               1,
					0.333333333333,  0.493114858533,  0.412891883689
				).finished(), "expected");
			// clang-format on

			CHECK(total_volume() == Approx(0.0816044878));
		};

		WHEN("displacement is solved using Eigen LDLT")
		{
			std::size_t constexpr max_steps = 13;
			size_t const step = Solver::LDLT::solve(mesh, attrib, max_steps);

			Test::write_ovm_mesh(mesh, attrib.x, fmt::format("two_elem_ldlt_{}", step));

			THEN("solution converges to deformed mesh")
			{
				WARN(fmt::format("Converged in {} steps", step));
				check_converges(max_steps, step);
			}
		}  // WHEN("displacement is solved")

		WHEN("displacement is solved Gauss-Seidel style")
		{
			std::size_t constexpr max_steps = 63;
			std::size_t step = Solver::Gauss::solve(mesh, attrib, max_steps);

			Test::write_ovm_mesh(mesh, attrib.x, fmt::format("two_elem_gauss_{}", step));

			THEN("solution converges to deformed mesh")
			{
				WARN(fmt::format("Converged in {} steps", step));
				check_converges(max_steps, step);
			}
		}
	}
}
