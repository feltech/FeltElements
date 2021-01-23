#include <FeltElements/Attributes.hpp>
#include <FeltElements/Derivatives.hpp>
#include <FeltElements/Solver.hpp>
// clang-format off
#include "FeltElements/internal/Format.hpp"
#include <catch2/catch.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/numeric.hpp>
// clang-format on
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

		WHEN("material properties attribute is constructed")
		{
			Attribute::MeshBody::MaterialProperties const attrib_material_properties{mesh};

			THEN("it is zero initialised")
			{
				Body::Material const & M = *attrib_material_properties;
				CHECK(M.rho == 0);
				CHECK(M.lambda == 0);
				CHECK(M.mu == 0);
			}
		}
		WHEN("body forces attribute is constructed")
		{
			Attribute::MeshBody::Forces const attrib_forces{mesh};

			THEN("it is zero initialised")
			{
				Body::Forces const & M = *attrib_forces;
				CHECK(M.p == 0);
				using Tensor::Func::all_of;
				CHECK(all_of(M.F_by_m == 0));
			}
		}

		WHEN("surface traction attribute is constructed")
		{
			Attribute::Surface::Traction attrib_traction{mesh};

			THEN("attributes are initialised to zero")
			{
				for (auto halffaceh : boost::make_iterator_range(mesh.halffaces()))
				{
					Node::Force zero;
					zero.zeros();
					Test::check_equal(attrib_traction[halffaceh], "traction", zero, "zero");
				}
			}
		}

		WHEN("element nodal force attributes are constructed")
		{
			Attribute::Cell::NodalForces const attrib_forces{mesh};

			THEN("properties are default initialised")
			{
				auto itcellh = mesh.cells_begin();
				Element::Forces T = attrib_forces[*itcellh];
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
				Element::NodePositions x1 = attrib_x.for_element(attrib_vtxhs[*itcellh]);
				itcellh++;
				Element::NodePositions x2 = attrib_x.for_element(attrib_vtxhs[*itcellh]);

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
				Element::NodePositions x1 = attrib_x.for_element(attrib_vtxhs[*itcellh]);
				itcellh++;
				Element::NodePositions x2 = attrib_x.for_element(attrib_vtxhs[*itcellh]);

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

		WHEN("surface elements attribute is constructed")
		{
			Attribute::Cell::VertexHandles const attrib_vtxhs{mesh};
			Attribute::Cell::Boundary const attrib_boundary{mesh, attrib_vtxhs};
			Attribute::Vertex::MaterialPosition const attrib_X{mesh};

			THEN("it is initialised to the surface vertices of the mesh")
			{
				using Triangle = Tensor::Matrix<4, 3>;

				std::vector<Triangle> trianglesFromVtxhs;

				for (auto cellh : boost::make_iterator_range(mesh.cells()))
				{
					auto const & cell_boundary = attrib_boundary[cellh];
					auto const & cell_vtxhs = attrib_vtxhs[cellh];
					// 3 sides of each tetrahedron are boundary elements.
					CHECK(cell_boundary.size() == 3);

					for (auto const & vtxhidxs : cell_boundary)
					{
						using Tensor::Func::all;
						using Tensor::Func::cross;
						using Tensor::Func::fseq;
						using Tensor::Func::norm;

						BoundaryElement::Vtxhs vtxhs;
						vtxhs[0] = cell_vtxhs[vtxhidxs[0]];
						vtxhs[1] = cell_vtxhs[vtxhidxs[1]];
						vtxhs[2] = cell_vtxhs[vtxhidxs[2]];
						BoundaryElement::NodePositions const & x = attrib_X.for_element(vtxhs);
						Triangle triangle = 0;
						Tensor::Vector<3> const vtx0 = x(0, all);
						Tensor::Vector<3> const vtx1 = x(1, all);
						Tensor::Vector<3> const vtx2 = x(2, all);
						Tensor::Vector<3> normal = cross(vtx1 - vtx0, vtx2 - vtx0);
						normal /= norm(normal);

						triangle(fseq<0, 3>(), all) = x;
						triangle(3, all) = normal;
						trianglesFromVtxhs.push_back(triangle);
					}
				}

				std::vector<Triangle> trianglesFromX;

				for (auto cellh : boost::make_iterator_range(mesh.cells()))
				{
					auto const & cell_boundary = attrib_boundary[cellh];
					auto const & cell_vtxhs = attrib_vtxhs[cellh];
					Element::BoundaryNodePositions const xs =
						attrib_X.for_elements(cell_vtxhs, cell_boundary);

					for (auto const & x : xs)
					{
						using Tensor::Func::all;
						using Tensor::Func::cross;
						using Tensor::Func::fseq;
						using Tensor::Func::norm;

						Triangle triangle = 0;
						Tensor::Vector<3> const vtx0 = x(0, all);
						Tensor::Vector<3> const vtx1 = x(1, all);
						Tensor::Vector<3> const vtx2 = x(2, all);
						Tensor::Vector<3> normal = cross(vtx1 - vtx0, vtx2 - vtx0);
						normal /= norm(normal);

						triangle(fseq<0, 3>(), all) = x;
						triangle(3, all) = normal;
						trianglesFromX.push_back(triangle);
					}
				}

				// clang-format off
				std::vector<Triangle> expected{{
					{0.000000, 0.000000, 0.000000},
					{1.000000, 0.000000, 0.000000},
					{0.000000, 0.000000, 1.000000},
					{0.000000, -1.000000, 0.000000}
				}, {
					{0.000000, 0.000000, 0.000000},
					{0.000000, 0.000000, 1.000000},
					{0.000000, 0.500000, 0.500000},
					{-1.000000, 0.000000, 0.000000}
				}, {
					{0.000000, 0.000000, 1.000000},
					{1.000000, 0.000000, 0.000000},
					{0.000000, 0.500000, 0.500000},
					{0.577350, 0.577350, 0.577350}
				}, {
					{0.000000, 1.000000, 0.000000},
					{1.000000, 0.000000, 0.000000},
					{0.000000, 0.000000, 0.000000},
					{0.000000, 0.000000, -1.000000}
				}, {
					{0.000000, 1.000000, 0.000000},
					{0.000000, 0.500000, 0.500000},
					{1.000000, 0.000000, 0.000000},
					{0.577350, 0.577350, 0.577350}
				}, {
					{0.000000, 1.000000, 0.000000},
					{0.000000, 0.000000, 0.000000},
					{0.000000, 0.500000, 0.500000},
					{-1.000000, 0.000000, 0.000000}
				}};
				// clang-format on

				CHECK(trianglesFromVtxhs.size() == expected.size());
				CHECK(trianglesFromX.size() == expected.size());

				for (std::size_t tri_idx = 0; tri_idx < trianglesFromVtxhs.size(); tri_idx++)
					Test::check_equal(
						trianglesFromVtxhs[tri_idx], "from vtxhs", expected[tri_idx], "expected");

				for (std::size_t tri_idx = 0; tri_idx < trianglesFromX.size(); tri_idx++)
					Test::check_equal(
						trianglesFromX[tri_idx], "from x", expected[tri_idx], "expected");
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
		Cellh cellh{0};

		WHEN("vertex index mapping is fetched")
		{
			Attribute::Cell::VertexHandles const attrib_vtxhs{mesh};
			auto const & vtxhs = attrib_vtxhs[cellh];

			THEN("mapping is expected")
			{
				CHECK(vtxhs == Element::Vtxhs{Vtxh{0}, Vtxh{1}, Vtxh{2}, Vtxh{3}});
			}

			AND_WHEN("material node position tensor is constructed")
			{
				Attribute::Vertex::MaterialPosition const attrib_X{mesh};
				Element::NodePositions const X = attrib_X.for_element(vtxhs);

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

				AND_WHEN("material volume is calculated directly")
				{
					Scalar V = Derivatives::V(X);

					THEN("volume is correct")
					{
						CHECK(V == 1.0 / 6.0);
					}

					AND_WHEN("material volume is calculated from local coord transform")
					{
						Scalar v = Derivatives::v(X);

						THEN("volume is correct")
						{
							CHECK(v == 1.0 / 6.0);
						}
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

				AND_WHEN("boundary position tensor is constructed")
				{
					Attribute::Cell::Boundary const attrib_boundary{mesh, attrib_vtxhs};
					Element::BoundaryVtxhIdxs const & vtxh_idxs = attrib_boundary[cellh];
					Element::BoundaryNodePositions const & s =
						attrib_X.for_elements(vtxhs, vtxh_idxs);

					THEN("areas of boundary faces are calculated correctly")
					{
						CHECK(Derivatives::A(s[0]) == 0.5);
						CHECK(Derivatives::A(s[1]) == Approx(0.8660254038));
						CHECK(Derivatives::A(s[2]) == 0.5);
						CHECK(Derivatives::A(s[3]) == 0.5);
					}
				}
			}  // AND_WHEN("material node position tensor is constructed")
		}
	}

	GIVEN("two-element mesh")
	{
		FeltElements::Mesh mesh = Test::load_ovm_mesh(file_name_two);

		AND_WHEN("material volume is calculated directly and by local transform")
		{
			Attribute::Cell::VertexHandles const attrib_vtxhs{mesh};
			Attribute::Vertex::MaterialPosition const attrib_X{mesh};

			auto const & vtxhs = attrib_vtxhs[Cellh{0}];
			Element::NodePositions const X = attrib_X.for_element(vtxhs);

			Scalar V = Derivatives::V(X);
			Scalar v = Derivatives::v(X);

			THEN("volumes agree")
			{
				CHECK(V == v);
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
		auto const & vtxh0 = attrib.vtxhs[cellh];
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

			AND_WHEN("spatial volume is calculated from local transform")
			{
				Scalar v = Derivatives::v(x);

				THEN("volume is correct")
				{
					CHECK(v == 1.0 / 12.0);
				}

				AND_WHEN("another node is deformed")
				{
					x(1, 2) += 0.1;

					THEN("spatial volume calculated from local transform equals naive calculation")
					{
						Scalar V = Derivatives::V(x);
						v = Derivatives::v(x);
						CHECK(V == v);
					}
				}
			}
		}
	}
}

SCENARIO("Coordinate derivatives in undeformed mesh")
{
	THEN("Levi-Civita alternating tensor is correct")
	{
		CHECK(Derivatives::levi_civita(0, 0, 0) == 0);
		CHECK(Derivatives::levi_civita(0, 0, 1) == 0);
		CHECK(Derivatives::levi_civita(0, 0, 2) == 0);
		CHECK(Derivatives::levi_civita(0, 1, 0) == 0);
		CHECK(Derivatives::levi_civita(0, 1, 1) == 0);
		CHECK(Derivatives::levi_civita(0, 1, 2) == 1);
		CHECK(Derivatives::levi_civita(0, 2, 0) == 0);
		CHECK(Derivatives::levi_civita(0, 2, 1) == -1);
		CHECK(Derivatives::levi_civita(0, 2, 2) == 0);
		CHECK(Derivatives::levi_civita(1, 0, 0) == 0);
		CHECK(Derivatives::levi_civita(1, 0, 1) == 0);
		CHECK(Derivatives::levi_civita(1, 0, 2) == -1);
		CHECK(Derivatives::levi_civita(1, 1, 0) == 0);
		CHECK(Derivatives::levi_civita(1, 1, 1) == 0);
		CHECK(Derivatives::levi_civita(1, 1, 2) == 0);
		CHECK(Derivatives::levi_civita(1, 2, 0) == 1);
		CHECK(Derivatives::levi_civita(1, 2, 1) == 0);
		CHECK(Derivatives::levi_civita(1, 2, 2) == 0);
		CHECK(Derivatives::levi_civita(2, 0, 0) == 0);
		CHECK(Derivatives::levi_civita(2, 0, 1) == 1);
		CHECK(Derivatives::levi_civita(2, 0, 2) == 0);
		CHECK(Derivatives::levi_civita(2, 1, 0) == -1);
		CHECK(Derivatives::levi_civita(2, 1, 1) == 0);
		CHECK(Derivatives::levi_civita(2, 1, 2) == 0);
		CHECK(Derivatives::levi_civita(2, 2, 0) == 0);
		CHECK(Derivatives::levi_civita(2, 2, 1) == 0);
		CHECK(Derivatives::levi_civita(2, 2, 2) == 0);
	}

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
		auto mesh = Test::load_ovm_mesh(file_name_one);
		Attributes const attrib{mesh};
		auto const & X = attrib.X.for_element(attrib.vtxhs[Cellh{0}]);

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

			BoundaryElement::NodePositions const & Xs1 =
				attrib.X.for_element(attrib.surface_vtxh->at(0));
			BoundaryElement::NodePositions const & Xs2 =
				attrib.X.for_element(attrib.surface_vtxh->at(1));
			BoundaryElement::NodePositions const & Xs3 =
				attrib.X.for_element(attrib.surface_vtxh->at(2));
			BoundaryElement::NodePositions const & Xs4 =
				attrib.X.for_element(attrib.surface_vtxh->at(3));

			auto const & dX1_by_dS = Derivatives::dX_by_dS(Xs1);
			auto const & dX2_by_dS = Derivatives::dX_by_dS(Xs2);
			auto const & dX3_by_dS = Derivatives::dX_by_dS(Xs3);
			auto const & dX4_by_dS = Derivatives::dX_by_dS(Xs4);

			THEN("derivative is correct")
			{
				CHECK(attrib.surface_vtxh->size() == 4);
				// clang-format off
				Test::check_equal(Xs1, "Xs1", {
					{0.000000, 0.000000, 0.000000},
					{1.000000, 0.000000, 0.000000},
					{0.000000, 0.000000, 1.000000}
				});
				Test::check_equal(Xs2, "Xs2", {
					{0.000000, 0.000000, 1.000000},
					{1.000000, 0.000000, 0.000000},
					{0.000000, 1.000000, 0.000000}
				});
				Test::check_equal(Xs3, "Xs3", {
					{0.000000, 0.000000, 0.000000},
					{0.000000, 1.000000, 0.000000},
					{1.000000, 0.000000, 0.000000}
				});
				Test::check_equal(Xs4, "Xs4", {
					{0.000000, 0.000000, 0.000000},
					{0.000000, 0.000000, 1.000000},
					{0.000000, 1.000000, 0.000000}
				});
				Test::check_equal(dX1_by_dS, "dX1_by_dS", {
					{1.000000, 0.000000},
					{0.000000, 0.000000},
					{0.000000, 1.000000}
				});
				Test::check_equal(dX2_by_dS, "dX2_by_dS", {
					{1.000000, 0.000000},
					{0.000000, 1.000000},
					{-1.000000, -1.000000}
				});
				Test::check_equal(dX3_by_dS, "dX3_by_dS", {
					{0.000000, 1.000000},
					{1.000000, 0.000000},
					{0.000000, 0.000000}
				});
				Test::check_equal(dX4_by_dS, "dX4_by_dS", {
					{0.000000, 0.000000},
					{0.000000, 1.000000},
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
					Test::check_equal(n1, "n1", {0.0, -1.0, 0.0});
					Test::check_equal(
						n2, "n2", {1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0)});
					Test::check_equal(n3, "n3", {0.0, 0.0, -1.0});
					Test::check_equal(n4, "n4", {-1.0, 0.0, 0.0});
				}
			}

			AND_WHEN("n*da is calculated using Levi Cevita tensor")
			{
				using Tensor::Idxs;
				using Tensor::Func::all;
				using Tensor::Func::einsum;
				using Tensor::Func::fix;
				enum
				{
					i,
					j,
					k
				};

				Tensor::Vector<3> n1da = einsum<Idxs<i, j, k>, Idxs<j>, Idxs<k>>(
					Derivatives::levi_civita,
					Tensor::Vector<3>{dX1_by_dS(all, fix<0>)},
					Tensor::Vector<3>{dX1_by_dS(all, fix<1>)});
				Tensor::Vector<3> n2da = einsum<Idxs<i, j, k>, Idxs<j>, Idxs<k>>(
					Derivatives::levi_civita,
					Tensor::Vector<3>{dX2_by_dS(all, fix<0>)},
					Tensor::Vector<3>{dX2_by_dS(all, fix<1>)});
				Tensor::Vector<3> n3da = einsum<Idxs<i, j, k>, Idxs<j>, Idxs<k>>(
					Derivatives::levi_civita,
					Tensor::Vector<3>{dX3_by_dS(all, fix<0>)},
					Tensor::Vector<3>{dX3_by_dS(all, fix<1>)});
				Tensor::Vector<3> n4da = einsum<Idxs<i, j, k>, Idxs<j>, Idxs<k>>(
					Derivatives::levi_civita,
					Tensor::Vector<3>{dX4_by_dS(all, fix<0>)},
					Tensor::Vector<3>{dX4_by_dS(all, fix<1>)});

				THEN("n*da is correct")
				{
					Test::check_equal(n1da, "n1da", {0.0, -1.0, 0.0});
					Test::check_equal(n2da, "n2da", {1.0, 1.0, 1.0});
					Test::check_equal(n3da, "n3da", {0.0, 0.0, -1.0});
					Test::check_equal(n4da, "n4da", {-1.0, 0.0, 0.0});
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
			Element::NodePositions expected{
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
					Scalar const J = Derivatives::det_dx_by_dX(F);

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
					double const J = Derivatives::det_dx_by_dX(F);

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
		auto mesh = Test::load_ovm_mesh(file_name_one);
		auto attrib = Attributes{mesh};
		auto const & vtxhs = attrib.vtxhs[Cellh{0}];
		auto const & X = attrib.X.for_element(vtxhs);
		auto x = attrib.x.for_element(vtxhs);

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
			auto const J = Derivatives::det_dx_by_dX(F);
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

		WHEN("uniform pressure traction force is calculated")
		{
			Scalar constexpr p = 10.0;
			Node::Force t = Derivatives::t(
				p, Derivatives::dX_by_dS(attrib.x.for_element(attrib.surface_vtxh->at(0))));

			THEN("pressure component is as expected")
			{
				Test::check_equal(t, "t", {0.0, -5.0, 0.0});
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
				auto const J = Derivatives::det_dx_by_dX(F);
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
			auto mesh = Test::load_ovm_mesh(file_name_one);
			Attributes attrib{mesh};
			Cellh const cellh{0};
			auto const & vtxhs = attrib.vtxhs[cellh];
			auto const & X = attrib.X.for_element(vtxhs);
			auto const & x = attrib.x.for_element(vtxhs);

			auto const & dN_by_dX = Derivatives::dN_by_dX(X);
			auto const & F = Derivatives::dx_by_dX(x, dN_by_dX);
			auto const J = Derivatives::det_dx_by_dX(F);
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

			AND_GIVEN("surface shape derivative")
			{
				Element::BoundaryVtxhIdxs const & S_to_Vs = attrib.boundary_faces_vtxh_idxs[cellh];
				Element::BoundaryNodePositions const & Ss = attrib.X.for_elements(vtxhs, S_to_Vs);

				AND_GIVEN("pressure is zero")
				{
					constexpr Scalar p = 0.0;

					WHEN("external pressure tangent stiffness tensor is calculated")
					{
						auto const & Kp = Derivatives::Kp(Ss, S_to_Vs, p);

						THEN("pressure component is zero")
						{
							Test::check_equal(
								Kp,
								"Kp",
								{
									// clang-format off
									{
										{{0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}},
										{{0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}},
										{{0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}}
									}, {
										{{0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}},
										{{0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}},
										{{0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}}
									}, {
										{{0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}},
										{{0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}},
										{{0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}}
									}, {
										{{0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}},
										{{0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}},
										{{0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000}}
									}
									// clang-format on
								});
						}
					}
				}

				AND_GIVEN("pressure is one")
				{
					constexpr Scalar p = 1.0;

					WHEN("external pressure tangent stiffness tensor is calculated")
					{
						auto const & Kp = Derivatives::Kp(Ss, S_to_Vs, p);

						THEN("pressure component is correct")
						{
							Test::check_equal(
								Kp,
								"Kp",
								{
									// clang-format off
									{
										{{0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, -0.166667}, {0.000000, 0.166667, 0.166667}, {0.000000, -0.166667, 0.000000}},
										{{0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, -0.166667}, {-0.166667, 0.000000, 0.000000}, {0.166667, 0.000000, 0.166667}},
										{{0.000000, 0.000000, 0.000000}, {0.166667, 0.166667, 0.000000}, {-0.166667, 0.000000, 0.000000}, {0.000000, -0.166667, 0.000000}}
									}, {
										{{0.000000, 0.000000, 0.166667}, {0.000000, 0.000000, 0.000000}, {0.000000, -0.061004, -0.288675}, {0.000000, 0.061004, -0.061004}},
										{{0.000000, 0.000000, 0.166667}, {0.000000, 0.000000, 0.000000}, {0.061004, 0.000000, -0.061004}, {-0.061004, 0.000000, -0.288675}},
										{{-0.166667, -0.166667, 0.000000}, {0.000000, 0.000000, 0.000000}, {0.288675, 0.061004, 0.000000}, {0.061004, 0.288675, 0.000000}}
									}, {
										{{0.000000, -0.166667, -0.166667}, {0.000000, 0.061004, 0.288675}, {0.000000, 0.000000, 0.000000}, {0.000000, 0.288675, 0.061004}},
										{{0.166667, 0.000000, 0.000000}, {-0.061004, 0.000000, 0.061004}, {0.000000, 0.000000, 0.000000}, {-0.288675, 0.000000, -0.061004}},
										{{0.166667, 0.000000, 0.000000}, {-0.288675, -0.061004, 0.000000}, {0.000000, 0.000000, 0.000000}, {-0.061004, 0.061004, 0.000000}}
									}, {
										{{0.000000, 0.166667, 0.000000}, {0.000000, -0.061004, 0.061004}, {0.000000, -0.288675, -0.061004}, {0.000000, 0.000000, 0.000000}},
										{{-0.166667, 0.000000, -0.166667}, {0.061004, 0.000000, 0.288675}, {0.288675, 0.000000, 0.061004}, {0.000000, 0.000000, 0.000000}},
										{{0.000000, 0.166667, 0.000000}, {-0.061004, -0.288675, 0.000000}, {0.061004, -0.061004, 0.000000}, {0.000000, 0.000000, 0.000000}}
									}
									// clang-format on
								});
						}
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
			auto const J = Derivatives::det_dx_by_dX(F);
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

			WHEN("neo-Hookean Cauchy stress tensor is calculated")
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

//SCENARIO("Single step solution parts")
//{
//	GIVEN("one element mesh and basic material properties")
//	{
//		auto mesh = Test::load_ovm_mesh(file_name_one);
//		Attributes attrib{mesh};
//
//		WHEN("pressure stiffness is solved")
//		{
//			Cellh cellh{0};
//			Scalar const p = 1;
//			auto const & cell_vtxhs = attrib.vtxhs[cellh];
//			auto const & boundary_faces_idxs = attrib.boundary_faces_vtxh_idxs[cellh];
//			auto const & boundary_faces_x = attrib.x.for_elements(cell_vtxhs, boundary_faces_idxs);
//
//			auto const & Kp = Derivatives::Kp(boundary_faces_x, boundary_faces_idxs, p);
//
//			Element::Forces R = 0;
//			for (const auto & [s, idxs] :
//				 boost::range::combine(boundary_faces_x, boundary_faces_idxs))
//			{
//				const Node::Force & t = (1.0 / BoundaryElement::num_nodes) *
//					Derivatives::t(p, Derivatives::dX_by_dS(s));
//				for (Tensor::Index const node_idx : idxs) R(node_idx, Tensor::Func::all) -= t;
//			}
//
//			using Stiffness = Eigen::Matrix<Scalar, 12, 12, Eigen::RowMajor>;
//			using StiffnessMap = Eigen::Map<Stiffness const, EIGEN_FASTOR_ALIGN>;
//			using Residual = Eigen::Matrix<Scalar, 12, 1>;
//			using ResidualMap = Eigen::Map<Residual const, EIGEN_FASTOR_ALIGN>;
//			using Displacement = Eigen::Matrix<Scalar, 12, 1>;
//
//			StiffnessMap const map_K{Kp.data(), 12, 12};
//			ResidualMap const map_R{R.data(), 12};
//			Residual const mat_R = -map_R;
//			Stiffness const mat_K = map_K;
//
//			INFO("K")
//			INFO(mat_K)
//			INFO("R")
//			INFO(mat_R)
//
//			CHECK(std::abs(mat_K.determinant()) > 0.00001);
//
//			Displacement const mat_u = mat_K.inverse() * mat_R;
//
//			THEN("displacements are as expected")
//			{
//				Displacement expected;
//				expected << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
//				CHECK(mat_u == expected);
//			}
//		}
//	}
//}

SCENARIO("Solution of a single element")
{
	GIVEN("single element mesh and material properties")
	{
		using namespace FeltElements;
		using Tensor::Func::all;
		using Tensor::Func::last;
		using Tensor::Func::seq;

		// Material properties: https://www.azom.com/properties.aspx?ArticleID=920
		double mu;	// = 0.4;	 // Shear modulus: 0.0003 - 0.02
					//		double const E = 1;	 // Young's modulus: 0.001 - 0.05
		// Lame's first parameter: https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
		double lambda;	// = (mu * (E - 2 * mu)) / (3 * mu - E);
		mu = lambda = 4;

		auto mesh = Test::load_ovm_mesh(file_name_one);
		Attributes attrib{mesh};
		(*attrib.material).lambda = lambda;
		(*attrib.material).mu = mu;
		auto const & vtxhs = attrib.vtxhs[Cellh{0}];
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
			Solver::update_elements_stiffness_and_residual(mesh, attrib);

			THEN("attributes hold correct solutions")
			{
				Cellh const cellh{0};
				auto const & x = attrib.x.for_element(attrib.vtxhs[cellh]);
				Scalar const v = Derivatives::V(x);
				auto const & dN_by_dX = attrib.dN_by_dX[cellh];

				auto const & F = Derivatives::dx_by_dX(x, dN_by_dX);
				auto const FFt = Derivatives::b(F);
				Scalar const J = Derivatives::det_dx_by_dX(F);

				auto const & dx_by_dL = Derivatives::dX_by_dL(x);
				auto const & dL_by_dx = Derivatives::dL_by_dX(dx_by_dL);
				auto dN_by_dx = Derivatives::dN_by_dX(dL_by_dx);

				auto const & sigma = Derivatives::sigma(J, FFt, lambda, mu);

				auto const & c = Derivatives::c(J, lambda, mu);
				Element::Forces const & T = Derivatives::T(dN_by_dx, v, sigma);

				auto const & Kc = Derivatives::Kc(dN_by_dx, v, c);
				auto const & Ks = Derivatives::Ks(dN_by_dx, v, sigma);
				Element::Stiffness const & K = Kc + Ks;

				Test::check_equal(attrib.R[cellh], "T (attribute)", T, "T (check)");
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
				Scalar const J = Derivatives::det_dx_by_dX(F);

				auto const & dx_by_dL = Derivatives::dX_by_dL(x);
				auto const & dL_by_dx = Derivatives::dL_by_dX(dx_by_dL);
				auto dN_by_dx = Derivatives::dN_by_dX(dL_by_dx);

				auto const & sigma = Derivatives::sigma(J, FFt, lambda, mu);

				auto const & c = Derivatives::c(J, lambda, mu);
				Element::Forces T = Derivatives::T(dN_by_dx, v, sigma);
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
				Scalar const J = Derivatives::det_dx_by_dX(F);

				auto const & dx_by_dL = Derivatives::dX_by_dL(x);
				auto const & dL_by_dx = Derivatives::dL_by_dX(dx_by_dL);
				auto dN_by_dx = Derivatives::dN_by_dX(dL_by_dx);

				auto const & sigma = Derivatives::sigma(J, FFt, lambda, mu);

				auto const & c = Derivatives::c(J, lambda, mu);
				Element::Forces T = Derivatives::T(dN_by_dx, v, sigma);
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

				Element::NodePositions u = 0;

				for (Tensor::Index a = 0; a < Element::NodePositions::dimension(0); a++)
				{
					using namespace Tensor;
					using Func::einsum;
					Tensor::Map<3> const & T_a{&T(a, 0)};
					Tensor::Map<3, 4, 3> const & K_a{&K(a, 0, 0, 0)};
					Tensor::Matrix<3> const & K_a_a = K(a, all, a, all);
					u(a, all) = inv(K_a_a) % (-T_a - einsum<Idxs<i, b, j>, Idxs<b, j>>(K_a, u));
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

void check_solvers(
	std::string_view const & file_name_prefix,
	Mesh & mesh,
	Attributes & attrib,
	std::size_t const max_ldlt_steps,
	std::size_t const max_guass_steps,
	Scalar const expected_volume,
	std::vector<Scalar> const & expected_positions)
{
	auto const total_volume = [&attrib, &mesh]() {
		auto const add = [&attrib](auto const total, auto const & cellh) {
			return total + Derivatives::V(attrib.x.for_element(attrib.vtxhs[cellh]));
		};
		return boost::accumulate(boost::make_iterator_range(mesh.cells()), 0.0, add);
	};

	auto const material_volume = total_volume();
	Test::write_ovm_mesh(mesh, attrib.x, fmt::format("{}_initial_deformed", file_name_prefix));

	CAPTURE(attrib.material->rho, attrib.material->lambda, attrib.material->mu, material_volume);

	auto const rows = static_cast<Eigen::Index>(mesh.n_vertices());
	Eigen::Index constexpr cols = 3;

	INFO("Mesh material vertices")
	Solver::LDLT::EigenMapOvmVertices mat_vtxs{mesh.vertex(Vtxh{0}).data(), rows, cols};
	INFO(mat_vtxs)

	INFO("Mesh initial spatial vertices")
	Solver::LDLT::EigenMapTensorVertices const & mat_x{attrib.x[Vtxh{0}].data(), rows, cols};
	INFO(mat_x)

	auto const check_converges = [&mat_x, &total_volume, expected_volume, &expected_positions](
									 auto const max_step, auto const final_step) {
		CHECK(final_step < max_step);
		Test::check_equal(
			mat_x,
			"x",
			Solver::LDLT::EigenMapOvmVertices{
				expected_positions.data(), mat_x.rows(), mat_x.cols()},
			"expected");

		CHECK(total_volume() == Approx(expected_volume).template margin(Test::epsilon));
	};

	WHEN("displacement is solved using Eigen LDLT")
	{
		auto const max_steps = max_ldlt_steps;
		size_t const step = Solver::LDLT::solve(mesh, attrib, max_steps);

		Test::write_ovm_mesh(mesh, attrib.x, fmt::format("{}_ldlt_{}", file_name_prefix, step));

		THEN("solution converges to deformed mesh")
		{
			if (step < max_steps - 1)
				WARN(fmt::format("LDLT converged in {} steps", step));
			check_converges(max_steps, step);
		}
	}  // WHEN("displacement is solved")

	WHEN("displacement is solved Gauss-Seidel style")
	{
		auto const max_steps = max_guass_steps;
		std::size_t step = Solver::Gauss::solve(mesh, attrib, max_steps);

		Test::write_ovm_mesh(mesh, attrib.x, fmt::format("{}_gauss_{}", file_name_prefix, step));

		THEN("solution converges to deformed mesh")
		{
			if (step < max_steps - 1)
				WARN(fmt::format("Gauss-Seidel converged in {} steps", step));
			check_converges(max_steps, step);
		}
	}
}

//SCENARIO("Solution of one element")
//{
//	GIVEN("one element mesh and basic material properties")
//	{
//		auto mesh = Test::load_ovm_mesh(file_name_one);
//		Attributes attrib{mesh};
//		// Silicon rubber material properties: https://www.azom.com/properties.aspx?ArticleID=920
//		constexpr Scalar E = 0.01 * 1e9;   // Young's modulus: 0.001 - 0.05 GPa
//		attrib.material->rho = 2 * 1e3;	   // Density: 1.1 - 2.3 Mg/m3
//		attrib.material->mu = 0.01 * 1e9;  // Shear modulus: 0.0003 - 0.02 GPa
//		// -- Lame's first parameter: https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
//		attrib.material->lambda =
//			(attrib.material->mu * (E - 2 * attrib.material->mu)) / (3 * attrib.material->mu - E);
//
//		AND_GIVEN("no constraints but pressure")
//		{
////			constexpr Scalar atm = 101325;	// Earth atmospheric pressure (Pa = N/m^2)
//			attrib.forces->p = -1;
//
//			check_solvers(
//				"single_no_constraint",
//				mesh,
//				attrib,
//				100,
//				100,
//				1.0 / 6.0,
//				{
//					// clang-format off
//					0,         0,         0,
//					0,         0,         1,
//					1,         0,         0,
//					0,         1,         0
//					// clang-format on
//				});
//		}
//	}
//}

SCENARIO("Solution of two elements")
{
	GIVEN("two element mesh and basic material properties")
	{
		auto mesh = Test::load_ovm_mesh(file_name_two);
		Attributes attrib{mesh};
		// Silicon rubber material properties: https://www.azom.com/properties.aspx?ArticleID=920
		constexpr Scalar E = 0.01 * 1e9;   // Young's modulus: 0.001 - 0.05 GPa
		attrib.material->rho = 2 * 1e3;	   // Density: 1.1 - 2.3 Mg/m3
		attrib.material->mu = 0.01 * 1e9;  // Shear modulus: 0.0003 - 0.02 GPa
		// -- Lame's first parameter: https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
		attrib.material->lambda =
			(attrib.material->mu * (E - 2 * attrib.material->mu)) / (3 * attrib.material->mu - E);
		attrib.forces->p = 0;					  // Normal pressure
		attrib.forces->F_by_m = {0.0, 0.0, 0.0};  // Force per unit mass.
		// Set boundary condition.
		for (auto vtxh : boost::make_iterator_range(mesh.vertices()))
		{
			Vtx const & vtx = mesh.vertex(vtxh);
			if (vtx == Vtx{0, 0, 0} || vtx == Vtx{1, 0, 0} || vtx == Vtx{0, 0, 1})
				attrib.fixed_dof[vtxh] = Node::Pos{1.0, 1.0, 1.0};
		}

		AND_GIVEN("high density material under self-weight")
		{
			attrib.material->rho *= 1000;				// Density.
			attrib.forces->F_by_m = {0.0, -9.81, 0.0};	// Force per unit mass.

			check_solvers(
				"self_weight",
				mesh,
				attrib,
				4,
				18,
				0.1037816713,
				{
					// clang-format off
				0,         0,         0,
				1,         0,         0,
				0,   0.673986, 0.0744648,
				0,         0,         1,
				0,  0.324621,  0.478153
					// clang-format on
				});
		}

		AND_GIVEN("a free vertex is displaced")
		{
			attrib.x[Vtxh{2}] += Node::Pos{0.2, 0.2, 0.2};

			check_solvers(
				"displaced_vertex",
				mesh,
				attrib,
				2,
				2,
				1.0 / 6.0,
				{
					// clang-format off
					0,         0,         0,
					1,         0,         0,
					0,         1,         0,
					0,         0,         1,
					0,       0.5,       0.5
					// clang-format on
				});
		}

		AND_GIVEN("no constraints and no deformation")
		{
			// Set boundary condition.
			for (auto vtxh : boost::make_iterator_range(mesh.vertices()))
				attrib.fixed_dof[vtxh] = Node::Pos{0.0, 0.0, 0.0};

			check_solvers(
				"no_constraint",
				mesh,
				attrib,
				1,
				1,
				1.0 / 6.0,
				{
					// clang-format off
					0,         0,         0,
					1,         0,         0,
					0,         1,         0,
					0,         0,         1,
					0,       0.5,       0.5
					// clang-format on
				});
		}

		AND_GIVEN("atmospheric pressure")
		{
			constexpr Scalar atm = 101325;	// Earth atmospheric pressure (Pa = N/m^2)
			attrib.forces->p = -10 * atm;

			check_solvers(
				"pressure",
				mesh,
				attrib,
				6,
				12,
				0.1484595104,
				{
					// clang-format off
				0,            0,            0,
				1,            0,            0,
				0,            0.929629,     0.0413998,
				0,            0,            1,
				0,            0.457987,     0.485926
					// clang-format on
				});
		}

		//		AND_GIVEN("no constraints but pressure")
		//		{
		//			// Set boundary condition.
		//			for (auto vtxh : boost::make_iterator_range(mesh.vertices()))
		//				attrib.fixed_dof[vtxh] = Node::Pos{0.0, 0.0, 0.0};
		//			constexpr Scalar atm = 101325;	// Earth atmospheric pressure (Pa = N/m^2)
		//			attrib.forces->p = -10 * atm;
		//
		//			check_solvers(
		//				"no_constraint",
		//				mesh,
		//				attrib,
		//				100,
		//				100,
		//				1.0 / 6.0,
		//				{
		//					// clang-format off
		//					0,         0,         0,
		//					1,         0,         0,
		//					0,         1,         0,
		//					0,         0,         1,
		//					0,       0.5,       0.5
		//					// clang-format on
		//				});
		//		}
	}
}
