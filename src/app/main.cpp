#include <sysexits.h>

#include <atomic>
#include <filesystem>
#include <future>
#include <iostream>

#include <fmt/chrono.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/program_options.hpp>

#include <FeltElements/Attributes.hpp>
#include <FeltElements/MeshDecorator.hpp>
#include <FeltElements/Solver.hpp>
#include <FeltElements/internal/Format.hpp>	 // Format vectors
#include "tap_adaptor.hpp"
#include "vec_semantic.hpp"

constexpr const char * version = "1.0.0";
constexpr const std::chrono::seconds::rep secs_per_dump = 300;

// Requirements for range-v3 to loop over OpenVolumeMesh VertexIter:
// static_assert(ranges::semiregular<OpenVolumeMesh::VertexIter>, "Not semiregular");
// static_assert(ranges::detail::dereferenceable_<OpenVolumeMesh::VertexIter>, "Not derefabley");
// static_assert(ranges::detail::integer_like_<ranges::iter_difference_t<OpenVolumeMesh::VertexIter>>,
// "signed_integer_like_");
// static_assert(ranges::detail::signed_integer_like_<ranges::iter_difference_t<OpenVolumeMesh::VertexIter>>,
// "signed_integer_like_"); static_assert(ranges::weakly_incrementable<OpenVolumeMesh::VertexIter>,
// "weakly_incrementable");
// static_assert(ranges::input_or_output_iterator<OpenVolumeMesh::VertexIter>,
// "input_or_output_iterator"); static_assert(ranges::sentinel_for<OpenVolumeMesh::VertexIter,
// OpenVolumeMesh::VertexIter>, "sized_sentinel_for");

namespace CommandLine
{
constexpr std::string_view kMatrix = "matrix";
constexpr std::string_view kGauss = "gauss";
constexpr std::array Solvers{kMatrix, kGauss};

using ConstraintPlane = FeltElements::Tensor::Vector<4>;

}  // namespace CommandLine

namespace FeltElements
{
void execute(
	std::string const & input_file_path,
	std::string const & output_file_path,
	Body::Material const & material,
	Body::Forces const & forces,
	CommandLine::ConstraintPlane const & constraint_plane,
	std::string_view const solver,
	Solver::Params const params)
{
	spdlog::info("Reading mesh file '{}'", input_file_path);
	Mesh mesh = MeshDecorator::fromFile(input_file_path);

	spdlog::info("Constructing initial mesh attributes");
	Attributes attrs{mesh};
	*attrs.material = material;
	*attrs.forces = forces;

	MeshDecorator interface{mesh, attrs};

	spdlog::info("Calculating fixed vertices");
	auto const norm = constraint_plane(FeltElements::Tensor::Func::fseq<0, 3>());
	auto const height = constraint_plane(3);
	std::size_t num_constrained_vtxs = 0;
	for (auto const & vtxh : interface.vertices)
	{
		Tensor::ConstMap<3> const vec{mesh.vertex(vtxh).data()};
		if (Tensor::Func::inner(norm, vec) > height)
			continue;
		attrs.fixed_dof[vtxh] = {1, 1, 1};
		num_constrained_vtxs++;
	}

	spdlog::info("    constrained {}/{} vertices", num_constrained_vtxs, mesh.n_vertices());

	spdlog::info("Starting {} solver", solver);

	if (solver == CommandLine::kMatrix)
		spdlog::info("Eigen matrix solver using {} threads", Eigen::nbThreads());

	auto const solver_fn = (solver == CommandLine::kMatrix) ? &FeltElements::Solver::Matrix::solve
															: &FeltElements::Solver::Gauss::solve;

	Solver::Stats stats{};

	auto solver_future =
		std::async(std::launch::async, solver_fn, std::ref(mesh), std::ref(attrs), params, &stats);

	std::size_t last_step = 0;
	auto last_dump = std::chrono::system_clock::now();

	while (true)
	{
		if (solver_future.wait_for(std::chrono::milliseconds(1000)) == std::future_status::ready)
			break;

		// Logging.
		std::size_t curr_step = stats.step_counter.load();
		if (curr_step > last_step)
		{
			spdlog::info(
				"Step = {}; increment = {}; delta = {:e}",
				curr_step,
				stats.force_increment_counter.load(),
				stats.max_norm.load());
			last_step = curr_step;
		}

		// Dump intermediate mesh.
		auto now = std::chrono::system_clock::now();
		if (std::chrono::duration_cast<std::chrono::seconds>(now - last_dump).count() >
			secs_per_dump)
		{
			stats.pause.test_and_set(boost::memory_order_acquire);
			std::scoped_lock lock{stats.running};
			now = std::chrono::system_clock::now();
			std::string const file_name = fmt::format(
				"intermediate_{:%Y-%m-%d_%H.%M.%S}_{}.ovm",
				fmt::localtime(std::chrono::system_clock::to_time_t(now)),
				curr_step);
			spdlog::info("Dumping mesh to '{}'", file_name);
			interface.toFile(file_name);
			stats.pause.clear(boost::memory_order_release);
			stats.pause.notify_one();
			last_dump = now;
		}
	}

	try
	{
		Scalar const max_norm = solver_future.get();
		spdlog::info(
			"Finished in {}/{} steps with residual = {}",
			stats.step_counter.load(),
			params.num_steps,
			max_norm);
	}
	catch (std::exception & e)
	{
		spdlog::error("Failed after {} steps: {}", stats.step_counter.load(), e.what());
	}

	spdlog::info("Writing output to '{}'", output_file_path);

	interface.toFile(output_file_path);
}
}  // namespace FeltElements

int main(int argc, char * argv[])
{
	namespace po = boost::program_options;
	using FeltElements::Scalar;

	spdlog::set_default_logger(spdlog::stderr_color_mt("stderr"));

	std::string input_file_path;
	std::string output_file_path;
	std::string solver;
	FeltElements::Solver::Params params;

	FeltElements::Body::Material material{};
	Scalar K;
	FeltElements::Body::Forces forces{};
	CommandLine::ConstraintPlane constraint_plane;

	po::options_description main_desc("Usage: feltelements [input-file] [output-file] [options]");
	po::options_description cmdline_desc("Generic options");
	po::positional_options_description posarg_desc;
	po::options_description config_desc("Configuration options");

	cmdline_desc.add_options()("help,h", "Show help message and quit")(
		"version,v", "Show version and quit");

	posarg_desc.add("mesh-file", 1);
	posarg_desc.add("output-file", 1);

	config_desc.add_options()(
		"input-file,i",
		po::value(&input_file_path)->value_name("STRING")->required(),
		"Mesh file path to load")(
		"output-file,o",
		po::value(&output_file_path)->value_name("STRING")->required(),
		"Output file path")(
		"solver,s",
		po::value(&solver)->value_name("matrix|gauss")->required(),
		"Solver to use: Eigen matrix or Gauss-Seidel")(
		"max-steps,n",
		po::value(&params.num_steps)->value_name("NATURAL")->required(),
		"Maximum number of steps to run the solver")(
		"force-increments,q",
		po::value(&params.num_force_increments)->value_name("NATURAL")->default_value(1),
		"Maximum number of steps to run the solver")(
		"density,r",
		po::value(&material.rho)->value_name("REAL")->required(),
		"Density of the material (kg/m^3)")(
		"bulk-modulus,K",
		po::value(&K)->value_name("REAL")->required(),
		"Bulk modulus of the material")(
		"shear-modulus,m",
		po::value(&material.mu)->value_name("REAL")->required(),
		"Shear modulus of the material")(
		"pressure,p",
		po::value(&forces.p)->value_name("REAL")->required(),
		"Uniform pressure (Pa = N/m^2)")(
		"body-force,f",
		po::vec_value<3>(&forces.F_by_m)->required()->value_name("x y z"),
		"Body force per unit mass R^3 (N/kg)")(
		"constraint-plane,c",
		po::vec_value<4>(&constraint_plane)->required()->value_name("x y z d"),
		"Plane such that vertices on the negative side will be fixed");

	main_desc.add(cmdline_desc).add(config_desc);

	try
	{
		po::variables_map vm;
		po::store(
			po::command_line_parser(argc, argv).options(main_desc).positional(posarg_desc).run(),
			vm);

		if (vm.count("help"))
		{
			std::cout << main_desc;
			return 0;
		}
		if (vm.count("version"))
		{
			std::cout << version << std::endl;
			return 0;
		}

		po::notify(vm);

		if (!boost::algorithm::any_of_equal(CommandLine::Solvers, solver))
			throw po::validation_error{
				po::validation_error::invalid_option_value,
				"solver",
				solver,
				boost::program_options::command_line_style::allow_long};
	}
	catch (po::error & ex)
	{
		fmt::print("\nError: {}.\n\n", ex.what());
		std::cout << main_desc;
		return EX_USAGE;
	}

	material.lambda = FeltElements::Body::Material::lames_first(K, material.mu);

	spdlog::info(
		"Arguments:\nmesh = {}\nsolver = {}\nmax steps = {}\nforce increments = {}\n"
		"rho = {:.5}\nK = {:.5}\nlambda = {:.5}\nmu = {:.5}\np = {:.5}\nf/m = {}\n"
		"constraint plane = {}",
		input_file_path,
		solver,
		params.num_steps,
		params.num_force_increments,
		material.rho,
		K,
		material.lambda,
		material.mu,
		forces.p,
		forces.F_by_m,
		constraint_plane);

	try
	{
		FeltElements::execute(
			input_file_path, output_file_path, material, forces, constraint_plane, solver, params);
	}
	catch (std::filesystem::filesystem_error & e)
	{
		spdlog::error(e.what());
		return e.code().value();
	}

	return EX_OK;
}
