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
#include <FeltElements/MeshFacade.hpp>
#include <FeltElements/Solver.hpp>
#include <FeltElements/internal/Format.hpp>	 // Format vectors
#include "vec_semantic.hpp"

constexpr const char * version = "1.0.0";
constexpr const std::chrono::seconds secs_per_dump{5};

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
struct MeshDumper : MeshIO
{
	using Clock = std::chrono::system_clock;
	using Period = std::chrono::seconds;
	using Time = std::chrono::time_point<Clock>;

	Solver::Base & solver;
	std::string final_path;
	Period period_secs;

	Time last{Clock::now()};

	void maybe_dump(std::size_t const curr_step)
	{
		Time now = Clock::now();
		if (std::chrono::duration_cast<Period>(now - last) < period_secs)
			return;

		Solver::Base::Unpauser pause_lock = solver.pauser.scoped_pause();

		now = Clock::now();
		std::string const file_name = fmt::format(
			"intermediate_{:%Y-%m-%d_%H.%M.%S}_{}.ovm",
			fmt::localtime(Clock::to_time_t(now)),
			curr_step);
		spdlog::info("Dumping mesh to '{}'", file_name);

		toFile(file_name);

		last = now;
	}

	~MeshDumper()
	{
		spdlog::info("Writing final output to '{}'", final_path);
		toFile(final_path);
	}
};

int execute(
	std::string const & input_file_path,
	std::string const & output_file_path,
	Body::Material const & material,
	Body::Forces const & forces,
	CommandLine::ConstraintPlane const & constraint_plane,
	std::string_view const solver_type,
	Solver::Params const params)
{
	spdlog::info("Reading mesh file '{}'", input_file_path);
	Mesh mesh = MeshIO::fromFile(input_file_path);

	spdlog::info("Constructing initial mesh attributes");

	Attributes attrs{mesh};
	*attrs.material = material;
	*attrs.forces = forces;

	spdlog::info("Calculating fixed vertices");

	auto const norm = constraint_plane(Tensor::Func::fseq<0, 3>());
	auto const height = constraint_plane(3);
	std::size_t num_constrained_vtxs = 0;
	// TODO: add to FixedDOF as a method
	for (auto const & vtxh : MeshIters{mesh, attrs}.vertices)
	{
		Tensor::ConstMap<3> const vec{mesh.vertex(vtxh).data()};
		if (Tensor::Func::inner(norm, vec) > height)
			continue;
		attrs.fixed_dof[vtxh] = {1, 1, 1};
		num_constrained_vtxs++;
	}
	spdlog::info("    constrained {}/{} vertices", num_constrained_vtxs, mesh.n_vertices());

	spdlog::info("Starting {} solver", solver_type);

	if (solver_type == CommandLine::kMatrix)
		spdlog::info("Eigen matrix solver using {} threads", Eigen::nbThreads());

	std::unique_ptr<Solver::Base> solver;
	if (solver_type == CommandLine::kMatrix)
		solver = std::make_unique<Solver::Matrix>(mesh, attrs, params);
	else
		solver = std::make_unique<Solver::Gauss>(mesh, attrs, params);

	// RAII dump mesh file even if we die with an exception.
	MeshDumper dumper{{mesh, attrs}, *solver, output_file_path, secs_per_dump};

	auto solver_future = std::async(std::launch::async, &Solver::Base::solve, solver.get());

	std::size_t last_step = 0;

	while (solver_future.wait_for(std::chrono::seconds(1)) != std::future_status::ready)
	{
		// Logging.
		std::size_t const curr_step = solver->stats.step_counter.load();
		if (curr_step > last_step)
		{
			spdlog::info(
				"Step = {}; increment = {}; delta = {:e}",
				curr_step,
				solver->stats.force_increment_counter.load(std::memory_order_relaxed),
				solver->stats.max_norm.load(std::memory_order_relaxed));
			last_step = curr_step;
		}

		// Dump intermediate mesh.
		dumper.maybe_dump(curr_step);
	}

	try
	{
		solver_future.get();  // Will throw any held exception
		spdlog::info(
			"Finished in {}/{} steps with residual = {}",
			solver->stats.step_counter.load(),
			params.num_steps,
			solver->stats.max_norm);
		return EX_OK;
	}
	catch (std::exception & e)
	{
		spdlog::error(
			"Failed after {}/{} steps with residual {}: {}",
			solver->stats.step_counter.load(),
			params.num_steps,
			solver->stats.max_norm.load(),
			e.what());
		return EX_SOFTWARE;
	}
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
		return FeltElements::execute(
			input_file_path, output_file_path, material, forces, constraint_plane, solver, params);
	}
	catch (std::filesystem::filesystem_error & e)
	{
		spdlog::error(e.what());
		return e.code().value();
	}
}
