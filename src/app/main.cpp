#include <fmt/format.h>

#include <FeltElements/Attributes.hpp>
#include <FeltElements/internal/Format.hpp>  // Format vectors
#include <boost/program_options.hpp>
#include <iostream>

#include "vec_semantic.hpp"

constexpr const char * kVersion = "1.0.0";

int main(int argc, char * argv[])
{
	namespace po = boost::program_options;
	using FeltElements::Scalar;
	Scalar E;
	FeltElements::Body::Material material{};
	FeltElements::Body::Forces forces{};

	po::options_description main_desc("Usage");
	po::options_description cmdline_desc("Generic options");
	po::options_description config_desc("Configuration options");
	cmdline_desc.add_options()("help,h", "Show help message and quit.")(
		"version,v", "Show version and quit.");
	config_desc.add_options()(
		"density,r",
		po::value(&material.rho)->value_name("REAL")->required(),
		"Density of the material (kg/m^3)")(
		"youngs-modulus,E",
		po::value(&E)->value_name("REAL")->required(),
		"Young's modulus of the material")(
		"shear-modulus,m",
		po::value(&material.mu)->value_name("REAL")->required(),
		"Shear modulus of the material")(
		"pressure,p",
		po::value(&forces.p)->value_name("REAL")->required(),
		"Uniform pressure (Pa = N/m^2)")(
		"body-force,f",
		po::vec_value<3>(&forces.F_by_m)->required(),
		"Body force per unit mass R^3 (N/kg)");

	main_desc.add(cmdline_desc).add(config_desc);

	try
	{
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, main_desc), vm);

		if (vm.count("help"))
		{
			std::cout << main_desc;
			return 0;
		}
		if (vm.count("version"))
		{
			std::cout << kVersion << std::endl;
			return 0;
		}
		po::notify(vm);
	}
	catch (po::error & ex)
	{
		fmt::print("\nError: {}.\n\n", ex.what());
		std::cout << main_desc;
		return 1;
	}

	material.lambda = FeltElements::Body::Material::lames_first(E, material.mu);

	fmt::print(
		"rho = {:.5}\nE = {:.5}\nlambda = {:.5}\np = {:.5}\nf/m = {}\n",
		material.rho,
		E,
		material.lambda,
		forces.p,
		forces.F_by_m);

	return 0;
}
