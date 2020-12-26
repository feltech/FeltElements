#pragma once

#include "Typedefs.hpp"

namespace FeltElements::Material
{
struct Properties
{
	/// Density.
	Scalar rho;
	/// Shear modulus / Lame's second parameter.
	Scalar mu;
	/// Lame's first parameter.
	Scalar lambda;
	/// Pressure
	Scalar p;
	/// Body force per unit mass (i.e. acceleration).
	Node::Force F_by_m;
};
}
