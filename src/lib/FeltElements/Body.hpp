#pragma once

#include "Typedefs.hpp"

namespace FeltElements::Body
{
struct Material
{
	/// Density.
	Scalar rho;
	/// Shear modulus / Lame's second parameter.
	Scalar mu;
	/// Lame's first parameter.
	Scalar lambda;

	/// Compute Lame's first parameter from bulk modulus and shear modulus.
	static constexpr Scalar lames_first(Scalar const K, Scalar const mu)
	{
		return K - 2 * mu / 3;
	}
};

struct Forces
{
	/// Pressure
	Scalar p;
	/// Body force per unit mass (i.e. acceleration).
	Node::Force F_by_m;
};
}  // namespace FeltElements::Body
