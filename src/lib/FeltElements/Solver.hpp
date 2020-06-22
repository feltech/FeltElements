#pragma once
#include "Attributes.hpp"

namespace FeltElements::Solver
{
void update_elements_stiffness_and_internal_forces(
	Mesh const& mesh, FeltElements::Attributes& attributes, Scalar const lambda, Scalar const mu);
}
