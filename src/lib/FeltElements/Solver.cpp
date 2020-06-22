#include "Solver.hpp"
namespace FeltElements::Solver
{
void update_elements_stiffness_and_internal_forces(
	Mesh const& mesh, FeltElements::Attributes& attributes, Scalar const lambda, Scalar const mu)
{
	for (auto itcellh = mesh.cells_begin(); itcellh != mesh.cells_end(); itcellh++)
	{
		auto const & cellh = *itcellh;
		auto const & cell_vtxhs = attributes.vtxh[cellh];
		auto [K, T] = Derivatives::KT(
			attributes.x.for_element(cell_vtxhs), attributes.dN_by_dX[cellh], lambda, mu);
		attributes.T[*itcellh] = T;
		attributes.K[*itcellh] = K;
	}
}
}
