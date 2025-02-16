#include "BoundaryConditions.hpp"


void reflective(pde_state& u, const pde_state& grid, std::vector<size_t> ids, double value) {
	u[ids[0]] = -u[ids[1]];
}

void gradient(pde_state& u, const pde_state& grid, std::vector<size_t> ids, double value) {
	double dx = grid[ids[1]] - grid[ids[0]];
	u[ids[0]] = u[ids[1]] + value * dx;
}

void extrapolated(pde_state& u, const pde_state& grid, std::vector<size_t> ids, double value) {
	u[ids[0]] = 2 * u[ids[1]] - u[ids[2]];
}

void constant_value(pde_state& u, const pde_state& grid, std::vector<size_t> ids, double value) {
	// constant value boundary condition
	u[ids[0]] = value;
}

BCFunction agnosticBCs(const std::vector<size_t>& ids, validBCs bc_id, double bc_value)
{
	// develop boundary conditions agnostic to whether it is applied to the upper or lower bound
	std::function<void(pde_state& u, const pde_state& grid)> bc;

	if (bc_id == validBCs::REFLECTIVE)
	{
		// applying reflective boundary condition
		bc = std::bind(reflective, std::placeholders::_1, std::placeholders::_2, ids, bc_value);
	}
	else if (bc_id == validBCs::GRADIENT) 
	{
		// gradient boundary condition
		bc = std::bind(gradient, std::placeholders::_1, std::placeholders::_2, ids, bc_value);
	}
	else if (bc_id == validBCs::EXTRAPOLATED)
	{
		// extrapolated boundary condition
		bc = std::bind(extrapolated, std::placeholders::_1, std::placeholders::_2, ids, bc_value);
	}
	else if (bc_id == validBCs::CONSTANT)
	{
		// constant value boundary condition
		bc = std::bind(constant_value, std::placeholders::_1, std::placeholders::_2, ids, bc_value);
	}
	else
	{
		throw std::invalid_argument("Boundary Condition has not been implemented. Valid options are: gradient:<value>, extrapolated, constant:<value>.");
	}
	return bc;
}