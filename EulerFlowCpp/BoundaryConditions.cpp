#include "BoundaryConditions.hpp"


static void reflective(pde_state& u, const pde_state& grid, ident ids, double value) {
	u[ids[0]] = -u[ids[1]];
}

static void gradient(pde_state& u, const pde_state& grid, ident ids, double value) {
	double dx = grid[ids[1]] - grid[ids[0]];
	u[ids[0]] = u[ids[1]] + value * dx;
}

static void extrapolated(pde_state& u, const pde_state& grid, ident ids, double value) {
	u[ids[0]] = 2 * u[ids[1]] - u[ids[2]];
}

static void constant_value(pde_state& u, const pde_state& grid, ident ids, double value) {
	// constant value boundary condition
	u[ids[0]] = value;
}

BCFunction agnosticBCs(const ident& ids, validBCs bc_id, double bc_value)
{
	// develop boundary conditions agnostic to whether it is applied to the upper or lower bound
	std::function<void(pde_state& u, const pde_state& grid)> bc;

	if (bc_id == validBCs::REFLECTIVE) {
		// applying reflective boundary condition
		bc = std::bind(reflective, std::placeholders::_1, std::placeholders::_2, ids, bc_value);
	}
	else if (bc_id == validBCs::GRADIENT) {
		// gradient boundary condition
		bc = std::bind(gradient, std::placeholders::_1, std::placeholders::_2, ids, bc_value);
	}
	else if (bc_id == validBCs::EXTRAPOLATED) {
		// extrapolated boundary condition
		bc = std::bind(extrapolated, std::placeholders::_1, std::placeholders::_2, ids, bc_value);
	}
	else if (bc_id == validBCs::CONSTANT) {
		// constant value boundary condition
		bc = std::bind(constant_value, std::placeholders::_1, std::placeholders::_2, ids, bc_value);
	}
	else {
		throw std::invalid_argument("Boundary Condition has not been implemented. Valid options are: gradient:<value>, extrapolated, constant:<value>.");
	}
	return bc;
}


void testBoundaryConditions() {
	// creating an empty array to test out boundary conditions
	pde_state u0, grid;
	size_t nEle = 10;
	for (size_t i = 0; i < nEle; i++) {
		u0.push_back(std::pow(double(i), 2));
		grid.push_back(double(i));
	}
	std::array<size_t, 3> left_ids = { 0, 1, 2 };
	std::array<size_t, 3> right_ids = { nEle - 1, nEle - 2, nEle - 3 };

	// using std::functions, create boundary conditions
	BCFunction bc_left = agnosticBCs(left_ids, validBCs::CONSTANT, 10.0);
	BCFunction bc_right = agnosticBCs(right_ids, validBCs::EXTRAPOLATED);

	// print out the grid and u0 prior to applying BCs
	std::cout << "grid (m), u0 (K) " << std::endl;
	for (size_t i = 0; i < nEle; i++) {
		std::cout << grid[i] << ", " << u0[i] << std::endl;
	}

	// print out grid and u0 post-appling bcs
	bc_left(u0, grid);
	bc_right(u0, grid);

	std::cout << "Post-BC application:\ngrid (m), u0 (K) " << std::endl;
	for (size_t i = 0; i < nEle; i++) {
		std::cout << grid[i] << ", " << u0[i] << std::endl;
	}
}