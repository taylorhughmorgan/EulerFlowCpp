#include "SedovEuler.hpp"

/*************** EULER SOLUTION CLASS ********************************/
EulerSol::EulerSol(pde_state& m_grid,
	bc_map m_BCs,
	size_t m_order,
	std::array<double, 2> m_alpha,
	std::array<double, 2> m_beta,
	double m_gamma)
{
	// Initialize Euler Solution
	gamma = m_gamma;
	order = m_order;
	grid = m_grid;
	alpha = m_alpha;
	beta = m_beta;
	size = grid.size();

	// array of grid ^ order
	grid_to_order.resize(size);
	std::transform(grid.begin(), grid.end(), grid_to_order.begin(), [m_order](double x) { return std::pow(x, m_order); });
}

void EulerSol::createICs(pde_state& rho0, pde_state& v0, pde_state& p0, pde_state& W0)
{
	// Convert primative variables into the initial conditions - W0
	size_t out_size = rho0.size() + v0.size() + p0.size();
	W0.resize(out_size);
	// using equations of state, calculate internal energy
}

void EulerSol::conv2Primatives(pde_state& W, pde_state& rho, pde_state& u, pde_state& E, pde_state& p)
{
	// Convert W result to primative values -> rho, u, E, and P
}

void EulerSol::operator()(const pde_state& x, pde_state& dxdt, const double t)
{
	// Defining the system of equations
	pde_state rho_r(x.begin(), x.begin() + size);					// first block contains rho
	pde_state rho_U_r(x.begin() + size, x.begin() + 2 * size);		// second block contains rho * U
	pde_state rho_E_r(x.begin() + 2 * size, x.begin() + 3 * size);	// third block contains rho * E
	pde_state rho(rho_r.size()), rho_U(rho_U_r.size()), rho_E(rho_E_r.size());

	// dividing by grid ^ order to get rho, rho * U, and rho * E
	std::transform(rho_r.begin(), rho_r.end(), grid_to_order.begin(), rho.begin(), std::divides<double>());
	std::transform(rho_U_r.begin(), rho_U_r.end(), grid_to_order.begin(), rho_U.begin(), std::divides<double>());
	std::transform(rho_E_r.begin(), rho_E_r.end(), grid_to_order.begin(), rho_E.begin(), std::divides<double>());

	// getting primatives u and E by dividing rho and U
	pde_state u(rho.size()), E(rho.size());
}

/*********************** SEDOV BLAST SOLUTION CLASS ***********************/
SedovBlast::SedovBlast(double ScaleLen__m,	// length scale
	double DomainLen__m,		// size of the domain
	double RExpl__m,			// radius of explosion
	double PExpl__Pa,			// pressure of explosion
	double tFinal__s,			// final simulation time
	double rho0__kgpm3 = 1.225, // ambient air density, kg / m ^ 3
	double P0__Pa = 101325,		// ambient air pressure, Pa
	size_t order = 0,			// order of the equations, 0 = cartesian, 1 - cylindrical, 2 = spherical
	double gamma = 1.4,			// ratio of specific heats, N / A
	size_t minNGridPts = 500	// minimum number of grid points
)
{
	//
}

void SedovBlast::solve()
{
	//ODEs = std::make_unique<EulerSol>()
}