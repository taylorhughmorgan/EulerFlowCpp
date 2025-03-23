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

	// setting up the ghost grid
	dr = grid[1] - grid[0];
	ghost_grid.resize(size + 2);
	ghost_grid[0] = grid[0] - dr;
	ghost_grid.back() = grid.back() + dr;
	for (size_t i = 1; i < ghost_grid.size() - 1; i++) {
		ghost_grid[i] = grid[i - 1];
	}
	// initializing ghost grids: rho, U, and E
	ghost_size = ghost_grid.size();
	ghost_rho.resize(ghost_size);
	ghost_U.resize(ghost_size);
	ghost_E.resize(ghost_size);
	ghost_p.resize(ghost_size);
	ghost_H.resize(ghost_size);

	std::fill(ghost_rho.begin(), ghost_rho.end(), 1.0);
	std::fill(ghost_U.begin(), ghost_U.end(), 0.0);
	std::fill(ghost_E.begin(), ghost_E.end(), 0.0);

	// initialize the state, flux, and source vector variables
	for (size_t iW = 0; iW < W.size(); iW++) {
		W[iW] = pde_state(ghost_size, 0.0);
		F[iW] = pde_state(ghost_size, 0.0);
		S[iW] = pde_state(ghost_size, 0.0);
		Qj[iW] = pde_state(size, 0.0);
		dissipation[iW] = pde_state(size, 0.0);
		residual[iW] = pde_state(size, 0.0);
	}

	// array of grid ^ order
	grid_to_order.resize(size);
	std::transform(grid.begin(), grid.end(), grid_to_order.begin(), [m_order](double x) { return std::pow(x, m_order); });

	// defining boundary conditions
	std::array<size_t, 3> left_ids = { 0, 1, 2 };
	std::array<size_t, 3> right_ids = { size - 1, size - 2, size - 3 };

	rhoLowerBC = agnosticBCs(left_ids, validBCs::GRADIENT);
	rhoUpperBC = agnosticBCs(right_ids, validBCs::GRADIENT);
	uLowerBC = agnosticBCs(left_ids, validBCs::REFLECTIVE);
	uUpperBC = agnosticBCs(right_ids, validBCs::GRADIENT);
	ELowerBC = agnosticBCs(left_ids, validBCs::GRADIENT);
	EUpperBC = agnosticBCs(right_ids, validBCs::GRADIENT);
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

void EulerSol::operator()(const pde_state& W, pde_state& dWdt, const double t)
{
	// Defining the system of equations
	for (size_t i = 0; i < size; i++) {
		rho[i] = W[i] / grid_to_order[i];
		rho_U[i] = W[i + size] / grid_to_order[i];
		rho_E[i] = W[i + 2 * size] / grid_to_order[i];

		// getting primatives u and E by dividing rho and U
		u[i] = rho_U[i] / rho[i];
		E[i] = rho_E[i] / rho[i];

		// fill ghost cells starting with index 1
		ghost_rho[i + 1] = rho[i];
		ghost_E[i + 1] = E[i];
		ghost_E[i + 1] = u[i];
	}

	// apply boundary conditions
	rhoLowerBC(ghost_rho, ghost_grid);
	rhoUpperBC(ghost_rho, ghost_grid);
	uLowerBC(ghost_U, ghost_grid);
	uUpperBC(ghost_U, ghost_grid);
	ELowerBC(ghost_E, ghost_grid);
	EUpperBC(ghost_E, ghost_grid);

	// Apply Equations of State on Ghost Grid
	for (size_t iG = 0; iG < ghost_size; iG++)
	{
		ghost_p[iG] = ghost_rho[iG] * (gamma - 1.0) * (ghost_E[iG] - 0.5 * pow(ghost_U[iG], 2));
		ghost_H[iG] = ghost_p[iG] + ghost_p[iG] / ghost_rho[iG];
		ghost_E[iG] = sqrt(gamma * ghost_p[iG] / ghost_rho[iG]);

		// develop W - state vector variable
		this->W[0][iG] = ghost_rho[iG] * pow(ghost_grid[iG], order);
		this->W[1][iG] = ghost_rho[iG] * ghost_U[iG] * pow(ghost_grid[iG], order);
		this->W[2][iG] = ghost_rho[iG] * ghost_E[iG] * pow(ghost_grid[iG], order);

		// develop F - flux vector variable
		this->F[0][iG] = ghost_rho[iG] * ghost_U[iG] * pow(ghost_grid[iG], order);
		this->F[1][iG] = (ghost_rho[iG] * pow(ghost_U[iG], 2) + ghost_p[iG]) * pow(ghost_grid[iG], order);
		this->F[2][iG] = ghost_rho[iG] * ghost_U[iG] * ghost_H[iG] * pow(ghost_grid[iG], order);

		// develop S - source term variable
		this->S[1][iG] = order * ghost_p[iG] * pow(ghost_grid[iG], order);
	}

	// calculate second-order Euler flux and dissipation flux
	JST_2ndOrderEulerFlux(F, Qj);
	JST_DissipationFlux(dissipation, this->W, ghost_p, ghost_U, ghost_cs, alpha, beta);

	// calculating residuals
	for (size_t iflux = 0; iflux < residual.size(); iflux++) {
		for (size_t i = 0; i < residual[i].size(); i++) {
			residual[iflux][i] = S[iflux][i] - 1.0 / dr * (Qj[iflux][i] - dissipation[iflux][i]);
		}
	}

	// loading residual into dWdt
	for (size_t i = 0; i < size; i++) {
		dWdt[i] = residual[0][i];
		dWdt[i + size] = residual[1][i];
		dWdt[i + 2 * size] = residual[2][i];
	}
}


/*********************** SEDOV BLAST SOLUTION CLASS ***********************/
SedovBlast::SedovBlast(double ScaleLen__m,	// length scale
	double DomainLen__m,// size of the domain
	double RExpl__m,	// radius of explosion
	double PExpl__Pa,	// pressure of explosion
	double tFinal__s,	// final simulation time
	double rho0__kgpm3, // ambient air density, kg / m ^ 3
	double P0__Pa,		// ambient air pressure, Pa
	size_t order,		// order of the equations, 0 = cartesian, 1 - cylindrical, 2 = spherical
	double gamma,		// ratio of specific heats, N / A
	size_t minNGridPts	// minimum number of grid points
)
{
	//
}

void SedovBlast::solve()
{
	//ODEs = std::make_unique<EulerSol>()
}