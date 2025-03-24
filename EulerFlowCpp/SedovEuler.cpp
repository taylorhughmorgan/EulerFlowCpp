#include "SedovEuler.hpp"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;

/**************** OBSERVER CLASS *************************************/
struct push_back_state_and_time
{
	std::vector< pde_state >& m_states;
	std::vector< double >& m_times;

	push_back_state_and_time(std::vector< pde_state >& states, std::vector< double >& times)
		: m_states(states), m_times(times) {
	}

	void operator()(const pde_state& x, double t)
	{
		m_states.push_back(x);
		m_times.push_back(t);
	}
};

/*************** EULER SOLUTION CLASS ********************************/
EulerSol::EulerSol(
	pde_state& m_grid,
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

	rhoLowerBC = agnosticBCs(left_ids, m_BCs[Fields::RHO][0]);
	rhoUpperBC = agnosticBCs(right_ids, m_BCs[Fields::RHO][1]);
	uLowerBC = agnosticBCs(left_ids, m_BCs[Fields::VEL][0]);
	uUpperBC = agnosticBCs(right_ids, m_BCs[Fields::VEL][1]);
	ELowerBC = agnosticBCs(left_ids, m_BCs[Fields::ENERGY][0]);
	EUpperBC = agnosticBCs(right_ids, m_BCs[Fields::ENERGY][1]);
}

void EulerSol::createICs(const pde_state& rho0, const pde_state& v0, const pde_state& p0, pde_state& W0)
{
	// Convert primative variables into the initial conditions - W0
	size_t out_size = rho0.size() + v0.size() + p0.size();
	W0.resize(out_size);
	// using equations of state, calculate internal energy
	pde_state E0(rho0.size());

	for (size_t i = 0; i < rho0.size(); i++) {
		E0[i] = p0[i] / (rho0[i] * (gamma - 1.0)) + 0.5 * pow(v0[i], 2);
		W0[i] = rho0[i] * grid_to_order[i];
		W0[i + size] = rho0[i] * v0[i] * grid_to_order[i];
		W0[i + size] = rho0[i] * E0[i] * grid_to_order[i];
	}
}

void EulerSol::conv2Primatives(const pde_state& W_out, pde_state& rho_out, pde_state& u_out, pde_state& E_out, pde_state& p_out)
{
	// Convert W result to primative values -> rho, u, E, and P
	// size the arrays correctly
	rho_out.resize(this->size);
	u_out.resize(this->size);
	E_out.resize(this->size);
	p_out.resize(this->size);
	for (size_t i = 0; i < this->size; i++) {
		rho_out[i] = W_out[i] / grid_to_order[i];
		double rho_U = W_out[i + size] / grid_to_order[i];
		double rho_E = W_out[i + 2 * size] / grid_to_order[i];
		
		u_out[i] = rho_U / rho_out[i];
		E_out[i] = rho_E / rho_out[i];
		p_out[i] = rho_out[i] * (gamma - 1) * (E_out[i] - 0.5 * pow(u_out[i], 2));
	}
}

void EulerSol::operator()(const pde_state& x, pde_state& dxdt, const double t)
{
	// Defining the system of equations
	for (size_t i = 0; i < size; i++) {
		rho[i] = x[i] / grid_to_order[i];
		rho_U[i] = x[i + size] / grid_to_order[i];
		rho_E[i] = x[i + 2 * size] / grid_to_order[i];

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
		W[0][iG] = ghost_rho[iG] * pow(ghost_grid[iG], order);
		W[1][iG] = ghost_rho[iG] * ghost_U[iG] * pow(ghost_grid[iG], order);
		W[2][iG] = ghost_rho[iG] * ghost_E[iG] * pow(ghost_grid[iG], order);

		// develop F - flux vector variable
		F[0][iG] = ghost_rho[iG] * ghost_U[iG] * pow(ghost_grid[iG], order);
		F[1][iG] = (ghost_rho[iG] * pow(ghost_U[iG], 2) + ghost_p[iG]) * pow(ghost_grid[iG], order);
		F[2][iG] = ghost_rho[iG] * ghost_U[iG] * ghost_H[iG] * pow(ghost_grid[iG], order);

		// develop S - source term variable
		S[1][iG] = order * ghost_p[iG] * pow(ghost_grid[iG], order);
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
		dxdt[i] = residual[0][i];
		dxdt[i + size] = residual[1][i];
		dxdt[i + 2 * size] = residual[2][i];
	}
}


/*********************** SEDOV BLAST SOLUTION CLASS ***********************/
SedovBlast::SedovBlast(
	double m_ScaleLen__m,		// length scale
	double m_DomainLen__m,		// size of the domain
	double m_RExpl__m,			// radius of explosion
	double m_PExpl__Pa,			// pressure of explosion
	double m_tFinal__s,			// final simulation time
	double m_rho0__kgpm3, // ambient air density, kg / m ^ 3
	double m_P0__Pa,		// ambient air pressure, Pa
	size_t m_order,			// order of the equations, 0 = cartesian, 1 - cylindrical, 2 = spherical
	double m_gamma,			// ratio of specific heats, N / A
	size_t m_minNGridPts	// minimum number of grid points
)
{
	// Convert the parameters of the Sedov Blast to nondimensional form, for speed and numerical stability.
	this->ScaleLen__m = m_ScaleLen__m;
	this->DomainLen__m = m_DomainLen__m;
	this->RExpl__m = m_RExpl__m;
	this->PExpl__Pa = m_PExpl__Pa;
	this->tFinal__s = m_tFinal__s;
	this->rho0__kgpm3 = m_rho0__kgpm3;
	this->P0__Pa = m_P0__Pa;
	this->order = m_order;
	this->gamma = m_gamma;
	this->minNGridPts = m_minNGridPts;

	// Calculate the internal energy of the explosion
	size_t n = order + 1;
	EExpl__J = pow(pi, double(n) / 2.0) * boost::math::tgamma(double(n) / 2.0 + 1.0) * pow(RExpl__m, n);

	// dimensionless parameters: scale using rho0, P0, and diameter
	double UScale = sqrt(P0__Pa / rho0__kgpm3);
	double lenStar = DomainLen__m / ScaleLen__m;
	double rExplStar = RExpl__m / ScaleLen__m;
	double pExplStar = PExpl__Pa / P0__Pa;
	double tFinStar = tFinal__s * UScale / ScaleLen__m;

	// set up the radial grid, we want at least 10 points for the explosion
	nGridPts = size_t(std::ceil(lenStar / rExplStar * 10.0));
	nGridPts = std::max(nGridPts, minNGridPts);

	double rMinStar = std::min(rExplStar / 10.0, lenStar / 100.0);
	linspace<double>(rMinStar, lenStar, nGridPts, grid);
	linspace<double>(0.0, tFinStar, minNGridPts, times);

	dr = grid[1] - grid[0];
	dt = times[1] - times[0];

	// setting the initial conditions - setting rho0Star and p0Star to 1.0 for dimensionless variables
	rho0Star.resize(nGridPts);
	p0Star.resize(nGridPts);
	v0Star.resize(nGridPts);
	std::fill(rho0Star.begin(), rho0Star.end(), 1.0);
	std::fill(p0Star.begin(), p0Star.end(), 1.0);
	std::fill(v0Star.begin(), v0Star.end(), 0.0);

	for (size_t i = 0; i < nGridPts; i++) {
		if (grid[i] < rExplStar) {
			p0Star[i] = pExplStar;
		}
		else break;
	}

	// time and grid in dimensional/metric scale
	for (size_t i = 0; i < nGridPts; i++) {
		rGrid__m.push_back(grid[i] * ScaleLen__m);
		tGrid__s.push_back(times[i] * ScaleLen__m / UScale);
	}
}


void SedovBlast::solve()
{
	//
	std::array<double, 2> alpha = { 0.5, 0.5 };
	std::array<double, 2> beta = { 0.25, 0.5 };
	bc_map boundary_conditions;

	// system of equations
	EulerSol ODEs(grid, boundary_conditions, order, alpha, beta, gamma);
	// intial conditions
	pde_state W0Star;
	ODEs.createICs(rho0Star, p0Star, v0Star, W0Star);

	std::cout << "Solving the Euler Equation as a system of ODES. \n" <<
		"t_range = [" << 0.0 << ", " << times.back() << "](dimensionless) \n" <<
		"nGridPts = " << nGridPts << std::endl <<
		"r_range = [" << grid[0] << "," << grid.back() << "](dimensionless)";

	// define observer
	std::vector<pde_state> states_sol;
	pde_state times_sol;
	push_back_state_and_time observer(states_sol, times_sol);

	// define numerical stepper
	auto stepper = runge_kutta_cash_karp54<pde_state>();

	integrate_times(stepper, ODEs, W0Star, times.begin(), times.end(), dt, observer);

	std::cout << "Solution reached!\n";
}