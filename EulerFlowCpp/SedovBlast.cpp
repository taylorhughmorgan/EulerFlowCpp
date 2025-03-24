#include "SedovBlast.hpp"
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