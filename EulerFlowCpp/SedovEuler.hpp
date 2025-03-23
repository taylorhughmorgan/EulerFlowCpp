#pragma once
#include "JST.hpp"
#include <memory>
#include <algorithm>

typedef std::map<std::string, std::vector<std::string>> bc_map;

class EulerSol
{
/*Solve Euler's system of equations describing the behavior of inviscid, compressible flow by reducing to a system of ordinary differential equations. */
private:
	pde_state rho, rho_U, rho_E, u, E;
	BCFunction rhoUpperBC, rhoLowerBC, uUpperBC, uLowerBC, EUpperBC, ELowerBC;
	// paramters on the ghost grid
	pde_state ghost_grid, ghost_rho, ghost_U, ghost_E, ghost_p, ghost_H, ghost_cs;
	flux_state W, F, S, Qj, dissipation, residual;
public:
	pde_state grid, grid_to_order;
	double gamma;
	size_t order, size, ghost_size;
	std::array<double, 2> alpha, beta;
	double dr;

	// default constructor
	EulerSol(pde_state& m_grid,	// computational grid
		bc_map m_BCs,			// 
		size_t m_order = 0,		// order of equations, 0=cartesian, 1=cylindrical/polar, 2=spherical
		std::array<double, 2> m_alpha = { 0.5, 0.5 }, // dissipative flux terms for spatial differencing
		std::array<double, 2> m_beta = { 0.25, 0.5 }, // dissipative flux terms for spatial differencing
		double m_gamma = 1.4	// ratio of specific heats
		);

	// create initial conditions
	void createICs(pde_state& rho0, pde_state& v0, pde_state& p0, pde_state& W0);
	// convert output to primatives
	void conv2Primatives(pde_state& W, pde_state& rho, pde_state& u, pde_state& E, pde_state& p);
	// function call for ODE integrator
	void operator()(const pde_state& W, pde_state& dWdt, const double t);
};

class SedovBlast
{
public:
	// system of equations
	std::unique_ptr<EulerSol> ODEs;

	// default constructor
	SedovBlast(double ScaleLen__m,	// length scale
		double DomainLen__m,		// size of the domain
		double RExpl__m,			// radius of explosion
		double PExpl__Pa,			// pressure of explosion
		double tFinal__s,			// final simulation time
		double rho0__kgpm3 = 1.225, // ambient air density, kg / m ^ 3
		double P0__Pa = 101325,		// ambient air pressure, Pa
		size_t order = 0,			// order of the equations, 0 = cartesian, 1 - cylindrical, 2 = spherical
		double gamma = 1.4,			// ratio of specific heats, N / A
		size_t minNGridPts = 500	// minimum number of grid points
	);
	// solve the system of equations over the time domain
	void solve();
};