#pragma once
#include "JST.hpp"
#include "math_utils.hpp"
#include <memory>
#include <algorithm>

typedef std::map<std::string, validBCs> bc_map;

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
	EulerSol(
		pde_state& m_grid,		// computational grid
		bc_map m_BCs,			// 
		size_t m_order = 0,		// order of equations, 0=cartesian, 1=cylindrical/polar, 2=spherical
		std::array<double, 2> m_alpha = { 0.5, 0.5 }, // dissipative flux terms for spatial differencing
		std::array<double, 2> m_beta = { 0.25, 0.5 }, // dissipative flux terms for spatial differencing
		double m_gamma = 1.4	// ratio of specific heats
		);

	// create initial conditions
	void createICs(const pde_state& rho0, const pde_state& v0, const pde_state& p0, pde_state& W0);
	// convert output to primatives
	void conv2Primatives(const pde_state& W_out, pde_state& rho_out, pde_state& u_out, pde_state& E_out, pde_state& p_out);
	// function call for ODE integrator
	void operator()(const pde_state& W, pde_state& dWdt, const double t);
};

class SedovBlast
{
private:
	double dr;
public:
	double ScaleLen__m, DomainLen__m, RExpl__m, PExpl__Pa, tFinal__s, rho0__kgpm3, P0__Pa, gamma, EExpl__J;
	size_t order, minNGridPts, nGridPts;
	pde_state grid, times, rGrid__m, tGrid__s;
	pde_state rho0Star, p0Star, v0Star;

	// default constructor
	SedovBlast(
		double m_ScaleLen__m,		// length scale
		double m_DomainLen__m,		// size of the domain
		double m_RExpl__m,			// radius of explosion
		double m_PExpl__Pa,			// pressure of explosion
		double m_tFinal__s,			// final simulation time
		double m_rho0__kgpm3 = 1.225, // ambient air density, kg / m ^ 3
		double m_P0__Pa = 101325,		// ambient air pressure, Pa
		size_t m_order = 0,			// order of the equations, 0 = cartesian, 1 - cylindrical, 2 = spherical
		double m_gamma = 1.4,			// ratio of specific heats, N / A
		size_t m_minNGridPts = 500	// minimum number of grid points
	);
	// solve the system of equations over the time domain
	void solve();
};