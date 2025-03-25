#pragma once
#include "EulerFlow.hpp"
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

class SedovBlast
{
private:
	double dr, dt;
	double UScale, lenStar;
public:
	double ScaleLen__m, DomainLen__m, RExpl__m, PExpl__Pa, tFinal__s, rho0__kgpm3, P0__Pa, gamma, EExpl__J;
	size_t order, minNGridPts, nGridPts;
	pde_state grid, times, rGrid__m, tGrid__s;
	pde_state rho0Star, p0Star, v0Star;
	// solution variables
	std::vector<pde_state> rho_sol, p_sol, E_sol, u_sol;
	// system of equations
	std::unique_ptr<EulerSol> ODEs;


	// default constructor
	SedovBlast(
		double m_ScaleLen__m,			// length scale
		double m_DomainLen__m,			// size of the domain
		double m_RExpl__m,				// radius of explosion
		double m_PExpl__Pa,				// pressure of explosion
		double m_tFinal__s,				// final simulation time
		double m_rho0__kgpm3 = 1.225,	// ambient air density, kg / m ^ 3
		double m_P0__Pa = 101325,		// ambient air pressure, Pa
		size_t m_order = 0,				// order of the equations, 0 = cartesian, 1 - cylindrical, 2 = spherical
		double m_gamma = 1.4,			// ratio of specific heats, N / A
		size_t m_minNGridPts = 500		// minimum number of grid points
	);
	// solve the system of equations over the time domain
	void solve();
	// save the results to the file
	void save(std::string foutname, json inputs);
};