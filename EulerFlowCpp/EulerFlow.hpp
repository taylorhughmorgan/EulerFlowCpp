#pragma once
#include "JST.hpp"
#include "math_utils.hpp"
#include <memory>
#include <algorithm>
#include <map>

// typedef for constructing boundary conditions
enum class Fields { RHO, VEL, ENERGY };
typedef std::map<Fields, std::array<validBCs,2>> bc_map;

class EulerSol
{
/*Solve Euler's system of equations describing the behavior of inviscid, compressible flow by reducing to a system of ordinary differential equations. */
private:
	pde_state rho, rho_U, rho_E, u, E;
	BCFunction rhoUpperBC, rhoLowerBC, uUpperBC, uLowerBC, EUpperBC, ELowerBC;
	// paramters on the ghost grid
	pde_state ghost_grid, ghost_rho, ghost_U, ghost_E, ghost_p, ghost_H, ghost_cs;
	flux_state W, F, S, Qj, Dj, residual;
	JST_DissipFlux DissipFlux;
public:
	pde_state grid, grid_to_order;
	double gamma;
	size_t order, size, ghost_size;
	std::array<double, 2> alpha, beta;
	double dr;

	// default constructor
	EulerSol(
		pde_state& m_grid,		// computational grid
		bc_map m_BCs = { {Fields::RHO, {validBCs::GRADIENT, validBCs::GRADIENT} },
						{Fields::VEL, {validBCs::REFLECTIVE, validBCs::GRADIENT} },
						{Fields::ENERGY, {validBCs::GRADIENT, validBCs::GRADIENT} } }, // 
		size_t m_order = 0,		// order of equations, 0=cartesian, 1=cylindrical/polar, 2=spherical
		std::array<double, 2> m_alpha = { 0.5, 0.5 }, // dissipative flux terms for spatial differencing
		std::array<double, 2> m_beta = { 0.25, 0.5 }, // dissipative flux terms for spatial differencing
		double m_gamma = 1.4	// ratio of specific heats
		);

	// create initial conditions
	void createICs(const pde_state& rho0, const pde_state& v0, const pde_state& p0, pde_state& W0);
	// convert output to primatives
	void conv2Primatives(const pde_state& W_out, pde_state& rho_out, pde_state& u_out, pde_state& E_out, pde_state& p_out);
	// convert all output to primatives
	void convAllPrimatives(const std::vector<pde_state>& W_out, std::vector<pde_state>& rho_out, std::vector<pde_state>& u_out, std::vector<pde_state>& E_out, std::vector<pde_state>& p_out);
	// function call for ODE integrator
	void operator()(const pde_state& W, pde_state& dWdt, const double t);
};
