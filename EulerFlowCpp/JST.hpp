#pragma once
#include "BoundaryConditions.hpp"
#include "map"

// Fluxes
enum class Fluxes { F_j, Q_j, D_j };

/* Numerical, finite-differencing methods used by JST (Jameson-Shmidt-Turkel */
void d3dx3_fwd(pde_state& y, pde_state& d3ydx3);

void d3dx3_bckwrd(pde_state& y, pde_state& d3ydx3);

void JST_2ndOrderEulerFlux(std::map<Fluxes, pde_state>& fluxVector, std::map<Fluxes, pde_state>& dissipationFlux);

void JSD_DissipationFlux(
	std::map<Fluxes, pde_state>& D_j, // dissipation flux
	std::map<Fluxes, pde_state>& W,	 // state vector variable
	pde_state& p,					 // pressure
	pde_state& u,					 // velocity
	pde_state& cs,					 // speed of sound
	std::array<double, 2> alpha,	 // dissipative flux terms for spatial differencing
	std::array<double, 2> beta		 // dissipative flux terms for spatial differencing
);