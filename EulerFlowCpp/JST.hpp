#pragma once
#include "BoundaryConditions.hpp"

// Fluxes
enum class Fluxes { F_j, Q_j, D_j };
typedef std::array<pde_state, 3> flux_state;

/* Numerical, finite-differencing methods used by JST (Jameson-Shmidt-Turkel */
void d3dx3_fwd(pde_state& y, pde_state& d3ydx3);

void d3dx3_bckwrd(pde_state& y, pde_state& d3ydx3);

void JST_2ndOrderEulerFlux(flux_state& fluxVector, flux_state& dissipationFlux);

void JST_DissipationFlux(
	flux_state& D_j, // dissipation flux
	flux_state& W,	 // state vector variable
	pde_state& p,					 // pressure
	pde_state& u,					 // velocity
	pde_state& cs,					 // speed of sound
	std::array<double, 2>& alpha,	 // dissipative flux terms for spatial differencing
	std::array<double, 2>& beta		 // dissipative flux terms for spatial differencing
);