#pragma once
#include "BoundaryConditions.hpp"

/* Numerical, finite-differencing methods used by JST (Jameson-Shmidt-Turkel */
void d3dx3_fwd(pde_state& y, pde_state& d3ydx3);

void d3dx3_bckwrd(pde_state& y, pde_state& d3ydx3);