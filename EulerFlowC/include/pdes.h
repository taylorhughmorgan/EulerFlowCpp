#ifndef PDES_H
#define PDES_H
#include "boundary_conditions.h"
#include <gsl/gsl_odeiv2.h>

// coordinate system definition: spherical, polar/cylindrical, and cartesian
typedef enum Coords {Spherical, Polar, Cartesian} CoordSys;

// 1D heat equation 

int heat_eqn(size_t n_grid_pts, // number of grid pts 
    double alpha,               // heat diffusivity
    double grid_size,           // grid size
    double tfinal,              // final time
    double out_every,           // save time 
    BC1D * bc,                  // boundary conditions
    gsl_odeiv2_step_type stepper, // time step integrator
    double * y                  // initial condition
);

#endif