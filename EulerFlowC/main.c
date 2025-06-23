/*
Author: Hugh Morgan
Description: solve the Euler Equations in C.
Date: 2025-06-17
*/
#include <stdio.h>
#include <math.h>
#include "pdes.h"

int main() {
    printf("Hello World!\n");
    size_t n_grid_pts = 40;
    double alpha = 1.0;
    double grid_size = 1.0;
    double t_final = 10.0;
    double out_every = 0.2;
    double y0[n_grid_pts];

    // define initial conditions
    for (size_t i = 0; i < n_grid_pts; ++i) {
        if (i > (n_grid_pts/3) && i < (2*n_grid_pts/3)) {
            y0[i] = 1;
        }
        else y0[i] = 0;
    }
    // defining boundary conditions
    BC1D bc;
    bc.left_bc  = &fixed_gradient_bc;
    bc.right_bc = &fixed_gradient_bc;
    bc.left_bc_val  = 0.0;
    bc.right_bc_val = 0.0;

    int res = heat_eqn(n_grid_pts, alpha, grid_size, t_final, out_every, &bc, y0);
    return 0;
}