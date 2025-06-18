/*
Author: Hugh Morgan
Description: solve the Euler Equations in C.
Date: 2025-06-17
*/
#include <math.h>
#include "include/fluxes.h"


int main() {
    printf("Hello World!\n");
    FluxScheme lax = { .compute_flux = lax_friedrich_flux, .name = "Lax-Friedrichs" };
    FluxScheme roe = { .compute_flux = roe_flux, .name = "Roe" };

    // define a 1D grid
    const size_t n_cells = 100;
    Grid1D * grid = Grid1D_alloc(n_cells);
    State * cells = State_alloc(n_cells);
    State * fluxes = State_alloc(n_cells - 1);

    // set up the grid and cell values
    double x_max = 2 * M_PI;
    double dx = x_max / (double)(n_cells - 1);
    for (size_t i = 0; i < n_cells; ++i) {
        grid->data[i] = dx * (double)i;
        cells[i].density = sin(grid->data[i]);
        cells[i].momentum = cos(grid->data[i]);
        cells[i].energy = -sin(grid->data[i]);
    }

    // compute lax-friedrich flux
    compute_fluxes(&lax, cells, n_cells, fluxes);
    // compute roe fluxes
    compute_fluxes(&roe, cells, n_cells, fluxes);

    printf("Computing Fluxes Complete!\n");
    
    // freeing memory
    Grid1D_free(grid);
    State_free(cells);
    State_free(fluxes);

    return 0;
}