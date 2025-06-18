/*
Author: Hugh Morgan
Description: solve the Euler Equations in C.
Date: 2025-06-17
*/
#include <math.h>
#include "include/grids.h"

#define NCELLS 100

typedef struct {
    double density;
    double momentum;
    double energy;
} State;

typedef struct {
    State (*compute_flux)(const State* left, const State* right);
    const char* name;
} FluxScheme;

State lax_friedrich_flux(const State* left, const State* right) {
    // Placeholder implementation of 
    State flux;
    flux.density = 0.5 * (left->density + right->density);
    flux.momentum = 0.5 * (left->momentum + right->momentum);
    flux.energy = 0.5 * (left->energy + right->energy);
    return flux;
}

State roe_flux(const State* left, const State* right) {
    // Placeholder implementation of Roe flux - real version would be more complex
    State flux;
    flux.density = 0.5 * (left->density + right->density);
    flux.momentum = 0.5 * (left->momentum + right->momentum);
    flux.energy = 0.5 * (left->energy + right->energy);
    return flux;
}

void compute_fluxes(const FluxScheme* scheme, const State* cells, size_t n_cells, State* fluxes_out) {
    // compute fluxes
    for (size_t i = 0; i < n_cells - 1; ++i) {
        fluxes_out[i] = scheme->compute_flux(&cells[i], &cells[i + 1]);
    }
}

int main() {
    printf("Hello World!\n");
    FluxScheme lax = { .compute_flux = lax_friedrich_flux, .name = "Lax-Friedrichs" };
    FluxScheme roe = { .compute_flux = roe_flux, .name = "Roe" };

    // define a 1D grid
    Grid1D * grid = Grid1D_alloc(NCELLS);
    State cells[NCELLS];
    State fluxes[NCELLS - 1];

    // set up the grid and cell values
    double x_max = 2 * M_PI;
    double dx = x_max / (double)(NCELLS - 1);
    for (size_t i = 0; i < NCELLS; ++i) {
        grid->data[i] = dx * (double)i;
        cells[i].density = sin(grid->data[i]);
        cells[i].momentum = cos(grid->data[i]);
        cells[i].energy = -sin(grid->data[i]);
    }

    // compute lax-friedrich flux
    compute_fluxes(&lax, cells, NCELLS, fluxes);
    // compute roe fluxes
    compute_fluxes(&roe, cells, NCELLS, fluxes);

    printf("Computing Fluxes Complete!\n");
    // freeing memory
    Grid1D_free(grid);
    
    return 0;
}