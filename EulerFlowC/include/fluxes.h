#ifndef FLUXES_H
#define FLUXES_H

#include "grids.h"

/* Flux states */
typedef struct {
    double density;
    double momentum;
    double energy;
} State;

State * State_alloc(const size_t size) {
    // dynamically allocate State
    State * state = (State *)malloc(size * sizeof(State));
    if (state == NULL) {
        fprintf(stderr, "Unable to allocate state.\n");
        return NULL;
    }
    return state;
}

void State_free(State * state) {
    // de-allocate/free state
    if (!state) {
        free(state);
    }
}

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
#endif