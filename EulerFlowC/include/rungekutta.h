/* structs for general-purpose runge-kutta integrator */
#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include <stdlib.h>
#include <stddef.h>
#include <math.h>

// define ODE function pointer
typedef void(*ODEfunction)(double t, const double* y, double* dydt, size_t dim);

// define integrator
typedef struct Integrator {
    void (*step)(struct Integrator*, ODEfunction, double* y, double t, double dt, size_t dim);
    const char* name;
} Integrator;

void rk2_step(Integrator* rk, ODEfunction f, double* y, double t, double dt, size_t dim) {
    // Second-order Runge-Kutta integrator (Heun's method)
    // allocate memory dynamically
    double* k1 = (double*)malloc(dim * sizeof(double));
    double* k2 = (double*)malloc(dim * sizeof(double));
    double* y_temp = (double*)malloc(dim * sizeof(double));

    if (!k1 || !k2 || !y_temp) {
        // handle allocation failure
        free(k1); free(k2); free(y_temp);
        return;
    }

    // first step
    f(t, y, k1, dim);
    for (size_t i = 0; i < dim; ++i) {
        y_temp[i] = y[i] + dt * k1[i];
    }

    // second step
    f(t + dt, y_temp, k2, dim);
    for (size_t i = 0; i < dim; ++i) {
        y[i] += 0.5 * dt * (k1[i] + k2[i]);
    }
    // free memory
    free(k1);
    free(k2);
    free(y_temp);
}

void rk4_step(Integrator* rk, ODEfunction f, double* y, double t, double dt, size_t dim) {
    // Fourth-order Runge-Kutta method
    // dynamically allocate memory
    double* k1 = (double*)malloc(dim * sizeof(double));
    double* k2 = (double*)malloc(dim * sizeof(double));
    double* k3 = (double*)malloc(dim * sizeof(double));
    double* k4 = (double*)malloc(dim * sizeof(double));
    double* y_temp = (double*)malloc(dim * sizeof(double));
    
    if (!k1 || !k2 || !k3 || !k4 || !y_temp) {
        // handle allocation failure
        free(k1); free(k2); free(k3); free(k4); free(y_temp);
        return;
    }

    // first step
    f(t, y, k1, dim);
    for (size_t i = 0; i < dim; ++i) {
        y_temp[i] = y[i] + 0.5 * dt * k1[i];
    }
    // second step
    f(t + 0.5 * dt, y_temp, k2, dim);
    for (size_t i = 0; i < dim; ++i) {
        y_temp[i] = y[i] + 0.5 * dt * k2[i];
    }
    // third step
    f(t + 0.5 * dt, y_temp, k3, dim);
    for (size_t i = 0; i < dim; ++i) {
        y_temp[i] = y[i] + dt * k3[i];
    }
    // fourth step
    f(t + 0.5 * dt, y_temp, k4, dim);

    // weighted solution
    for (size_t i = 0; i < dim; ++i) {
        y[i] += dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }

    // freeing memory
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(y_temp);
}

void integrate(Integrator* rk, ODEfunction f, double* y, double t0, double tf, double dt, size_t dim) {
    // function to integrate from t0 to t_f
    size_t nsteps = (size_t)ceil((tf - t0) / dt);
    double t = t0;
    for (size_t it = 0; it < nsteps; ++it) {
        rk->step(rk, f, y, t, dt, dim);
        t += dt;
    }
}
#endif