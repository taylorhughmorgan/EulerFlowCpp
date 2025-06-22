/* Boundary Conditions generic to an 1D PDE */
#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include <stdio.h>

// declare boundary condition typedef
typedef enum {LEFT, RIGHT} SIDE;
typedef void (*BoundaryCondition)(double*, const double, const double, const size_t);

typedef struct {
    BoundaryCondition left_bc;
    BoundaryCondition right_bc;
    double left_bc_val;
    double right_bc_val;
} BC1D;

// define fixed value or Dirclet boundary condition
void fixed_value_bc(double* arr, const double val, const double dx, const size_t idx) {
    (void)dx;
    arr[idx] = val;
}

// fixed gradient or Nuemann boundary condition
void fixed_gradient_bc(double* arr, const double val, const double dx, const size_t idx) {
    if (idx == 0) {
        arr[0] = arr[1] - dx * val;
    }
    else {
        arr[idx] = arr[idx - 1] + dx * val; 
    }
}

// extrapolated boundary condition
void extrapolated_bc(double* arr, const double val, const double dx, const size_t idx) {
    if (idx == 0) {
        arr[0] = -arr[2] + 2 * arr[1] + dx * val;
    }
    else {
        arr[idx] = -arr[idx - 2] + 2 * arr[idx - 1] + dx * val;
    }
}

#endif