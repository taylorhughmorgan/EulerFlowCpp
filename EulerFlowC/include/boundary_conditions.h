/* Boundary Conditions generic to an 1D PDE */
#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include <stdio.h>

// declare boundary condition typedef
typedef enum {LEFT, RIGHT} SIDE;
typedef void (*BoundaryCondition)(double*, const double, const double, const size_t);

typedef struct {
    BoundaryCondition left_bc;  // left boundary condition
    BoundaryCondition right_bc; // right boundary condition
    double left_bc_val;         // value for left boundary condition
    double right_bc_val;        // value for right boundary condition
} BC1D;

// define fixed value or Dirclet boundary condition
void fixed_value_bc(double* arr, const double val, const double dx, const size_t idx);

// fixed gradient or Nuemann boundary condition
void fixed_gradient_bc(double* arr, const double val, const double dx, const size_t idx);

// extrapolated boundary condition
void extrapolated_bc(double* arr, const double val, const double dx, const size_t idx);
#endif