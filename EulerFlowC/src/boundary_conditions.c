#include "boundary_conditions.h"

// define fixed value or Dirclet boundary condition
void fixed_value_bc(double* arr, const double val, const double dx, const size_t idx) {
    (void)dx; // avoid unused parameter warning
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