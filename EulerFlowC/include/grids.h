/*
Define 1D grids
*/
#ifndef GRIDS_H
#define GRIDS_H

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>

typedef struct {
    size_t size;    /* number of grid points */
    double* data;   /* contiguous data buffer, length == size */
} Grid1D;

// constructor and destructor for Grid1D object
Grid1D * Grid1D_alloc(const size_t size);

void Grid1D_free(Grid1D * grid);
#endif