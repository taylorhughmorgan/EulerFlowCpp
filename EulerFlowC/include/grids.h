/*
Define 1D grids
*/
#ifndef GRIDS_H
#define GRIDS_H

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    size_t size;    /* number of grid points */
    double* data;   /* contiguous data buffer, length == size */
} Grid1D;


Grid1D * Grid1D_alloc(const size_t size) {
    // allocate grid struct
    Grid1D * grid = (Grid1D *)malloc(sizeof(Grid1D));
    if (grid == NULL) {
        fprintf(stderr, "Unable to allocate Grid1D!\n");
        free(grid);
        return NULL;
    }
    grid->size = size;

    // allocate array size
    grid->data = (double *)malloc(size * sizeof(double));
    if (grid->data == NULL) {
        fprintf(stderr, "Unable to allocate data to Grid1D.\n");
        free(grid->data);
        return NULL;
    }
    return grid;
}

void Grid1D_free(Grid1D * grid) {
    // de-allocate/free memory
    if (!grid) {
        free(grid->data);
        free(grid);
    }
}
#endif