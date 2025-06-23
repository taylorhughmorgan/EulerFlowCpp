#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "pdes.h"

// structure containing simulation parameters
typedef struct {
    size_t size;        // size of the problem
    double alpha;       // diffusion coefficient
    double dx;          // grid size
    double dt;          // time step
} PARAMS;


// function for applying BCs
void apply_BC(BC1D * bc, double* y, PARAMS* sim_params) {
    bc->left_bc(y, 0.0, sim_params->dx, 0);
    bc->right_bc(y, 0.0, sim_params->dx, sim_params->size - 1);    
}

int func(double t, const double y[], double dydt[], void *params)
{
    /* RHS of the heat equation:
        \partial U / \partial t =  */
    (void)(t); /* avoid unused parameter warning */
    PARAMS sim_params = *(PARAMS *)params;
    double coef = sim_params.alpha * sim_params.dt / (sim_params.dx * sim_params.dx);

    for (size_t i = 0; i < sim_params.size - 1; ++i) {
        dydt[i] = coef * (y[i - 1] - 2 * y[i] + y[i + 1]);
    }
    return GSL_SUCCESS;
}

void save_result(FILE * fptr, const double t, const double* y_sol, size_t size) {
    // save the results to a file
    fprintf(fptr, "%.5e, ", t);
    for (size_t i = 0; i < size; ++i) {
        fprintf(fptr, "%.5e, ", y_sol[i]);
    }
    fprintf(fptr, "\n");
}

int heat_eqn(size_t n_grid_pts, // number of grid pts 
    double alpha,               // heat diffusivity
    double grid_size,           // grid size
    double tfinal,              // final time
    double out_every,           // save time 
    BC1D* bc,                   // boundary conditions
    double* y                   // initial condition
)
{
    // intializing the thermodynamic parameters and grid size
    PARAMS sim_params;
    sim_params.alpha = alpha;
    sim_params.size = n_grid_pts;

    // setting up the grid
    double xmin = 0.0, xmax = grid_size;
    double t = 0.0;
    sim_params.dx = (xmax - xmin) / (double)sim_params.size;
    double courant = 0.5;
    double dtmax = courant * sim_params.dx * sim_params.dx / sim_params.alpha;
    sim_params.dt = dtmax;

    // define the ODE solver
    printf("Setting up the system of equations\n");
    gsl_odeiv2_system sys = {func, NULL, sim_params.size, &sim_params};

    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

    // define file output
    char filename[] = "heat_eqn_sol.csv";
    FILE * fptr = fopen(filename, "w");
    if (fptr == NULL) {
        printf("Error opening the file %s", filename);
        return -1;
    }
    else {
        printf("Setting up the output file\n");
        fprintf(fptr, "# Solution to the heat equation for alpha=%f on grid xmin=%f, xmax=%f\n", sim_params.alpha, xmin, xmax);
        fprintf(fptr, "time, ");
        for (size_t i = 0; i < sim_params.size; ++i) {
            fprintf(fptr, "y[%zu], ", i);
        }
        fprintf(fptr, "\n");
    }
    // save the initial state
    save_result(fptr, t, y, sim_params.size);
    
    // stepping through each time
    double output_tracker = 0.0;
    printf("Begin looping through solution\n");
    while (t < tfinal) {
        output_tracker += dtmax;
        double t_f = t + dtmax;
        int status = gsl_odeiv2_driver_apply(d, &t, t_f, y);
        
        if (status != GSL_SUCCESS) {
            printf("error, return value=%d\n", status);
            break;
        }
        // apply boundary conditions
        apply_BC(bc, y, &sim_params);
        
        // saving the result
        if (output_tracker >= out_every) {
            save_result(fptr, t, y, sim_params.size);
            printf("Solution reached and result saved for t=%f\n", t);
            output_tracker = 0.0;
        }
    }

    printf("Solution reached!\nSolution saved to %s", filename);
    fclose(fptr); // close the file output
    gsl_odeiv2_driver_free (d); // free the memory
    return GSL_SUCCESS;
}