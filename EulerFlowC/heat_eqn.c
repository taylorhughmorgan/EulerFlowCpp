#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

// structure containing thermodynamic parameters
typedef struct {
    size_t size;        // size of the problem
    double alpha;       // diffusion coefficient
    double dx;          // grid size
    double gamma;       // ratio of specific heats
} Thermo;


int func(double t, const double y[], double dydt[], void *params)
{
    /* RHS of the heat equation:
        \partial U / \partial t =  */
    (void)(t); /* avoid unused parameter warning */
    Thermo thermo = *(Thermo *)params;

    for (size_t i = 0; i < thermo.size - 1; ++i) {
        dydt[i] = thermo.alpha * (y[i - 1] - y[i] + y[i + 1]) / thermo.dx;
    }
    return GSL_SUCCESS;
}

int main () 
{
    Thermo thermo;
    thermo.gamma = 1.4;
    thermo.alpha = 1.0;
    thermo.size = 20;
    // setting up the grid
    double xmin = 0.0, xmax = 1.0;
    double t = 0.0, tfinal = 0.1;
    thermo.dx = (xmax - xmin) / (double)thermo.size;
    double courant = 0.25;
    double dtmax = courant * thermo.dx * thermo.dx / thermo.alpha;

    // define initial conditions
    double y[thermo.size];
    for (size_t i = 0; i < thermo.size; ++i) {
        if (i > 7 && i < 13) {
            y[i] = 1;
        }
        else y[i] = 0;
    }

    // define the ODE solver
    printf("Setting up the system of equations\n");
    gsl_odeiv2_system sys = {func, NULL, thermo.size, &thermo};

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
        fprintf(fptr, "# Solution to the heat equation for alpha=%f on grid xmin=%f, xmax=%f\n", thermo.alpha, xmin, xmax);
        fprintf(fptr, "time, ");
        for (size_t i = 0; i < thermo.size; ++i) {
            fprintf(fptr, "y[%zu], ", i);
        }
        fprintf(fptr, "\n");
    }

    // stepping through each time
    printf("Begin looping through solution\n");
    while (t < tfinal) {
        double t_f = t + dtmax;
        int status = gsl_odeiv2_driver_apply(d, &t, t_f, y);
        
        if (status != GSL_SUCCESS) {
            printf("error, return value=%d\n", status);
            break;
        }
        // saving the result
        fprintf(fptr, "%.5e, ", t);
        for (size_t i = 0; i < thermo.size; ++i) {
            fprintf(fptr, "%.5e, ", y[i]);
        }
        fprintf(fptr, "\n");
        printf("Solution reached for t=%f\n", t);
    }

    printf("Solution reached!\n");
    fclose(fptr); // close the file output
    gsl_odeiv2_driver_free (d); // free the memory
    return 0;
}