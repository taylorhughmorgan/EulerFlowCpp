#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

// structure containing thermodynamic parameters
typedef struct {
    double gamma;       // ratio of specific heats
    double mu_psis;     // kinematic viscosity, psi-s
    double rho0_kgpm3;  // ambient density kg/m^3
} Thermo;


int func (double t, const double y[], double f[], void *params)
{
    // rhs of the equations
    (void)(t); /* avoid unused parameter warning */
    Thermo thermo = *(Thermo *)params;
    f[0] = y[1];
    f[1] = -y[0] - thermo.mu_psis * y[1] * (y[0] * y[0] - 1);
    return GSL_SUCCESS;
}

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
    // jacobian of the RHS
    (void)(t); /* avoid unused parameter warning */
    double mu = *(double *)params;
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 1, 1.0);
    gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
    gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    return GSL_SUCCESS;
}

int main () 
{
    Thermo thermo;
    thermo.gamma = 1.4;
    thermo.mu_psis = 10;
    thermo.rho0_kgpm3 = 1.225;
    gsl_odeiv2_system sys = {func, jac, 2, &thermo};

    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
    double t = 0.0, t1 = 100.0;
    double y[2] = { 1.0, 0.0 };

    // define file output
    char filename[] = "results_ode_sol.csv";
    FILE * fptr = fopen(filename, "w");
    if (fptr == NULL) {
        printf("Error opening the file %s", filename);
        return -1;
    }
    else {
        fprintf(fptr, "# Solution to the second-order nonlinear Van der Pol oscillator equation\n");
        fprintf(fptr, "time, y[0], y[1],\n");
    }

    // stepping through each time
    for (size_t i = 1; i <= 100; i++) {
        double ti = i * t1 / 100.0;
        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
        
        if (status != GSL_SUCCESS) {
            printf("error, return value=%d\n", status);
            break;
        }
        fprintf(fptr, "%.5e, %.5e, %.5e,\n", t, y[0], y[1]);
    }
    fclose(fptr); // close the file output
    gsl_odeiv2_driver_free (d); // free the memory
    return 0;
}