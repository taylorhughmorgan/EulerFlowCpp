/* Solving the Compressible Euler equations using GSL's gsl_odeiv2_step_rk8pd */
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_statistics.h>
#include "pdes.h"
#include "grids.h"

void array_copy(const gsl_block * src, double * dest) {
    // copy values from one array to another
    for (size_t i = 0; i < src->size; ++i) {
        dest[i] = src->data[i];
    }
}

// structure for containing grids and ghost grids
typedef struct
{
    gsl_block * rho;
    gsl_block * rho_U;
    gsl_block * rho_E;
    gsl_block * u;
    gsl_block * E;
    gsl_block * p;
    gsl_block * ghost_rho;
    gsl_block * ghost_u;
    gsl_block * ghost_E;
} Fields;

void alloc_fields(Fields* fields, size_t size) {
    // initialize the standard and ghost grids
    fields->rho     = gsl_block_alloc(size);
    fields->rho_U   = gsl_block_alloc(size);
    fields->rho_E   = gsl_block_alloc(size);
    fields->u       = gsl_block_alloc(size);
    fields->E       = gsl_block_alloc(size);
    fields->p       = gsl_block_alloc(size);
    fields->ghost_rho = gsl_block_alloc(size + 2);
    fields->ghost_u = gsl_block_alloc(size + 2);
    fields->ghost_E = gsl_block_alloc(size + 2);
}

void free_fields(Fields* fields) {
    // free the memory
    gsl_block_free(fields->rho);
    gsl_block_free(fields->rho_U);
    gsl_block_free(fields->rho_E);
    gsl_block_free(fields->u);
    gsl_block_free(fields->E);
    gsl_block_free(fields->p);
    gsl_block_free(fields->ghost_rho);
    gsl_block_free(fields->ghost_u);
    gsl_block_free(fields->ghost_E);
}

// structure containing grid and field variables
typedef struct {
    double gamma;
    double mu;
    double dx;
    double dt;
    double max_wave_speed;
    size_t size, ghost_size;
    gsl_block * grid;
    gsl_block * ghost_grid;
    gsl_block * grid_to_order;
    Fields * fields;
} Params;

void alloc_params(Params * params, size_t size) {
    // allocate grids and fields in params
    params->size = size;
    params->grid = gsl_block_alloc(size);
    params->ghost_grid = gsl_block_alloc(size + 2);
    params->grid_to_order = gsl_block_alloc(size);
    alloc_fields(params->fields, size);
}

void free_params(Params * params) {
    // free params and fields
    gsl_block_free(params->grid);
    gsl_block_free(params->ghost_grid);
    gsl_block_free(params->grid_to_order);
    free_fields(params->fields);
}

// function for applying BCs
void apply_BC(BC1D * bc, double* y, Params* sim_params) {
    bc->left_bc(y, 0.0, sim_params->dx, 0);
    bc->right_bc(y, 0.0, sim_params->dx, sim_params->size - 1);    
}

int func(double t, const double y[], double dydt[], void *params)
{
    /* RHS of the heat equation:
        \partial U / \partial t =  */
    (void)(t); /* avoid unused parameter warning */
    // parsing input parameters
    Params self = *(Params *)params;

    // instead of temporary variables, access pre-allocated fields
    double * rho   = &(self.fields->rho->data);
    double * rho_U = &(self.fields->rho_U->data);
    double * rho_E = &(self.fields->rho_E->data);
    double * U     = &(self.fields->u->data);
    double * E     = &(self.fields->E->data);
    // defining ghost grid:
    double ghost_rho[self.size + 2], ghost_U[self.size + 2], ghost_E[self.size + 2]; // density, velocity, internal energy
    double ghost_p[self.size + 2], ghost_H[self.size + 2], ghost_cs[self.size + 2];  // pressure, enthalpy, sound speed
    // defining state and flux vectors
    double W1[self.size + 2], W2[self.size + 2], W3[self.size + 2];
    double F1[self.size + 2], F2[self.size + 2], F3[self.size + 2];

    // unraveling rho, rho*U, and rho*E
    for (size_t i = 0; i < self.size; ++i) {
        rho[i] = y[i] / self.grid_to_order->data[i];
        rho_U[i] = y[i + self.size] / self.grid_to_order->data[i];
        rho_E[i] = y[i + 2 * self.size] / self.grid_to_order->data[i];

        // getting primatives U and E by dividing by rho
        U[i] = rho_U[i] / rho[i];
        E[i] = rho_E[i] / rho[i];

        // fill ghost cells starting with index=1
        ghost_rho[i + 1] = rho[i];
        ghost_U[i + 1] = U[i];
        ghost_E[i + 1] = E[i];
    }

    // apply boundary conditions

    // apply equations of state on ghost-grid
    for (size_t iG = 0; iG < self.ghost_size; ++iG) {
        ghost_p[iG] = ghost_rho[iG] * (self.gamma - 1.0) * (ghost_E[iG] - 0.5 * pow(ghost_U[iG], 2));
        ghost_H[iG] = ghost_E[iG] + ghost_p[iG] / ghost_rho[iG];
        ghost_cs[iG] = sqrt(self.gamma * ghost_p[iG] / ghost_rho[iG]);

        // develop W: state vector variable
        // develop F: flux vector variable
        // develop S: source term variable
    }

    // calculating second-order Euler flux and dissipation flux

    // calculating residuals

    // loading residuals into dy/dt
    for (size_t i = 0; i < self.size; ++i) {
        dydt[i] = 0.0;
        dydt[i + self.size] = 0.0;
        dydt[i + 2 * self.size] = 0.0;
    }

    return GSL_SUCCESS;
}

gsl_block * create_ICs(double * rho0, double * press0, double * vel0, size_t size) {
    // create initial conditions for flux-vector state
    gsl_block * y;
    y = gsl_block_alloc(3 * size);
    return y;
}

void calc_dt(Params * sim) {
    // calculate max sound speed, wave speed and time step
    double courant = 0.5;
    
    double * press = sim->fields->p->data;
    double * rho = sim->fields->rho->data;
    double max_sound_speed = sqrt( sim->gamma * gsl_stats_max(press, 1, sim->size) / gsl_stats_min(rho, 1, sim->size) );
    sim->max_wave_speed = abs( gsl_stats_max(sim->fields->u->data, 1, sim->size) );
    sim->dt = courant * sim->dx / sim->max_wave_speed;
}

int euler_eqns(CoordSys coord,  // coordinate system
    double * rho0,              // initial density
    double * press0,            // initial pressure
    double * vel0,              // initial velocity
    size_t grid_size,           // grid size
    double domain_len,          // length of the domain
    BC1D * rho_bc,              // density boundary condition
    BC1D * vel_bc,              // velocity boundary condition
    BC1D * press_bc,            // pressure boundary condition
    double tfinal              // final simulation time
) {
    // Integrate the Euler equations in 1D speherical, cartesial, or polar coordinates
    // intializing the thermodynamic parameters and grid size
    Params sim;
    alloc_params(&sim, grid_size);
    sim.gamma = 1.4;

    // setting up the grid
    double xmin = 0.0, xmax = domain_len;
    double t = 0.0;
    sim.dx = (xmax - xmin) / (double)sim.size;
    // calculating max wave speed, which drives time step
    array_copy(sim.fields->p, press0);
    array_copy(sim.fields->rho, rho0);
    array_copy(sim.fields->u, vel0);
    calc_dt(&sim);

    // creating initial conditions
    gsl_block * y = create_ICs(rho0, press0, vel0, sim.size);

    // define the ODE solver
    printf("Setting up the system of equations\n");
    gsl_odeiv2_system sys = {func, NULL, sim.size, &sim};
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);

    // define file output
    char filename[] = "euler_eqn_sol.csv";
    FILE * fptr = fopen(filename, "w");
    if (fptr == NULL) {
        printf("Error opening the file %s", filename);
        return -1;
    }
    else {
        printf("Setting up the output file\n");
        fprintf(fptr, "# Solution to the heat equation for gamma=%f on grid xmin=%f, xmax=%f\n", sim.gamma, xmin, xmax);
        fprintf(fptr, "time, ");
        for (size_t i = 0; i < sim.size; ++i) {
            fprintf(fptr, "y[%zu], ", i);
        }
        fprintf(fptr, "\n");
    }
    
    // stepping through each time
    double output_tracker = 0.0;

    printf("Begin looping through solution\n");
    while (t < tfinal) {
        output_tracker += sim.dt;
        double t_f = t + sim.dt;
        int status = gsl_odeiv2_driver_apply(d, &t, t_f, y);
        
        // update time step
        calc_dt(&sim);
        if (status != GSL_SUCCESS) {
            printf("error, return value=%d\n", status);
            break;
        }
    }

    printf("Solution reached!\nSolution saved to %s", filename);
    fclose(fptr); // close the file output
    gsl_odeiv2_driver_free (d); // free the memory
    free_params(&sim);
    return GSL_SUCCESS;
}