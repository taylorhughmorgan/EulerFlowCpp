// TestBoost.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "SedovBlast.hpp"

void testBoundaryConditions() {
    // creating an empty array to test out boundary conditions
    pde_state u0, grid;
    size_t nEle = 10;
    for (size_t i = 0; i < nEle; i++) {
        u0.push_back(std::pow(double(i), 2));
        grid.push_back(double(i));
    }
    std::array<size_t, 3> left_ids = { 0, 1, 2 };
    std::array<size_t, 3> right_ids = { nEle - 1, nEle - 2, nEle - 3 };

    // using std::functions, create boundary conditions
    BCFunction bc_left = agnosticBCs(left_ids, validBCs::CONSTANT, 10.0);
    BCFunction bc_right = agnosticBCs(right_ids, validBCs::EXTRAPOLATED);

    // print out the grid and u0 prior to applying BCs
    std::cout << "grid (m), u0 (K) " << std::endl;
    for (size_t i = 0; i < nEle; i++) {
        std::cout << grid[i] << ", " << u0[i] << std::endl;
    }

    // print out grid and u0 post-appling bcs
    bc_left(u0, grid);
    bc_right(u0, grid);

    std::cout << "Post-BC application:\ngrid (m), u0 (K) " << std::endl;
    for (size_t i = 0; i < nEle; i++) {
        std::cout << grid[i] << ", " << u0[i] << std::endl;
    }
}

int main()
{
    // testing boundary conditions
    testBoundaryConditions();

    // testing Sedov/Euler Solution
    double LenScale__m = 1;         // length scale of the problem
    double DomainLen__m = 10;       // size of the domain
    double PAmb__Pa = 101325;       // ambient air pressure
    double PExpl__Pa = 40 * PAmb__Pa; // Explosive pressure
    double RExpl__m = 3;            // radius of explosion
    double tFin__s = 0.010;         // final simulation time
    double rhoAmb__kgpm3 = 1.225;   // ambient air density
    size_t orders = 0;              // order of solution
    double gamma = 1.4;             // ratio of specific heats
    size_t minNGridPts = 500;       // minimum grid resolution

    SedovBlast Sedov(
        LenScale__m,	// length scale
        DomainLen__m,	// size of the domain
        RExpl__m,		// radius of explosion
        PExpl__Pa,		// pressure of explosion
        tFin__s,		// final simulation time
        rhoAmb__kgpm3,  // ambient air density, kg / m ^ 3
        PAmb__Pa,		// ambient air pressure, Pa
        orders,			// order of the equations, 0 = cartesian, 1 - cylindrical, 2 = spherical
        gamma,			// ratio of specific heats, N / A
        minNGridPts     // minimum number of grid points
    );
    Sedov.solve();

    return 0;
}
