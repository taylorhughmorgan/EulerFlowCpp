// TestBoost.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iterator>
#include <algorithm>
#include <cmath>
#include "SedovEuler.hpp"

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
    size_t orders = 2;              // order of solution
    double gamma = 1.4;             // ratio of specific heats

    SedovBlast Sedov(LenScale__m, DomainLen__m, RExpl__m, PExpl__Pa, tFin__s,
        PAmb__Pa, rhoAmb__kgpm3, orders, gamma);
    Sedov.solve();

    return 0;
}
