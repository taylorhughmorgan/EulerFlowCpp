// TestBoost.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "SedovBlast.hpp"

struct SedovParams {
    double LenScale__m;     // length scale of the problem
    double DomainLen__m;    // size of the domain
    double PAmb__Pa;        // ambient air pressure
    double PExpl__Pa;       // Explosive pressure
    double RExpl__m;        // radius of explosion
    double tFin__s;         // final simulation time
    double rhoAmb__kgpm3;   // ambient air density
    size_t orders;          // order of solution
    double gamma = 1.4;             // ratio of specific heats
    size_t minNGridPts = 500;       // minimum grid resolution
};

int main()
{
    // testing Sedov/Euler Solution
    // reading in input deck for sedov blast
    std::string finname = "sedov_input.json";
    std::ifstream input_file(finname);

    if (!input_file) {
        std::cerr << "Error opening file " << finname << std::endl;
        return 1;
    }

    json input_deck;
    input_file >> input_deck;

    SedovParams SP = {
        input_deck["sedov"]["Length Scale(m)"],
        input_deck["sedov"]["Domain Size(m)"],
        input_deck["sedov"]["Ambient Pressure(Pa)"],
        input_deck["sedov"]["Explosive Pressure(Pa)"],
        input_deck["sedov"]["Explosive Radius(m)"],
        input_deck["sedov"]["Final Time(s)"],
        input_deck["sedov"]["Ambient Density(kg/m^3)"],
        input_deck["sedov"]["order"],
        input_deck["sedov"]["gamma"],
        input_deck["sedov"]["minNGridPts"]
    };

    SedovBlast Sedov(
        SP.LenScale__m,	    // length scale
        SP.DomainLen__m,	// size of the domain
        SP.RExpl__m,		// radius of explosion
        SP.PExpl__Pa,		// pressure of explosion
        SP.tFin__s,		    // final simulation time
        SP.rhoAmb__kgpm3,   // ambient air density, kg / m ^ 3
        SP.PAmb__Pa,		// ambient air pressure, Pa
        SP.orders,			// order of the equations, 0 = cartesian, 1 - cylindrical, 2 = spherical
        SP.gamma,			// ratio of specific heats, N / A
        SP.minNGridPts      // minimum number of grid points
    );
    Sedov.solve(); // solve the system of equations and convert to SI units
    Sedov.save(input_deck["sedov"]["Output File"], input_deck);
    return 0;
}
