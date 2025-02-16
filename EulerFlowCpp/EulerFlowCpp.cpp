// TestBoost.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iterator>
#include <algorithm>
#include "BoundaryConditions.hpp"

int main()
{
    // creating an empty array to test out boundary conditions
    pde_state u0, grid;
    size_t nEle = 10;
    for (size_t i = 0; i < nEle; i++) {
        u0.push_back(i);
    }
    std::vector<size_t> left_ids = {0, 1};
    std::vector<size_t> right_ids = { nEle, nEle - 1 };
    // using std::functions 
    auto bc_left = agnosticBCs(left_ids, validBCs::CONSTANT, 10.0);
    auto bc_right = agnosticBCs(right_ids, validBCs::GRADIENT);

    return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
