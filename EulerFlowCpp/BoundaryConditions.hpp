#pragma once
#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <functional>

// type definition for PDE state
typedef std::vector<double> pde_state;

// function pointer for applying boundary conditions
using BCFunction = std::function<void(pde_state&, const pde_state&)>;

// list of valid Boundary Conditions
enum class validBCs { REFLECTIVE, GRADIENT, CONSTANT, EXTRAPOLATED };

BCFunction agnosticBCs(const std::vector<size_t>& ids, validBCs bc_id, double bc_value = 0.0);