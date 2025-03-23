#pragma once
#include <string>
#include <vector>
#include <array>
#include <stdexcept>
#include <iostream>
#include <functional>
#include <boost/math/constants/constants.hpp>

// define pi
static const double pi = boost::math::constants::pi<double>();

// type definition for PDE state
typedef std::vector<double> pde_state;

// type definition for bc_ids
typedef std::array<size_t, 3> ident;

// function pointer for applying boundary conditions
using BCFunction = std::function<void(pde_state&, const pde_state&)>;

// list of valid Boundary Conditions
enum class validBCs { REFLECTIVE, GRADIENT, CONSTANT, EXTRAPOLATED };

BCFunction agnosticBCs(const ident& ids, validBCs bc_id, double bc_value = 0.0);