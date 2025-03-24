#ifndef MATH_UTILS_HPP
#define MATH_UTILS_HPP
#include <vector>
#include <boost/math/constants/constants.hpp>

// define pi
static const double pi = boost::math::constants::pi<double>();

template <typename T>
void linspace(double start, double end, size_t num, std::vector<T>& linspaced);

#include "math_utils.cpp"

#endif