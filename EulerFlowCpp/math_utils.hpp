#ifndef MATH_UTILS_HPP
#define MATH_UTILS_HPP
#include <vector>

template <typename T>
void linspace(double start, double end, size_t num, std::vector<T>& linspaced);

#include "math_utils.cpp"

#endif