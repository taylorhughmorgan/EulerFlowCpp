#ifndef MATH_UTILS_CPP
#define MATH_UTILS_CPP
#include "math_utils.hpp"

template <typename T>
void linspace(double start, double end, size_t num, std::vector<T>& linspaced) {
    if (0 != num) {
        if (1 == num) {
            linspaced.push_back(static_cast<T>(start));
        }
        else {
            T delta = (end - start) / (num - 1);

            for (auto i = 0; i < (num - 1); ++i) {
                linspaced.push_back(static_cast<T>(start + delta * i));
            }
            // ensure that start and end are exactly the same as the input
            linspaced.push_back(static_cast<T>(end));
        }
    }
}

#endif