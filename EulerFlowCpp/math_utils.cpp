#include "math_utils.hpp"

void linspace(double start, double end, size_t num, std::vector<double>& linspaced) {
    if (0 != num) {
        if (1 == num) {
            linspaced.push_back(static_cast<double>(start));
        }
        else {
            double delta = (end - start) / (num - 1);

            for (auto i = 0; i < (num - 1); ++i) {
                linspaced.push_back(static_cast<double>(start + delta * i));
            }
            // ensure that start and end are exactly the same as the input
            linspaced.push_back(static_cast<double>(end));
        }
    }
}