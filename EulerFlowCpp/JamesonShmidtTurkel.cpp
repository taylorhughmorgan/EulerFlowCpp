#include "JamesonShmidtTurkel.hpp"

void d3dx3_fwd(pde_state& y, pde_state& d3ydx3) {
	/* Third-order forward spatial differencing */
	size_t ny = y.size();
	for (size_t i = 0; i < d3ydx3.size() - 1; i++) {
		d3ydx3[i] = y[i + 3] - 3 * y[i + 2] + 3 * y[i + 1] - y[i];
	}
	// backward difference on last cell
	d3ydx3.back() = y[ny - 1] - 3 * y[ny - 2] + 3 * y[ny - 3] - y[ny - 4];
}

void d3dx3_bckwrd(pde_state& y, pde_state& d3ydx3) {
	/* Third-order backward spatial differencing */
	// forward differencing on the first two cells
	d3ydx3[0] = y[3] - 3 * y[2] + 3 * y[1] - y[0];
	d3ydx3[1] = y[4] - 3 * y[3] + 3 * y[2] - y[1];
	for (size_t i = 2; i < d3ydx3.size(); i++) {
		d3ydx3[i] = y[i + 1] - 3 * y[i] + 3 * y[i - 1] - y[i - 2];
	}
}