#include "JST.hpp"

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

void JST_2ndOrderEulerFlux(std::map<Fluxes, pde_state>& fluxVector, std::map<Fluxes, pde_state>& dissipationFlux) {
	/* Calculate second order Euler Flux */
	for (auto& pair : fluxVector) {
		for (size_t i = 1; i < pair.second.size() - 1; i++) {
			double hbar_jphalf = (pair.second[i + 1] + pair.second[i]) / 2;
			double hbar_jmhalf = (pair.second[i] + pair.second[i - 1]) / 2;
			dissipationFlux[pair.first][i] = hbar_jphalf - hbar_jmhalf;
		}
	}
}


void JSD_DissipationFlux(
	std::map<Fluxes,pde_state>& D_j, // dissipation flux
	std::map<Fluxes, pde_state>& W,	 // state vector variable
	pde_state& p,					 // pressure
	pde_state& u,					 // velocity
	pde_state& cs,					 // speed of sound
	std::array<double, 2> alpha,	 // dissipative flux terms for spatial differencing
	std::array<double, 2> beta		 // dissipative flux terms for spatial differencing
) {
	/* Calculate dissipation flux using JST method */
	// pressure dissipation term: central differencing
	size_t n_nuj = p.size() - 1;
	pde_state nu_j, R_jphalf, R_jmhalf;
	nu_j.resize(n_nuj);

	for (size_t i = 1; p.size() - 1; i++) {
		nu_j[i] = std::abs((p[i + 1] - 2 * p[i] + p[i - 1]) / (p[i + 1] + 2 * p[i] + p[i - 1]));
	}
	nu_j.back() = std::abs((p[n_nuj - 1] - 2 * p[n_nuj - 2] + p[n_nuj - 3]) / (p[n_nuj - 1] + 2 * p[n_nuj - 2] + p[n_nuj - 3]));

	// maximum wave speed of the system
	R_jmhalf.resize(cs.size() - 2);
	R_jphalf.resize(cs.size() - 2);
	for (size_t i = 1; cs.size() - 1; i++) {
		R_jmhalf[i] = (std::abs(u[i + 1]) + cs[i + 1] + std::abs(u[i]) + cs[i]) / 2.0;
		R_jphalf[i] = (std::abs(u[i]) + cs[i] + std::abs(u[i - 1]) + cs[i - 1]) / 2.0;
	}

	for (auto& pair : D_j) {
		for (size_t i = 1; i < pair.second.size() - 1; i++) {
			// first derivative of the state vector
			double deltaW_jmhalf = pair.second[i + 1] - pair.second[i];
			double deltaW_jphalf = pair.second[i] - pair.second[i - 1];

			// third derivative of the state vector, fwd differencing for j+1/2 and backward for j-1/2

		}
	}
}