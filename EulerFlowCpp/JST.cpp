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

void JST_2ndOrderEulerFlux(flux_state& fluxVector, flux_state& dissipationFlux) {
	/* Calculate second order Euler Flux */
	for (size_t iflux = 0; iflux < fluxVector.size(); iflux++) {
		for (size_t i = 1; i < fluxVector[iflux].size() - 1; i++) {
			double hbar_jphalf = (fluxVector[iflux][i + 1] + fluxVector[iflux][i] ) / 2;
			double hbar_jmhalf = (fluxVector[iflux][i] + fluxVector[iflux][i - 1]) / 2;
			dissipationFlux[iflux][i - 1] = hbar_jphalf - hbar_jmhalf;
		}
	}
}


void JST_DissipationFlux(
	flux_state& D_j,				// dissipation flux
	flux_state& W,					// state vector variable
	pde_state& p,					// pressure
	pde_state& u,					// velocity
	pde_state& cs,					// speed of sound
	std::array<double, 2>& alpha,	// dissipative flux terms for spatial differencing
	std::array<double, 2>& beta		// dissipative flux terms for spatial differencing
) {
	/* Calculate dissipation flux using JST method */
	// nu_j - pressure dissipation term: central differencing
	size_t n_nuj = p.size();
	pde_state nu_j, R_jphalf, R_jmhalf;
	nu_j.resize(n_nuj);

	nu_j[0] = std::abs((p[2] - 2 * p[1] + p[0]) / (p[2] + 2 * p[1] + p[0]));
	for (size_t i = 1; i < p.size() - 1; i++) {
		nu_j[i] = std::abs((p[i + 1] - 2 * p[i] + p[i - 1]) / (p[i + 1] + 2 * p[i] + p[i - 1]));
	}
	nu_j.back() = std::abs((p[n_nuj - 1] - 2 * p[n_nuj - 2] + p[n_nuj - 3]) / (p[n_nuj - 1] + 2 * p[n_nuj - 2] + p[n_nuj - 3]));

	// maximum wave speed of the system
	R_jmhalf.resize(cs.size() - 2);
	R_jphalf.resize(cs.size() - 2);
	for (size_t i = 1; i < cs.size() - 1; i++) {
		R_jphalf[i - 1] = (std::abs(u[i + 1]) + cs[i + 1] + std::abs(u[i]) + cs[i]) / 2.0;
		R_jmhalf[i - 1] = (std::abs(u[i]) + cs[i] + std::abs(u[i - 1]) + cs[i - 1]) / 2.0;
	}

	for (size_t iflux = 0; iflux < W.size(); iflux++) {
		// third derivative of the state vector, fwd differencing for j+1/2 and backward for j-1/2
		pde_state delta3W_jphalf(W[iflux].size()-2), delta3W_jmhalf(W[iflux].size() - 2);
		d3dx3_fwd(W[iflux], delta3W_jphalf);
		d3dx3_bckwrd(W[iflux], delta3W_jmhalf);

		for (size_t i = 1; i < W[iflux].size() - 1; i++) {
			// first derivative of the state vector
			double deltaW_jphalf = W[iflux][i + 1] - W[iflux][i];
			double deltaW_jmhalf = W[iflux][i] - W[iflux][i - 1];

			// dissipative coefficients, S_(j-1/2) MIGHT NEED FIXING
			double S_jphalf = std::max(nu_j[i + 1], nu_j[i]); //std::max(nu_j[i], nu_j[i - 1]);
			double S_jmhalf = std::max(nu_j[i], nu_j[i - 1]);

			double eps2_jphalf = std::min(alpha[0], alpha[1] * S_jphalf);
			double eps2_jmhalf = std::min(alpha[0], alpha[1] * S_jmhalf);
			double eps4_jphalf = std::max(0.0, beta[0] - beta[1] * eps2_jphalf);
			double eps4_jmhalf = std::max(0.0, beta[0] - beta[1] * eps2_jmhalf);

			// artificial dissipation terms
			double d_jphalf = eps2_jphalf * R_jphalf[i - 1] * deltaW_jphalf - eps4_jphalf * R_jphalf[i - 1] * delta3W_jphalf[i - 1];
			double d_jmhalf = eps2_jmhalf * R_jmhalf[i - 1] * deltaW_jmhalf - eps4_jmhalf * R_jmhalf[i - 1] * delta3W_jmhalf[i - 1];
			D_j[iflux][i - 1] = d_jphalf - d_jmhalf;
		}
	}
}