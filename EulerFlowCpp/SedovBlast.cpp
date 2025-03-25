#include "SedovBlast.hpp"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/numeric/odeint.hpp>
#include <chrono>

using namespace boost::numeric::odeint;

// for writing the output to file
void write_to_csv(const std::string& filename, const std::vector<std::vector<double>>& data) {
	std::ofstream file(filename);
	if (!file.is_open()) {
		std::cerr << "Failed to open file: " << filename << std::endl;
		return;
	}

	for (const auto& row : data) {
		for (size_t i = 0; i < row.size(); ++i) {
			file << row[i];
			if (i < row.size() - 1) file << ","; // Separate values with commas
		}
		file << "\n"; // New row
	}

	file.close();
	std::cout << "Data written to " << filename << std::endl;
}

/**************** OBSERVER CLASS *************************************/
struct push_back_state_and_time
{
	/* In order to observe the solution during the integration steps all you have to do is to provide 
	a reasonable observer, which stores the intermediate steps in a container */
	std::vector< pde_state >& m_states;
	std::vector< double >& m_times;
	size_t m_grid_size;

	void det_stop_state(const pde_state& x) {
		// determine if a threshold has been reached and stop the integration
		// check if rho, rho*U and rho*E are greater than the previous iteration. or if the end points contain NaNs
		double rhoNew_end = x[m_grid_size - 1];
		double rhoOld_end = m_states.back()[m_grid_size - 1];

		if (rhoNew_end > rhoOld_end) {
			throw std::runtime_error("rho(new) > rho(old) @ end of domain.");
		}
		else if (std::isnan(rhoNew_end)) {
			throw std::runtime_error("NaN value identified at end of domain.");
		}
	}

	push_back_state_and_time(std::vector< pde_state >& states, std::vector< double >& times, const size_t grid_size)
		: m_states(states), m_times(times), m_grid_size(grid_size) {
	}

	void operator()(const pde_state& x, double t)
	{
		// append to tracked states
		m_states.push_back(x);
		m_times.push_back(t);
		// check if needing to abort simulation
		det_stop_state(x);
	}
};

/*********************** SEDOV BLAST SOLUTION CLASS ***********************/
SedovBlast::SedovBlast(
	double m_ScaleLen__m,	// length scale
	double m_DomainLen__m,	// size of the domain
	double m_RExpl__m,		// radius of explosion
	double m_PExpl__Pa,		// pressure of explosion
	double m_tFinal__s,		// final simulation time
	double m_rho0__kgpm3,	// ambient air density, kg / m ^ 3
	double m_P0__Pa,		// ambient air pressure, Pa
	size_t m_order,			// order of the equations, 0 = cartesian, 1 - cylindrical, 2 = spherical
	double m_gamma,			// ratio of specific heats, N / A
	size_t m_minNGridPts	// minimum number of grid points
)
{
	// Convert the parameters of the Sedov Blast to nondimensional form, for speed and numerical stability.
	this->ScaleLen__m = m_ScaleLen__m;
	this->DomainLen__m = m_DomainLen__m;
	this->RExpl__m = m_RExpl__m;
	this->PExpl__Pa = m_PExpl__Pa;
	this->tFinal__s = m_tFinal__s;
	this->rho0__kgpm3 = m_rho0__kgpm3;
	this->P0__Pa = m_P0__Pa;
	this->order = m_order;
	this->gamma = m_gamma;
	this->minNGridPts = m_minNGridPts;

	// Calculate the internal energy of the explosion
	size_t n = order + 1;
	EExpl__J = pow(pi, double(n) / 2.0) * boost::math::tgamma(double(n) / 2.0 + 1.0) * pow(RExpl__m, n);

	// dimensionless parameters: scale using rho0, P0, and diameter
	UScale = sqrt(P0__Pa / rho0__kgpm3);
	lenStar = DomainLen__m / ScaleLen__m;
	double rExplStar = RExpl__m / ScaleLen__m;
	double pExplStar = PExpl__Pa / P0__Pa;
	double tFinStar = tFinal__s * UScale / ScaleLen__m;

	// set up the radial grid, we want at least 10 points for the explosion
	nGridPts = size_t(std::ceil(lenStar / rExplStar * 10.0));
	nGridPts = std::max(nGridPts, minNGridPts);

	double rMinStar = std::min(rExplStar / 10.0, lenStar / 100.0);
	linspace<double>(rMinStar, lenStar, nGridPts, grid);
	linspace<double>(0.0, tFinStar, minNGridPts, times);

	dr = grid[1] - grid[0];
	dt = times[1] - times[0];

	// setting the initial conditions - setting rho0Star and p0Star to 1.0 for dimensionless variables
	rho0Star.resize(nGridPts);
	p0Star.resize(nGridPts);
	v0Star.resize(nGridPts);
	std::fill(rho0Star.begin(), rho0Star.end(), 1.0);
	std::fill(p0Star.begin(), p0Star.end(), 1.0);
	std::fill(v0Star.begin(), v0Star.end(), 0.0);

	for (size_t i = 0; i < nGridPts; i++) {
		if (grid[i] < rExplStar) {
			p0Star[i] *= pExplStar;
		}
		else break;
	}

	// time and grid in dimensional/metric scale
	for (size_t i = 0; i < nGridPts; i++) {
		rGrid__m.push_back(grid[i] * ScaleLen__m);
	}
}


void SedovBlast::solve()
{
	//
	std::array<double, 2> alpha = { 0.5, 0.5 };
	std::array<double, 2> beta = { 0.25, 0.5 };
	bc_map boundary_conds = { {Fields::RHO, {validBCs::GRADIENT, validBCs::GRADIENT} },
						{Fields::VEL, {validBCs::REFLECTIVE, validBCs::GRADIENT} },
						{Fields::ENERGY, {validBCs::GRADIENT, validBCs::GRADIENT} } };

	// system of equations
	ODEs = std::make_unique<EulerSol>(grid, boundary_conds, order, alpha, beta, gamma);
	// intial conditions
	pde_state W0Star;
	ODEs->createICs(rho0Star, v0Star, p0Star, W0Star);

	std::cout << "Solving the Euler Equation as a system of ODES. \n" <<
		"t_range = [" << 0.0 << ", " << times.back() << "](dimensionless) \n" <<
		"nGridPts = " << nGridPts << std::endl <<
		"r_range = [" << grid[0] << "," << grid.back() << "](dimensionless)\n";

	// define observer
	std::vector<pde_state> states_sol;
	pde_state times_sol;
	push_back_state_and_time observer(states_sol, times_sol, ODEs->size);

	// define numerical stepper
	auto stepper = make_controlled(1e-6, 1e-6, runge_kutta_cash_karp54<pde_state>()); //runge_kutta_cash_karp54<pde_state>();

	// time the execution
	auto start = std::chrono::high_resolution_clock::now();

	// run the integration
	try {
		integrate_times(stepper, (*ODEs), W0Star, times.begin(), times.end(), dt, observer);
	}
	catch (const std::runtime_error &e) {
		std::cout << "Integration stopped prematurely at " << times_sol.back() << ". Reason: " << e.what() << std::endl;
	}

	// stop the timer
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

	std::cout << "Solution reached after " << duration.count() << "miliseconds\n";

	// convert solution to primatives
	std::vector<pde_state> rhoStar_sol, pStar_sol, EStar_sol, uStar_sol;
	ODEs->convAllPrimatives(states_sol, rhoStar_sol, uStar_sol, EStar_sol, pStar_sol);
	std::cout << "Finished converting ODE solution to primative fields\n";

	// converting from dimensionless to dimensional fields
	rho_sol.resize(rhoStar_sol.size());
	p_sol.resize(pStar_sol.size());
	E_sol.resize(EStar_sol.size()); 
	u_sol.resize(uStar_sol.size());

	// looping through each state
	for (size_t iState = 0; iState < rhoStar_sol.size(); iState++)
	{
		rho_sol[iState].resize(rhoStar_sol[iState].size());
		p_sol[iState].resize(pStar_sol[iState].size());
		E_sol[iState].resize(EStar_sol[iState].size());
		u_sol[iState].resize(uStar_sol[iState].size());

		// converting each state value
		for (size_t i = 0; i < rhoStar_sol[iState].size(); i++) {
			rho_sol[iState][i] = rhoStar_sol[iState][i] * rho0__kgpm3;
			u_sol[iState][i] = uStar_sol[iState][i] * sqrt(P0__Pa / rho0__kgpm3);
			p_sol[iState][i] = pStar_sol[iState][i] * P0__Pa;
			E_sol[iState][i] = p_sol[iState][i] / (rho_sol[iState][i] * (gamma - 1.0)) + 0.5 * pow(u_sol[iState][i], 2);
		}
	}

	// creating the time grid
	for (size_t i = 0; i < times_sol.size(); i++) {
		tGrid__s.push_back(times_sol[i] * ScaleLen__m / UScale);
	}
	std::cout << "Converted from dimensionless to dimensional fields\n";
}


void SedovBlast::save(std::string foutname, json inputs)
{
	// save the output to a JSON file with the large arrays saved as CSV files
	std::string rho_fname = "rho_sol.csv";
	std::string press_fname = "press_sol.csv";
	std::string vel_fname = "vel_sol.csv";
	std::string energy_fname = "energy_sol.csv";
	
	nlohmann::ordered_json output = inputs;
	output["fields"]["density(kg/m^3)"] = rho_fname;
	output["fields"]["pressure(Pa)"] = press_fname;
	output["fields"]["velocity(m/s)"] = vel_fname;
	output["fields"]["energy(J)"] = energy_fname;
	output["grids"]["times(s)"] = tGrid__s;
	output["grids"]["grid(m)"] = rGrid__m;

	// writing to field data to files
	write_to_csv(rho_fname, rho_sol);
	write_to_csv(press_fname, p_sol);
	write_to_csv(vel_fname, u_sol);
	write_to_csv(energy_fname, E_sol);

	std::ofstream fout(foutname);
	if (fout.is_open()) {
		fout << output.dump(4); // dump with 4 spaces of indentation
		fout.close();
		std::cout << "Euler Solution written to " << foutname << std::endl;
	}
	else {
		std::cerr << "Error writting to file " << foutname << std::endl;
	}
}