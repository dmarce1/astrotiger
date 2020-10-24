/*
 * hydro_flux.hpp
 *
 *  Created on: Oct 19, 2020
 *      Author: dmarce1
 */

#ifndef ASTROTIGER_HYDRO_FLUX_HPP_
#define ASTROTIGER_HYDRO_FLUX_HPP_

#include <astrotiger/hydro_grid.hpp>
#include <astrotiger/options.hpp>
#include <vector>
#include <cmath>

template<class T>
inline double hydro_kinetic(const std::vector<T> &u) {
	T ekin = T(0.0);
	for (int dim = 0; dim < NDIM; dim++) {
		ekin += 0.5 * u[sx_i + dim] * u[sx_i + dim] / u[rho_i];
	}
	return ekin;
}

template<class T>
inline double hydro_pressure(const std::vector<T> &u) {
	using namespace std;
	auto eint = max(T(0), u[egas_i] - hydro_kinetic(u));
	if (eint < u[egas_i] * 0.001) {
		eint = std::pow(u[tau_i], opts.gamma);
	}
	return (opts.gamma - 1.0) * eint;
}

template<class T>
double hydro_flux(std::vector<T> &flux, const std::vector<T> &ul, const std::vector<T> &ur, int dim) {
	using namespace std;
	const auto pr = hydro_pressure(ur);
	const auto pl = hydro_pressure(ul);
	const auto rhorinv = 1.0 / ur[rho_i];
	const auto rholinv = 1.0 / ul[rho_i];
	const auto vr = ur[sx_i + dim] * rhorinv;
	const auto vl = ul[sx_i + dim] * rholinv;
	const auto ar = sqrt(opts.gamma * pr * rhorinv) + abs(vr);
	const auto al = sqrt(opts.gamma * pl * rholinv) + abs(vl);
	const auto a = max(ar, al);
	flux.resize(opts.nhydro);
	for (int f = 0; f < opts.nhydro; f++) {
		flux[f] = 0.5 * ((ur[f] * vr + ul[f] * vl)) - a * (ur[f] - ul[f]);
	}
	flux[sx_i + dim] += 0.5 * (pr + pl);
	flux[egas_i] += 0.5 * (pr * vr + pl * vl);
	return a;
}

#endif /* ASTROTIGER_HYDRO_FLUX_HPP_ */
