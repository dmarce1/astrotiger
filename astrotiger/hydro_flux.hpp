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
//	if (eint < u[egas_i] * 0.001) {
//		eint = std::pow(u[tau_i], opts.gamma);
//	}
	return (opts.gamma - 1.0) * eint;
}

template<class T>
void physical_flux(std::vector<T> &flux, const std::vector<T> u, T v, T p, int dim) {
	for (int f = 0; f < opts.nhydro; f++) {
		flux[f] = u[f] * v;
	}
	flux[sx_i + dim] += p;
	flux[egas_i] += p * v;
}

template<class T>
double hydro_flux(std::vector<T> &flux, const std::vector<T> &ul, const std::vector<T> &ur, int dim) {
	const auto pr = hydro_pressure(ur);
	const auto pl = hydro_pressure(ul);
	const auto rho_r = ur[rho_i];
	const auto rho_l = ul[rho_i];
	const auto rhorinv = 1.0 / rho_r;
	const auto rholinv = 1.0 / rho_l;
	const auto vr = ur[sx_i + dim] * rhorinv;
	const auto vl = ul[sx_i + dim] * rholinv;
	const auto al = std::sqrt(opts.gamma * std::max(pl / rho_l, 1.0e-20));
	const auto ar = std::sqrt(opts.gamma * std::max(pr / rho_r, 1.0e-20));
	const auto sr = std::max(vr + ar, vl + al);
	const auto sl = std::min(vr - ar, vl - al);

	if (0 <= sl) {
		physical_flux(flux, ul, vl, pl, dim);
	} else if (0 >= sr) {
		physical_flux(flux, ur, vr, pr, dim);
	} else {
		const auto s0 = ((pr - pl) + rho_l * vl * (sl - vl) - rho_r * vr * (sr - vr)) / (rho_l * (sl - vl) - rho_r * (sr - vr));
		static thread_local std::vector<T> u0(opts.nhydro);
		static thread_local std::vector<T> f0(opts.nhydro);
		if (0 < s0) {
			const auto rho0 = rho_l * (sl - vl) / (sl - s0);
			u0[rho_i] = rho0;
			u0[sx_i + dim] = rho0 * s0;
			u0[tau_i] = rho0 / rho_l * ul[tau_i];
			u0[pot_i] = rho0 / rho_l * ul[pot_i];
			for (int dim2 = 0; dim2 < NDIM; dim2++) {
				if (dim != dim2) {
					u0[sx_i + dim2] = rho0 * ul[sx_i + dim2] / rho_l;
				}
			}
			u0[egas_i] = rho0 * (ul[egas_i] / rho_l + (s0 - vl) * (s0 + pl / (rho_l * (sl - vl))));
			physical_flux(f0, ul, vl, pl, dim);
			for (int f = 0; f < opts.nhydro; f++) {
				flux[f] = f0[f] + (u0[f] - ul[f]) * sl;
			}
		} else {
			const auto rho0 = rho_r * (sr - vr) / (sr - s0);
			u0[rho_i] = rho0;
			u0[sx_i + dim] = rho0 * s0;
			u0[tau_i] = rho0 / rho_r * ur[tau_i];
			u0[pot_i] = rho0 / rho_r * ur[pot_i];
			for (int dim2 = 0; dim2 < NDIM; dim2++) {
				if (dim != dim2) {
					u0[sx_i + dim2] = rho0 * ur[sx_i + dim2] / rho_r;
				}
			}
			u0[egas_i] = rho0 * (ur[egas_i] / rho_r + (s0 - vr) * (s0 + pr / (rho_r * (sr - vr))));
			physical_flux(f0, ur, vr, pr, dim);
			for (int f = 0; f < opts.nhydro; f++) {
				flux[f] = f0[f] + (u0[f] - ur[f]) * sr;
			}
		}
	}
	return std::max(std::abs(sr), std::abs(sl));

}

//template<class T>
//double hydro_flux(std::vector<T> &flux, const std::vector<T> &ul, const std::vector<T> &ur, int dim) {
//	using namespace std;
//	const auto pr = hydro_pressure(ur);
//	const auto pl = hydro_pressure(ul);
//	const auto rhorinv = 1.0 / ur[rho_i];
//	const auto rholinv = 1.0 / ul[rho_i];
//	const auto vr = ur[sx_i + dim] * rhorinv;
//	const auto vl = ul[sx_i + dim] * rholinv;
//	const auto ar = sqrt(opts.gamma * pr * rhorinv) + abs(vr);
//	const auto al = sqrt(opts.gamma * pl * rholinv) + abs(vl);
//	const auto a = max(ar, al);
//	flux.resize(opts.nhydro);
//	for (int f = 0; f < opts.nhydro; f++) {
//		flux[f] = 0.5 * ((ur[f] * vr + ul[f] * vl)) - a * (ur[f] - ul[f]);
//	}
//	flux[sx_i + dim] += 0.5 * (pr + pl);
//	flux[egas_i] += 0.5 * (pr * vr + pl * vl);
//	return a;
//}

#endif /* ASTROTIGER_HYDRO_FLUX_HPP_ */
