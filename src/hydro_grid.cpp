/*
 * grid.cpp
 *
 *  Created on: Oct 14, 2020
 *      Author: dmarce1
 */

#include <astrotiger/hydro_grid.hpp>
#include <astrotiger/hydro_flux.hpp>
#include <astrotiger/multi_array.hpp>
#include <astrotiger/polytrope.hpp>
#include <astrotiger/cosmos.hpp>
#include <astrotiger/tree.hpp>
#include <astrotiger/rand.hpp>
#include <astrotiger/fileio.hpp>
#include <astrotiger/chemistry.hpp>
#include <vector>

hydro_grid::hydro_grid() {
}

energy_statistics hydro_grid::get_energy_statistics(const multi_array<double> &phi, const std::vector<multi_range> &exclude, double rho0) const {
	const auto a = cosmos_a();
	energy_statistics e;
	e.ekin = 0.0;
	e.epot = 0.0;
	const auto dv = std::pow(dx, NDIM);
	for (multi_iterator i(box.pad(-opts.hbw)); !i.end(); i++) {
		bool use = true;
		for (const auto &c : exclude) {
			if (c.contains(i)) {
				use = false;
				break;
			}
		}
		if (use) {
			e.ekin += U[egas_i][i] * dv * std::pow(a, NDIM * (1 - opts.gamma));
			e.epot += U[rho_i][i] * 0.5 * phi[i] * dv;
		}
	}
	return e;
}

void hydro_grid::to_array(multi_array<double> &a, const multi_range &bbox, int f, double w) const {
	for (multi_iterator i(bbox); !i.end(); i++) {
		a[i] = w * U[f][i] + (1.0 - w) * U0[f][i];
	}
}

void hydro_grid::store() {
	U0 = U;
}

void hydro_grid::store_flux() {
	F0 = F;
}

double hydro_grid::positivity_limit() const {
	double amax = 0.0;
	const double a = cosmos_a();
	for (multi_iterator i(box.pad(-opts.hbw)); !i.end(); i++) {
		double df_rho = 0.0;
		double df_tau = 0.0;
		for (int dim = 0; dim < NDIM; dim++) {
			multi_index ip1 = i;
			ip1[dim]++;
			df_rho += (F[dim][rho_i][ip1] - F[dim][rho_i][i]) * a;
			df_tau += (F[dim][tau_i][ip1] - F[dim][tau_i][i]) * a;
		}
		amax = std::max(amax, std::max(df_rho, 0.0) / U[rho_i][i]);
		amax = std::max(amax, std::max(df_tau, 0.0) / U[tau_i][i]);
	}
	return amax;
}

double hydro_grid::compute_flux(int rk) {
	amax = 0.0;
	const double a = cosmos_a();
	std::vector<multi_array<double>> VR(opts.nhydro);
	std::vector<multi_array<double>> VL(opts.nhydro);
	const auto urlbox = box.pad(2 - opts.hbw);
	const auto grad_box = box.pad(-opts.hbw + 1);
	std::vector<multi_array<double>> V(opts.nhydro);
	for (multi_iterator i(box); !i.end(); i++) {
		U[rho_i][i] /= std::pow(a, NDIM);
		if (opts.species) {
			for (int si = 0; si < opts.nspecies; si++) {
				U[spc_i + si][i] /= std::pow(a, NDIM);
			}
		}
		U[tau_i][i] /= std::pow(a, NDIM);
		for (int dim = 0; dim < NDIM; dim++) {
			U[sx_i + dim][i] /= std::pow(a, NDIM + 1);
		}
		U[egas_i][i] /= std::pow(a, (NDIM + 2));
		U[pot_i][i] = U[rho_i][i] * phi[i];
	}
	for (int i = 0; i < opts.nhydro; i++) {
		VR[i].resize(urlbox);
		VL[i].resize(urlbox);
		V[i].resize(box);
	}
	for (multi_iterator i(box); !i.end(); i++) {
		if (U[rho_i][i] < 0.0) {
			printf("compute_flux: rho is less than zero %e\n", U[rho_i][i]);
			abort();
		}
		if (U[tau_i][i] < 0.0) {
			printf("compute_flux: tau is less than zero %e\n", U[tau_i][i]);
			abort();
		}
		assert(U[rho_i][i] > 0.0);
		assert(U[tau_i][i] >= 0.0);
		const auto rhoinv = 1.0 / U[rho_i][i];
		double ek = 0.0;
		V[rho_i][i] = U[rho_i][i];
		if (opts.species) {
			for (int si = 0; si < opts.nspecies; si++) {
				V[spc_i + si][i] = U[spc_i + si][i];
			}
		}
		for (int dim = 0; dim < NDIM; dim++) {
			ek += 0.5 * std::pow(U[sx_i + dim][i], 2) * rhoinv;
			V[sx_i + dim][i] = U[sx_i + dim][i] * rhoinv;
		}
		V[egas_i][i] = U[egas_i][i] - ek;
		V[tau_i][i] = std::pow(U[tau_i][i], opts.gamma);
		V[pot_i][i] = phi[i];
	}
	multi_array<double> drho_dt(box);
	for (multi_iterator i(box); !i.end(); i++) {
		drho_dt[i] = 0.0;
	}
	for (int dim = 0; dim < NDIM; dim++) {
		for (int i = 0; i < opts.nhydro; i++) {
			auto this_box = fbox[dim];
			this_box.min[dim]--;
			for (multi_iterator j(this_box); !j.end(); j++) {
				const auto du = 0.5 * V[i].minmod_gradient(dim, j);
				VR[i][j] = V[i][j] - du;
				auto jp1 = j.index();
				jp1[dim]++;
				VL[i][jp1] = V[i][j] + du;
			}
		}
		std::vector<double> vr(opts.nhydro), vl(opts.nhydro), ur(opts.nhydro), ul(opts.nhydro), flux(opts.nhydro);
		for (multi_iterator j(fbox[dim]); !j.end(); j++) {
			for (int f = 0; f < opts.nhydro; f++) {
				vr[f] = VR[f][j];
				vl[f] = VL[f][j];
			}
			ur[rho_i] = vr[rho_i];
			ul[rho_i] = vl[rho_i];
			if (opts.species) {
				for (int si = 0; si < opts.nspecies; si++) {
					ur[spc_i + si] = vr[spc_i + si];
					ul[spc_i + si] = vl[spc_i + si];
				}
			}
			ur[pot_i] = vr[rho_i] * vr[pot_i];
			ul[pot_i] = vl[rho_i] * vl[pot_i];
			double ekr = 0.0;
			double ekl = 0.0;
			for (int dim = 0; dim < NDIM; dim++) {
				ekr += 0.5 * vr[rho_i] * std::pow(vr[sx_i + dim], 2);
				ekl += 0.5 * vl[rho_i] * std::pow(vl[sx_i + dim], 2);
				ur[sx_i + dim] = vr[sx_i + dim] * vr[rho_i];
				ul[sx_i + dim] = vl[sx_i + dim] * vl[rho_i];
			}
			ur[egas_i] = vr[egas_i] + ekr;
			ul[egas_i] = vl[egas_i] + ekl;
			ur[tau_i] = std::pow(vr[tau_i], 1.0 / opts.gamma);
			ul[tau_i] = std::pow(vl[tau_i], 1.0 / opts.gamma);
			const auto this_amax = hydro_flux(flux, ul, ur, dim);
			amax = std::max(amax, this_amax);
			auto jm = j.index();
			jm[dim]--;
			drho_dt[j] += flux[rho_i] / dx;
			drho_dt[jm] -= flux[rho_i] / dx;
			flux[egas_i] += flux[pot_i];
			flux[pot_i] = 0.0;
			flux[rho_i] *= std::pow(a, NDIM - 1);
			flux[tau_i] *= std::pow(a, NDIM - 1);
			if (opts.species) {
				for (int si = 0; si < opts.nspecies; si++) {
					flux[spc_i + si] *= std::pow(a, NDIM - 1);
				}
			}
			for (int dim2 = 0; dim2 < NDIM; dim2++) {
				flux[sx_i + dim2] *= std::pow(a, NDIM);
			}
			flux[egas_i] *= std::pow(a, (NDIM + 2) - 1);
			for (int f = 0; f < opts.nhydro; f++) {
				F[dim][f][j] = (1.0 - opts.beta[rk]) * F[dim][f][j] + opts.beta[rk] * flux[f];
			}
		}
	}
	if (opts.gravity) {
		for (multi_iterator j(box.pad(-opts.hbw)); !j.end(); j++) {
			vect<double> g;
			//		double sdotg = 0.0;
			for (int dim = 0; dim < NDIM; dim++) {
				auto jp = j.index();
				auto jm = j.index();
				jp[dim]++;
				jm[dim]--;
				g[dim] = 0.5 * (-phi[jp] + phi[jm]) / dx / a;
				//			sdotg += U[sx_i + dim][j] * g[dim];
			}
			std::vector<double> this_s(opts.nhydro, 0.0);
			for (int dim = 0; dim < NDIM; dim++) {
				this_s[sx_i + dim] += std::pow(a, NDIM + 1) * U[rho_i][j] * g[dim];
			}
//			this_s[egas_i] += std::pow(a, (NDIM + 2)) * sdotg;
			this_s[egas_i] -= std::pow(a, (NDIM + 2) - 1) * phi[j] * drho_dt[j];


			for (int f = 0; f < opts.nhydro; f++) {
				S[f][j] = (1.0 - opts.beta[rk]) * S[f][j] + opts.beta[rk] * this_s[f];
			}
		}
	}
	for (multi_iterator i(box); !i.end(); i++) {
		U[rho_i][i] *= std::pow(a, NDIM);
		U[tau_i][i] *= std::pow(a, NDIM);
		if (opts.species) {
			for (int si = 0; si < opts.nspecies; si++) {
				U[spc_i + si][i] *= std::pow(a, NDIM);
			}
		}
		for (int dim = 0; dim < NDIM; dim++) {
			U[sx_i + dim][i] *= std::pow(a, NDIM + 1);
		}
		U[egas_i][i] *= std::pow(a, (NDIM + 2));
	}
	return amax;
}

void hydro_grid::substep_update(int rk, double dt, double a0, double a1) {
	const double a = cosmos_a();
	const double onembeta = 1.0 - opts.beta[rk];
	const double lambda = dt / dx;
	for (multi_iterator i(box.pad(-opts.hbw)); !i.end(); i++) {
		for (int f = 0; f < opts.nhydro; f++) {
			auto &u = U[f][i];
			u = U0[f][i];
			for (int dim = 0; dim < NDIM; dim++) {
				multi_index ip1 = i;
				ip1[dim]++;
				u -= (F[dim][f][ip1] - F[dim][f][i]) * lambda;
			}
			u += S[f][i] * dt;
		}
	}
	for (multi_iterator i(box.pad(-opts.hbw)); !i.end(); i++) {
		if (U[rho_i][i] <= 0.0) {
			printf("substep_update: rho is less than zero %e\n", U[rho_i][i]);
			abort();
		}
		if (U[tau_i][i] < 0.0) {
			double ekin = 0.0;
			for (int dim = 0; dim < NDIM; dim++) {
				ekin += 0.5 * std::pow(U[sx_i + dim][i], 2) * std::pow(a, -NDIM - 2) / U[rho_i][i];
			}
			const double eint = U[egas_i][i] * std::pow(a, -(NDIM + 2)) - ekin;

			printf("substep_update %i: tau is less than zero %e %e\n", rk, U[tau_i][i], eint);
			if (eint < 0.0) {
				abort();
			} else {
				U[tau_i][i] = std::pow(a, NDIM) * std::pow(eint, 1.0 / opts.gamma);
			}
		}
	}
	if (rk == opts.nrk - 1) {
		const auto factor = dt / (1 << (NDIM - 1));
		for (int dim = 0; dim < NDIM; dim++) {
			for (int f = 0; f < opts.nhydro; f++) {
				for (multi_iterator i(fbox[dim]); !i.end(); i++) {
					if (i.index()[dim] % 2 == 0) {
						const auto j = i.index() / 2;
						Fc[dim][f][j] += F[dim][f][i] * factor;
					}
				}
			}
		}
	}
}

void hydro_grid::resize(double dx_, multi_range box_) {
	dx = dx_;
	box = box_;
	const int sz = opts.nhydro * box.volume();
	double *u0 = new double[sz];
	double *u1 = new double[sz];
	int j = 0;
	U0.resize(opts.nhydro);
	U.resize(opts.nhydro);
	S.resize(opts.nhydro);
	phi.resize(box);
	error.resize(box_.pad(1));
	for (int f = 0; f < opts.nhydro; f++) {
		U0[f].resize(box);
		U[f].resize(box);
		S[f].resize(box.pad(-opts.hbw));
	}
	for (int dim = 0; dim < NDIM; dim++) {
		F[dim].resize(opts.nhydro);
		Fc[dim].resize(opts.nhydro);
		for (int i = 0; i < opts.nhydro; i++) {
			fbox[dim] = box.pad(-opts.hbw);
			fcbox[dim] = fbox[dim].half();
			fbox[dim].max[dim]++;
			fcbox[dim].min[dim] = fbox[dim].min[dim] / 2;
			fcbox[dim].max[dim] = (fbox[dim].max[dim] - 1) / 2 + 1;
			F[dim][i].resize(fbox[dim]);
			Fc[dim][i].resize(fcbox[dim]);
		}
	}
	R.resize(box_.pad(opts.window));
	reset_coarse_flux_registers();
	if (opts.problem == "rt") {
		for (multi_iterator i(box); !i.end(); i++) {
			auto y = coord(i.index()[1]);
			if (y < 0.0) {
				y = -y;
			} else if (y > 1.0) {
				y = 2.0 - y;
			}
			phi[i] = 0.1 * y;
		}
	}
}

void hydro_grid::reset_coarse_flux_registers() {
	for (int dim = 0; dim < NDIM; dim++) {
		for (int f = 0; f < opts.nhydro; f++) {
			for (multi_iterator i(fcbox[dim]); !i.end(); i++) {
				Fc[dim][f][i] = 0.0;
			}
		}
	}
}

void hydro_grid::reset_flux_registers() {
	for (int dim = 0; dim < NDIM; dim++) {
		for (int f = 0; f < opts.nhydro; f++) {
			for (multi_iterator i(fbox[dim]); !i.end(); i++) {
				F[dim][f][i] = 0.0;
			}
		}
	}
	for (int f = 0; f < opts.nhydro; f++) {
		for (multi_iterator i(box.pad(-opts.hbw)); !i.end(); i++) {
			S[f][i] = 0.0;
		}
	}
}

hydro_grid::~hydro_grid() {
}

void hydro_grid::compute_refinement_criteria(const multi_array<int> &pcount) {
	const auto ibox = box.pad(-opts.hbw + opts.window);
	for (multi_iterator i(ibox); !i.end(); i++) {
		R[i] = false;
	}
	if (opts.particles) {
		for (multi_iterator i(box.pad(-opts.hbw)); !i.end(); i++) {
			R[i] = pcount[i] > 7;
		}
	} else if (opts.hydro) {
		multi_range dirs(multi_range(index_type(0)).pad(1));
		int dirmax = (std::pow(3, NDIM) - 1) / 2;
		int cnt = 0;
		for (multi_iterator i(dirs); !i.end(); i++) {
			if (cnt >= dirmax) {
				break;
			}
			cnt++;
			bool all_zero = true;
			const auto j = i.index();
			int n_one = 0;
			for (int dim = 0; dim < NDIM; dim++) {
				if (j[dim] != 0) {
					n_one++;
				}
			}
			if (n_one == 0) {
				continue;
			}
			for (multi_iterator k(box.pad(-opts.hbw)); !k.end(); k++) {
				if (!R[k] && U[rho_i][k] > 1.0e-6) {
					constexpr double tau = 0.1;
					for (int f = 0; f < opts.nhydro; f++) {
						const auto up = U[f][k.index() + j];
						const auto u0 = U[f][k];
						const auto um = U[f][k.index() - j];
						double c0;
						if (f < sx_i || f >= sx_i + NDIM) {
							c0 = tau * (std::abs(up) + std::abs(um) + 2.0 * std::abs(u0));
						} else {
							double vp = 0.0;
							double v0 = 0.0;
							double vm = 0.0;
							for (int dim = 0; dim < NDIM; dim++) {
								vm += std::max(0.0, U[egas_i][k.index() - j] / U[rho_i][k.index() - 1]);
								v0 += std::max(0.0, U[egas_i][k] / U[rho_i][k]);
								vp += std::max(0.0, U[egas_i][k.index() + j] / U[rho_i][k.index() + 1]);
							}
							c0 = tau * (std::sqrt(vp) + std::sqrt(vm) + 2.0 * std::sqrt(v0));
						}
						const auto num = std::abs(up + um - 2.0 * u0) / std::sqrt(n_one);
						const auto den = std::abs(up - u0) + std::abs(u0 - um) + c0;
						if (den > 0.0) {
							if (num / den > opts.refine_slope) {
								R[k] = true;
							}
						}
//					R[k] = true;
					}
				}
			}
		}
	}
}

void hydro_grid::initialize() {
	const auto a = cosmos_a();
	if (opts.problem == "cosmos") {
		const auto &parts = fileio_get_particles();
		for (multi_iterator i(box); !i.end(); i++) {
			U[rho_i][i] = 0.0;
			U[egas_i][i] = 0.0;
			for (int dim = 0; dim < NDIM; dim++) {
				U[sx_i + dim][i] = 0.0;
			}
		}
		multi_range pbox = multi_range(vect<int>(0)).pad(1);
		for (multi_iterator q(pbox); !q.end(); q++) {
			for (const auto &p : parts) {
				const auto x = p.x + vect<double>(q.index());
				multi_index j;
				for (int dim = 0; dim < NDIM; dim++) {
					j[dim] = int((x[dim] + 1.0 - 0.5 * dx) / dx) - int(1.0 / dx);
				}
				auto this_box = multi_range(j);
				for (int dim = 0; dim < NDIM; dim++) {
					this_box.max[dim] += 2;
					this_box.min[dim] -= 2;
				}
				if (box.intersection(this_box).volume()) {
					for (multi_iterator k(this_box); !k.end(); k++) {
						if (box.contains(k)) {
							double r = 0.0;
							for (int dim = 0; dim < NDIM; dim++) {
								r += std::pow(x[dim] - coord(k[dim]), 2);
							}
							r = std::sqrt(r);
							r = r / (2.0 * dx);
							double wt;
							if (r < 0.5) {
								wt = (1.0 - 6.0 * r * r * (1 - r));
							} else if (r < 1.0) {
								wt = 2.0 * std::pow(1.0 - r, 3.0);
							} else {
								wt = 0.0;
							}
							wt *= 8.0 / M_PI;
							const auto drho = wt * p.m * opts.omega_b / opts.omega_m / (std::pow(2.0 * a * dx, NDIM));
							U[rho_i][k] += drho;
							U[egas_i][k] += 0.5 * p.v.dot(p.v) * drho / (a * a);
							for (int dim = 0; dim < NDIM; dim++) {
								U[sx_i + dim][k] += p.v[dim] * drho / a;
							}
						}
					}
				}
			}
		}
		for (multi_iterator i(box); !i.end(); i++) {
			double ek = 0.0;
			for (int dim = 0; dim < NDIM; dim++) {
				ek += 0.5 * U[sx_i + dim][i] * U[sx_i + dim][i] / U[rho_i][i];
			}
			U[tau_i][i] = std::pow(U[egas_i][i] - ek, 1.0 / opts.gamma);
			if (opts.species) {
				const auto nH = 0.75 * U[rho_i][i];
				const auto nHe = 0.25 * U[rho_i][i];
				U[h_i][i] = 0.0 * nH;
				U[hp_i][i] = 1.0 * nH;
				U[he_i][i] = 0.0 * nHe;
				U[hep_i][i] = 0.0 * nHe;
				U[hepp_i][i] = 1.0 * nHe;
			}
		}
	} else {
		for (multi_iterator i(box); !i.end(); i++) {
			double r = 0.0;
			for (int dim = 0; dim < NDIM; dim++) {
				r += std::pow(coord(i[dim]) - 0.5, 2);
			}
			r = std::sqrt(r);
			for (int dim = 0; dim < NDIM; dim++) {
				U[sx_i + dim][i] = 0.0;
			}
			if (opts.problem == "sod") {
				double xsum = 0.0;
				for (int dim = 0; dim < 1; dim++) {
					xsum += coord(i[dim]) - 0.5;
				}
				U[sx_i][i] = 0.0;
				if (xsum > 0.0) {
					U[rho_i][i] = 1.0;
					U[egas_i][i] = 1.0;
				} else {
					U[rho_i][i] = 0.1;
					U[egas_i][i] = 0.125;
				}
				double ek = 0.0;
				for (int dim = 0.0; dim < NDIM; dim++) {
					ek += std::pow(U[sx_i + dim][i], 2) / U[rho_i][i];
				}
				U[tau_i][i] = std::pow(U[egas_i][i], 1.0 / opts.gamma);
				U[egas_i][i] -= ek;
			} else if (opts.problem == "blast") {
				U[rho_i][i] = 1.0;
				if (r < 3.0 * dx) {
					U[egas_i][i] = 1.0;
				} else {
					U[egas_i][i] = 1.0e-3;
				}
				U[tau_i][i] = std::pow(U[egas_i][i], 1.0 / opts.gamma);
			} else if (opts.problem == "kh") {
				const int dim = NDIM - 1;
				const double rnd = 2.0 * rand1() - 1.0;
				if (std::abs(coord(i.index()[dim]) - 0.5) > 0.25) {
					U[rho_i][i] = 1.0;
					U[sx_i][i] = +0.5 + rnd * 0.0001;
				} else {
					U[sx_i][i] = -0.5;
					U[rho_i][i] = 2.0 + 2.0 * rnd * 0.0001;
				}
				for (int dim = 1; dim < NDIM; dim++) {
					const double rnd = 2.0 * (double) (rand() + 0.5) / RAND_MAX - 1.0;
					U[sx_i + dim][i] = rnd * 0.001 * U[rho_i][i];
				}
				const auto p = 2.5;
				const auto ein = p / (opts.gamma - 1.0);
				double ek = 0.0;
				for (int dim = 0; dim < NDIM; dim++) {
					ek += 0.5 * std::pow(U[sx_i + dim][i], 2) / U[rho_i][i];
				}
				U[egas_i][i] = ein + ek;
				U[tau_i][i] = std::pow(ein, 1.0 / opts.gamma);
			} else if (opts.problem == "sphere") {
				for (int f = 0; f < opts.nhydro; f++) {
					U[f][i] = 0.0;
				}
				constexpr auto r0 = 0.1;
				constexpr auto rho0 = 1.0 / std::pow(r0, NDIM);
				if (r < 0.5 * r0) {
					U[rho_i][i] = rho0 * (1.0 - 6.0 * r / r0 * r / r0 * (1.0 - r / r0));
				} else if (r < r0) {
					U[rho_i][i] = rho0 * 2.0 * std::pow(1.0 - r / r0, 3);
				} else {
					U[rho_i][i] = 1.0e-6 / 4.0 / M_PI;
				}
				if ( NDIM == 3) {
					U[rho_i][i] *= 8.0 / M_PI;
				} else if ( NDIM == 2) {
					U[rho_i][i] *= 40.0 / (7.0 * M_PI);
				} else {
					U[rho_i][i] *= 4.0 / 3.0;
				}
				U[egas_i][i] = 1e-5;
			} else if (opts.problem == "rt") {
				const auto y = coord(i[1]);
				double p;
				if (y > 0.5) {
					U[rho_i][i] = 2.0 * (1.0 + rand1() * 0.000);
				} else {
					U[rho_i][i] = 1.0 * (1.0 + rand1() * 0.000);
				}
				if (std::abs(y - 0.5) < 0.25) {
					U[sy_i][i] = rand1() * 0.005 * (1.0 + std::cos(4.0 * M_PI * (y - 0.5))) * U[rho_i][i];
				}
				p = std::max(2.5 - 0.1 * y * U[rho_i][i], 1.0e-3);
				const auto ein = p / (opts.gamma - 1.0);
				double ek = 0.0;
				for (int dim = 0; dim < NDIM; dim++) {
					ek += 0.5 * std::pow(U[sx_i + dim][i], 2) / U[rho_i][i];
				}
				U[egas_i][i] = ein + ek;
				U[tau_i][i] = std::pow(ein, 1.0 / opts.gamma);
			} else if (opts.problem == "polytrope") {
				double r = 0.0;
				for (int dim = 0; dim < NDIM; dim++) {
					auto x = coord(i[dim]);
					x -= 0.5;
					r += std::pow(x, 2);
				}
				r = std::sqrt(r);
				const auto n = 1.5;
				const auto alpha = 0.02;
				const auto theta = lane_emden(r / alpha, dx / alpha / 2.0, n);
				assert(theta <= 1.0);
				auto &rho = U[rho_i][i];
				rho = std::max(1.0e-10, std::pow(theta, n));
				const auto vx = 1.0;
				const auto vy = 0.00;
				U[sx_i][i] = vx * rho;
				U[sy_i][i] = vy * rho;
				const auto K = 4.0 * M_PI * alpha * alpha * cosmos_a() * cosmos_a() / (n + 1);
				U[egas_i][i] = std::max(1.0e-12, K * std::pow(theta, 1.0 + n) / (opts.gamma - 1.0)); // * std::pow(unit * unit, opts.gamma);
				U[tau_i][i] = std::pow(U[egas_i][i], 1.0 / opts.gamma);
				U[egas_i][i] += 0.5 * (vx * vx + vy * vy) * rho;
			} else {
				printf("unknown problem %s\n", opts.problem.c_str());
				abort();
			}

		}
	}
	for (multi_iterator i(box); !i.end(); i++) {
		U[rho_i][i] *= std::pow(a, NDIM);
		if (opts.species) {
			for (int si = 0; si < opts.nspecies; si++) {
				U[spc_i + si][i] *= std::pow(a, NDIM);
			}
		}
		U[tau_i][i] *= std::pow(a, NDIM);
		for (int dim = 0; dim < NDIM; dim++) {
			U[sx_i + dim][i] *= std::pow(a, NDIM + 1);
		}
		U[egas_i][i] *= std::pow(a, (NDIM + 2));
	}
}

void hydro_grid::enforce_physical_bc(int dir) {
	const auto ibox = box.pad(-opts.hbw);
	multi_range bbox = box;
	const int dim = dir / 2;
	if (dir % 2 == 0) {
		bbox.min[dim] = ibox.min[dim] - opts.hbw;
		bbox.max[dim] = ibox.min[dim];
	} else if (dir % 2 == 1) {
		bbox.max[dim] = ibox.max[dim] + opts.hbw;
		bbox.min[dim] = ibox.max[dim];
	}

	for (multi_iterator i(bbox); !i.end(); i++) {
		int xi;
		if (dir % 2 == 0 && !opts.bnd[dir] == PERIODIC) {
			if (opts.bnd[dir] == OUTFLOW) {
				xi = ibox.min[dim];
			} else if (opts.bnd[dir] == REFLECTING) {
				xi = -i.index()[dim] + opts.hbw + 1;
			}
		} else {
			if (opts.bnd[dir] == OUTFLOW) {
				xi = ibox.max[dim] - 1;
			} else if (opts.bnd[dir] == REFLECTING) {
				xi = -i.index()[dim] + 2 * ibox.max[dim] - 2 * opts.hbw - 1;
			}
		}
		if (opts.bnd[dir] != PERIODIC) {
			auto j = i.index();
			j[dim] = xi;
			for (int f = 0; f < opts.nhydro; f++) {
				U[f][i] = U[f][j];
			}
			const int si = sx_i + dim;
			if (opts.bnd[dir] == OUTFLOW) {
				U[si][i] = 0.0;
			} else if (opts.bnd[dir] == REFLECTING) {
				U[si][i] = -U[si][i];
			}
		}
	}

}

void hydro_grid::update_energy() {
	const auto a = cosmos_a();
	for (multi_iterator i(box.pad(-opts.hbw)); !i.end(); i++) {
		double max_egas = 0.0;
		double ekin = 0.0;
		for (multi_iterator d(multi_range(0).pad(1)); !d.end(); d++) {
			max_egas = std::max(max_egas, U[egas_i][i]);
		}
		max_egas *= std::pow(a, -(NDIM + 2));
		for (int dim = 0; dim < NDIM; dim++) {
			ekin += 0.5 * std::pow(U[sx_i + dim][i], 2) * std::pow(a, -NDIM - 2) / U[rho_i][i];
		}
		const double eint = U[egas_i][i] * std::pow(a, -(NDIM + 2)) - ekin;
		if (eint > max_egas * 0.1) {
			U[tau_i][i] = std::pow(a, NDIM) * std::pow(eint, 1.0 / opts.gamma);
		}
	}
}

double hydro_grid::coord(index_type i) const {
	return (i + 0.5) * dx;
}

std::vector<multi_range> hydro_grid::refined_ranges(const std::vector<multi_range> &amr_boxes, const std::vector<multi_range> &forced) const {

	const auto ibox = box.pad(-opts.hbw);
	auto Rwin = R;
	for (multi_iterator i(ibox.pad(opts.window)); !i.end(); i++) {
		bool res = R[i];
		if (!res) {
			for (const auto &f : forced) {
				if (f.contains(i)) {
					res = true;
					break;
				}
			}
		}
		if (res) {
			Rwin[i] = true;
			const auto window = range<index_type>(i).pad(opts.window).intersection(box.pad(-opts.hbw));
			for (multi_iterator j(window); !j.end(); j++) {
				Rwin[j] = true;
			}
		}
	}

	std::vector<multi_range> boxes;
	std::vector<multi_range> tmp;
	std::vector<multi_range> finished;
	multi_range range;
	range.min = box.max;
	range.max = box.min;
	for (multi_iterator i(ibox); !i.end(); i++) {
		if (Rwin[i]) {
			for (int dim = 0; dim < NDIM; dim++) {
				range.min[dim] = std::min(range.min[dim], i[dim]);
				range.max[dim] = std::max(range.max[dim], i[dim] + 1);
			}
		}
	}
	if (range.volume()) {
		boxes.push_back(range);
		for (const auto &amr : amr_boxes) {
			tmp.resize(0);
			for (const auto &b : boxes) {
				const auto amrbox = amr.pad(2);
				if (amrbox.intersection(b).volume()) {
					auto tmp2 = b.subtract(amrbox);
					tmp.insert(tmp.end(), tmp2.begin(), tmp2.end());
				} else {
					tmp.push_back(b);
				}
			}
			boxes = std::move(tmp);
		}
		tmp.resize(0);
		for (const auto &b : boxes) {
			bool found = false;
			for (multi_iterator i(b); !i.end(); i++) {
				if (Rwin[i]) {
					found = true;
					break;
				}
			}
			if (found) {
				tmp.push_back(b);
			}
		}
		boxes = std::move(tmp);
		bool change;
		const auto max_vol = std::pow(opts.max_box / 2, NDIM);
		const auto min_vol = std::pow(opts.min_box / 2, NDIM);
		do {
			tmp.resize(0);
			for (const auto &b : boxes) {
				assert(b.volume());
				int count = 0;
				for (multi_iterator i(b); !i.end(); i++) {
					if (Rwin[i]) {
						count++;
					}
				}
				const auto this_vol = b.volume();
				const double efficiency = (double) count / (double) this_vol;

				if (this_vol > max_vol || /* (*/efficiency < opts.efficiency/*&& this_vol >= min_vol)*/) {
					auto new_boxes = b.split();
					assert(b.volume());
					assert(new_boxes.first.volume());
					assert(new_boxes.second.volume());
					range.min = new_boxes.first.max;
					range.max = new_boxes.first.min;
					for (multi_iterator i(new_boxes.first); !i.end(); i++) {
						if (Rwin[i]) {
							for (int dim = 0; dim < NDIM; dim++) {
								range.min[dim] = std::min(range.min[dim], i[dim]);
								range.max[dim] = std::max(range.max[dim], i[dim] + 1);
							}
						}
					}
					new_boxes.first = range;
					range.min = new_boxes.second.max;
					range.max = new_boxes.second.min;
					for (multi_iterator i(new_boxes.second); !i.end(); i++) {
						if (Rwin[i]) {
							for (int dim = 0; dim < NDIM; dim++) {
								range.min[dim] = std::min(range.min[dim], i[dim]);
								range.max[dim] = std::max(range.max[dim], i[dim] + 1);
							}
						}
					}
					new_boxes.second = range;
					if (new_boxes.first.volume()) {
						tmp.push_back(new_boxes.first);
					}
					if (new_boxes.second.volume()) {
						tmp.push_back(new_boxes.second);

					}
					assert(new_boxes.second.volume() || new_boxes.first.volume());
				} else {
					assert(b.volume());
					finished.push_back(b);
				}
			}
			change = tmp.size();
			boxes = std::move(tmp);
		} while (change);
		boxes = std::move(finished);
	}
	return std::move(boxes);
}

std::vector<double> hydro_grid::pack(multi_range bbox) const {
	std::vector<double> data;
	for (int f = 0; f < opts.nhydro; f++) {
		for (multi_iterator i(bbox); !i.end(); i++) {
			data.push_back(U[f][i]);
		}
	}
	return data;
}

std::vector<std::uint8_t> hydro_grid::pack_refinement(multi_range bbox) const {
//	printf("Packing %s\n", bbox.to_string().c_str());
	std::vector<std::uint8_t> data;
	for (multi_iterator i(bbox); !i.end(); i++) {
		data.push_back(R[i]);
	}
	return data;
}

std::vector<double> hydro_grid::pack_field(int f, multi_range bbox) const {
	std::vector<double> data;
	for (multi_iterator i(bbox); !i.end(); i++) {
		data.push_back(U[f][i]);
	}
	return data;
}

void hydro_grid::set_error_field(multi_array<double> &&field) {
	error = std::move(field);
}

double hydro_grid::compare_analytic(const std::vector<multi_range> &cboxes, multi_array<double> &results) const {
	const auto r0 = 0.1;
	const auto force = [](double r) {
		const auto r2 = r * r;
		const auto r3 = r2 * r;
		const auto r4 = r2 * r2;
		if ( NDIM == 3) {
			if (r < 0.5) {
				return (-32.0 / 3.0) * r + (192.0 / 5.0) * r3 - 32.0 * r4;
			} else if (r < 1.0) {
				const auto rinv = 1.0 / r;
				const auto r2inv = rinv * rinv;
				return (1.0 / 15.0) * r2inv - (64.0 / 3.0) * r + 48.0 * r2 - (192.0 / 5.0) * r3 + (32.0 / 3.0) * r4;
			} else {
				return -1.0 / (r * r);
			}
		} else if ( NDIM == 2) {
			if (r < 0.5) {
				return 2.0 * ((-40.0 / 7.0) * r + (120.0 / 7.0) * r3 - (96.0 / 7.0) * r4);
			} else if (r < 1.0) {
				const auto rinv = 1.0 / r;
				return 2.0 * ((1.0 / 7.0) * rinv - (80.0 / 7.0) * r + (160.0 / 7.0) * r2 - (120.0 / 7.0) * r3 + (32.0 / 7.0) * r4);
			} else {
				return -2.0 / r;
			}
		} else {
			if (r < 0.5) {
				return M_PI * ((-32.0 / 3.0) * r + (64.0 / 3.0) * r3 - (16.0) * r4);
			} else if (r < 1.0) {
				return M_PI * ((4.0 / 3.0) - (64.0 / 3.0) * r + (32.0) * r2 - (64.0 / 3.0) * r3 + (16.0 / 3.0) * r4);
			} else {
				return -4.0 * M_PI;
			}
		}
	};
	std::array<multi_index, NDIM> dir;
	for (int dim = 0; dim < NDIM; dim++) {
		dir[dim] = 0;
		dir[dim][dim] = 1;
	}
	double err_rms = 0.0;
	results.resize(box.pad(-opts.hbw));
	for (multi_iterator i(box.pad(-opts.hbw)); !i.end(); i++) {
		bool child = false;
		for (const auto &this_box : cboxes) {
			if (this_box.contains(i)) {
				child = true;
				break;
			}
		}
		if (!child) {
			std::array<double, NDIM> g;
			std::array<double, NDIM> ga;
			std::array<double, NDIM> dg;
			for (int dim = 0; dim < NDIM; dim++) {
				g[dim] = -(phi[i.index() + dir[dim]] - phi[i.index() - dir[dim]]) / (2.0 * dx);
			}
			double r = 0.0;
			for (int dim = 0; dim < NDIM; dim++) {
				r += std::pow(coord(i.index()[dim]) - 0.5, 2);
			}
			r = std::sqrt(r) / r0;
			const auto f = force(r) / std::pow(r0, NDIM - 1);
			for (int dim = 0; dim < NDIM; dim++) {
				ga[dim] = f * (coord(i.index()[dim]) - 0.5) / (r * r0);
				dg[dim] = g[dim] - ga[dim];
			}
			double dgabs = 0.0, gabs = 0.0;
			for (int dim = 0; dim < NDIM; dim++) {
				dgabs += dg[dim] * dg[dim];
			}
			for (int dim = 0; dim < NDIM; dim++) {
				gabs += ga[dim] * ga[dim];
			}
			if (gabs != 0.0) {
				dgabs = std::sqrt(dgabs);
				gabs = std::sqrt(gabs);
				const auto err = dgabs / gabs;
				results[i] = err;
				err_rms += err * err;
			} else {
				results[i] = 1e+10;
			}
		}
	}
	err_rms *= std::pow(dx, NDIM);
	return err_rms;
}

std::vector<std::vector<double>> hydro_grid::pack_output(const multi_array<std::uint8_t> &mask) const {
	std::vector<std::vector<double>> data;
	auto bbox = box.pad(-opts.hbw);
//	auto bbox = box;
	data.resize(opts.nhydro + 2);
	double units;
	const auto a = cosmos_a();
	for (int f = 0; f <= opts.nhydro + 1; f++) {
		if (f == rho_i || (f >= spc_i && f < spc_i + opts.nspecies)) {
			units = opts.code_to_g / std::pow(opts.code_to_cm * a, NDIM);
		} else if (f == tau_i) {
			units = std::pow(opts.code_to_g / std::pow(opts.code_to_cm, NDIM - 2) / std::pow(opts.code_to_s, 2), 1.0 / opts.gamma) / std::pow(a, NDIM);
		} else if (f >= sx_i && f <= sz_i) {
			units = opts.code_to_g / std::pow(opts.code_to_cm, NDIM - 1) / opts.code_to_s / std::pow(a, NDIM + 1);
		} else if (f == egas_i) {
			units = opts.code_to_g / std::pow(opts.code_to_cm, NDIM - 2) / std::pow(opts.code_to_s, 2) / std::pow(a, NDIM + 2);
		} else if (f == pot_i) {
			units = opts.code_to_g / std::pow(opts.code_to_cm, NDIM - 2) / std::pow(opts.code_to_s, 2);
		} else {
			units = 1.0;
		}
		const auto &u = f == opts.nhydro ? phi : (f == opts.nhydro + 1 ? error : U[f]);
		for (multi_iterator i(bbox); !i.end(); i++) {
			if (mask[i]) {
				data[f].push_back(u[i] * units);
			}
		}
	}
	return data;
}

std::vector<double> hydro_grid::pack_prolong(multi_range bbox, double w) const {
	std::vector<double> data;
	for (int f = 0; f < opts.nhydro; f++) {
		auto pro0 = U0[f].prolong(bbox);
		auto pro1 = U[f].prolong(bbox);
		for (multi_iterator i(bbox); !i.end(); i++) {
			data.push_back((1.0 - w) * pro0[i] + w * pro1[i]);
		}
	}
	return data;
}

std::vector<double> hydro_grid::pack_field_prolong(int f, multi_range bbox, double w) const {
	std::vector<double> data;
	if (w == 0.0 || w == 1.0) {
		auto pro = w == 0.0 ? U0[f].prolong(bbox) : U[f].prolong(bbox);
		for (multi_iterator i(bbox); !i.end(); i++) {
			data.push_back(pro[i]);
		}
	} else {
		auto pro0 = U0[f].prolong(bbox);
		auto pro1 = U[f].prolong(bbox);
		for (multi_iterator i(bbox); !i.end(); i++) {
			data.push_back((1.0 - w) * pro0[i] + w * pro1[i]);
		}
	}
	return data;
}

std::vector<double> hydro_grid::pack_restrict(multi_range bbox) const {
//	printf( "pack restrict %s\n", bbox.to_string().c_str());
	std::vector<double> data;
	for (int f = 0; f < opts.nhydro; f++) {
		auto res = U[f].restrict_(bbox);
		for (multi_iterator i(bbox); !i.end(); i++) {
			data.push_back(res[i]);
		}
	}
	return data;
}

std::vector<double> hydro_grid::pack_coarse_flux() {
	const auto inv = 1.0 / (1 << (NDIM - 1));
	std::vector<double> data;
	int size = 0;
	for (int dim = 0; dim < NDIM; dim++) {
		auto cbox = box.pad(-opts.hbw).half();
		cbox.max[dim]++;
		size += cbox.volume();
	}
	data.reserve(size);
	for (int dim = 0; dim < NDIM; dim++) {
		auto cbox = box.pad(-opts.hbw).half();
		cbox.max[dim]++;
		for (int f = 0; f < opts.nhydro; f++) {
			for (multi_iterator i(cbox); !i.end(); i++) {
				double this_flux = 0.0;
				auto this_box = multi_range(i.index()).double_();
				this_box.max[dim]--;
				for (multi_iterator j(this_box); !j.end(); j++) {
					this_flux += F[dim][f][j] * inv;
				}
				data.push_back(this_flux);
			}
		}
	}
	data.push_back(amax);
	return data;
}

double hydro_grid::unpack_fine_flux(const std::vector<double> &data, const multi_range &bbox) {
	int k = 0;
	double a;
	for (int dim = 0; dim < NDIM; dim++) {
		auto this_box = bbox;
		this_box.max[dim]++;
		for (int f = 0; f < opts.nhydro; f++) {
			for (multi_iterator i(this_box); !i.end(); i++) {
				if (fbox[dim].contains(i)) {
					assert(k < data.size());
					if (f == rho_i) {
						auto im = i.index();
						im[dim]--;
						if (box.pad(-opts.hbw).contains(i)) {
							S[egas_i][i] -= (data[k] - F[dim][f][i]) * phi[i] / dx * std::pow(cosmos_a(), (NDIM + 2) - NDIM);
						}
						if (box.pad(-opts.hbw).contains(im)) {
							S[egas_i][im] += (data[k] - F[dim][f][i]) * phi[im] / dx * std::pow(cosmos_a(), (NDIM + 2) - NDIM);
						}
					}
					F[dim][f][i] = data[k];
				}
				k++;
			}
		}
	}
	a = data[k];
	k++;
	assert(k == data.size());
	return a;
}

void hydro_grid::unpack_coarse_correction(const std::vector<double> &data, const multi_range &bbox, double dt, double a0, double a1) {
//	printf( "unpacking coarse_correcton\n");
	int k = 0;
	const auto lambda = dt / dx;
	const auto tau0 = U[tau_i];
	for (int dim = 0; dim < NDIM; dim++) {
		auto this_box = bbox;
		this_box.max[dim]++;
		for (int f = 0; f < opts.nhydro; f++) {
			for (multi_iterator i(this_box); !i.end(); i++) {
				assert(k < data.size());
				if (fbox[dim].contains(i)) {
					const auto flux = data[k] / dx - lambda * F0[dim][f][i];
					auto im = i.index();
					im[dim]--;
					if (box.contains(im)) {
						if (f == rho_i && F0[dim][f][i] != 0.0) {
							//						printf( "%e\n", data[k] / dx / lambda / F0[dim][f][i]);
						}
						U[f][im] -= flux;
						if (f == rho_i) {
							U[egas_i][im] += flux * phi[im] * std::pow(cosmos_a(), (NDIM + 2) - NDIM);
						}
					}
					if (box.contains(i)) {
						U[f][i] += flux;
						if (f == rho_i) {
							U[egas_i][i] -= flux * phi[i] * std::pow(cosmos_a(), (NDIM + 2) - NDIM);
						}
					}
				}
				k++;
			}
		}
	}
	const auto a = cosmos_a();
	for (multi_iterator i(box); !i.end(); i++) {
		if (U[tau_i][i] <= 0.0) {

			double ekin = 0.0;
			for (int dim = 0; dim < NDIM; dim++) {
				ekin += 0.5 * std::pow(U[sx_i + dim][i], 2) * std::pow(a, -NDIM - 2) / U[rho_i][i];
			}
			const double eint = U[egas_i][i] * std::pow(a, -(NDIM + 2)) - ekin;
			if (eint < 0.0) {
				U[tau_i][i] = std::max(0.1 * tau0[i], U[tau_i][i]);
			} else {
				U[tau_i][i] = std::pow(a, NDIM) * std::pow(eint, 1.0 / opts.gamma);
			}
		}
	}
	assert(k == data.size());
}

std::vector<double> hydro_grid::pack_coarse_correction() {
	std::vector<double> data;
	for (int dim = 0; dim < NDIM; dim++) {
		for (int f = 0; f < opts.nhydro; f++) {
			for (multi_iterator i(fcbox[dim]); !i.end(); i++) {
				data.push_back(Fc[dim][f][i]);
			}
		}
	}
	return data;
}

void hydro_grid::unpack_refinement(const std::vector<std::uint8_t> &data, multi_range bbox) {
	int k = 0;
	const auto rbox = box.pad(-opts.hbw);
	for (multi_iterator i(bbox); !i.end(); i++) {
		assert(k < data.size());
		R[i] = data[k];
		k++;
	}
}

void hydro_grid::unpack(const std::vector<double> &data, multi_range bbox) {
//	printf("unpacking %s\n", bbox.to_string().c_str());
	int k = 0;
	for (int f = 0; f < opts.nhydro; f++) {
		for (multi_iterator i(bbox); !i.end(); i++) {
			assert(k < data.size());
			U[f][i] = data[k];
			k++;
		}
	}
	for (multi_iterator i(bbox); !i.end(); i++) {
		if (U[rho_i][i] <= 0.0) {
			printf("1 rho is less than zero %e\n", U[rho_i][i]);
			abort();
		}
		if (U[tau_i][i] < 0.0) {
			double ek = 0.0;
			for (int dim = 0; dim < NDIM; dim++) {
				ek += 0.5 * U[sx_i + dim][i] * U[sx_i + dim][i] / U[rho_i][i];
			}
			const auto eint = U[egas_i][i] - ek;
			printf("1 tau is less than zero %e %e\n", U[tau_i][i], eint);
			abort();
		}
	}
	assert(k == data.size());
}

void hydro_grid::unpack_field(int f, const std::vector<double> &data, multi_range bbox) {
//	printf("unpacking %s\n", bbox.to_string().c_str());
	int k = 0;
	for (multi_iterator i(bbox); !i.end(); i++) {
		assert(k < data.size());
		U[f][i] = data[k];
		k++;
	}
}

statistics hydro_grid::get_statistics(const std::vector<multi_range> &child_ranges) const {
	statistics stats;
	const auto a = cosmos_a();
	stats.u.resize(opts.nhydro, 0.0);
	for (multi_iterator i(box.pad(-opts.hbw)); !i.end(); i++) {
		bool refined = false;
		for (const auto &c : child_ranges) {
			if (c.contains(i)) {
				refined = true;
				break;
			}
		}
		if (!refined) {
			for (int f = 0; f < opts.nhydro; f++) {
				stats.u[f] += std::pow(dx, NDIM) * U[f][i];
			}
		}
	}
	return stats;
}

std::vector<std::string> hydro_grid::field_names() {
	std::vector<std::string> names;
	names.push_back("rho");
	names.push_back("E");
	names.push_back("tau");
	names.push_back("pot");
	for (int dim = 0; dim < NDIM; dim++) {
		names.push_back(std::string("s") + char('x' + dim));
	}
	if (opts.species) {
		names.push_back(std::string("rho_H"));
		names.push_back(std::string("rho_Hp"));
		names.push_back(std::string("rho_He"));
		names.push_back(std::string("rho_Hep"));
		names.push_back(std::string("rho_Hepp"));
	}
	names.push_back("phi");
	names.push_back("error");
	return names;
}

