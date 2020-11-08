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
#include <astrotiger/tree.hpp>
#include <vector>

double rand1() {
	return (rand() + 0.5) / RAND_MAX;
}

hydro_grid::hydro_grid() {
}

void hydro_grid::to_array(multi_array<double> &a, const multi_range &bbox, int f, double w) const {
	for (multi_iterator i(bbox); !i.end(); i++) {
		a[i] = w * U[f][i] + (1.0 - w) * U0[f][i];
	}
}

void hydro_grid::store() {
	U0 = U;
}

double hydro_grid::compute_flux(int rk) {
	double amax = 0.0;
	std::vector<multi_array<double>> VR(opts.nhydro);
	std::vector<multi_array<double>> VL(opts.nhydro);
	const auto urlbox = box.pad(2 - opts.hbw);
	const auto grad_box = box.pad(-opts.hbw + 1);
	std::vector<multi_array<double>> V(opts.nhydro);
	for (int i = 0; i < opts.nhydro; i++) {
		VR[i].resize(urlbox);
		VL[i].resize(urlbox);
		V[i].resize(box);
	}
	for (multi_iterator i(box); !i.end(); i++) {
		if( U[rho_i][i] < 0.0 )
		printf( "%e %i %i %s\n", U[rho_i][i], i.index()[0], i.index()[1], box.to_string().c_str());
		assert(U[rho_i][i] > 0.0);
		assert(U[tau_i][i] >= 0.0);
		const auto rhoinv = 1.0 / U[rho_i][i];
		double ek = 0.0;
		V[rho_i][i] = U[rho_i][i];
		for (int dim = 0; dim < NDIM; dim++) {
			ek += 0.5 * std::pow(U[sx_i + dim][i], 2) * rhoinv;
			V[sx_i + dim][i] = U[sx_i + dim][i] * rhoinv;
		}
		V[egas_i][i] = U[egas_i][i] - ek;
		V[tau_i][i] = std::pow(U[tau_i][i], opts.gamma);

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
			amax = std::max(amax, hydro_flux(flux, ul, ur, dim));
			for (int f = 0; f < opts.nhydro; f++) {
				F[dim][f][j] = (1.0 - opts.beta[rk]) * F[dim][f][j] + opts.beta[rk] * flux[f];
			}
		}
		for (multi_iterator j(box.pad(-opts.hbw)); !j.end(); j++) {
			if (opts.gravity) {
				auto jp = j.index();
				auto jm = j.index();
				jp[dim]++;
				jm[dim]--;
				S[sx_i + dim][j] = (1.0 - opts.beta[rk]) * S[sx_i + dim][j] + opts.beta[rk] * 0.5 * U[rho_i][j] * (-phi[jp] + phi[jm]) / dx;
			}
		}
	}
	return amax;
}

void hydro_grid::substep_update(int rk, double dt) {
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

void hydro_grid::compute_drho_dt() {
	const double lambda = 1.0 / dx;
	for (multi_iterator i(box.pad(-opts.hbw)); !i.end(); i++) {
		auto& drho = drho_dt[i];
		drho = 0.0;
		for (int dim = 0; dim < NDIM; dim++) {
			multi_index ip1 = i;
			ip1[dim]++;
			drho -= (F[dim][rho_i][ip1] - F[dim][rho_i][i]) * lambda;
		}
	}
}

void hydro_grid::store_flux() {
	F0 = F;
}

void hydro_grid::resize(double dx_, multi_range box_) {
	dx = dx_;
	box = box_;
	const int sz = opts.nhydro * box.volume();
	double *u0 = new double[sz];
	double *u1 = new double[sz];
	int j = 0;
	U0.resize(opts.nhydro);
	drho_dt.resize(box);
	U.resize(opts.nhydro);
	S.resize(opts.nhydro);
	phi.resize(box_.pad(opts.hbw));
	error.resize(box_.pad(opts.hbw));
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

void hydro_grid::compute_refinement_criteria() {
	const auto ibox = box.pad(-opts.hbw + opts.window);
	for (multi_iterator i(ibox); !i.end(); i++) {
		R[i] = false;
	}
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
			if (!R[k]) {
				constexpr double eps = 0.1;
				for (int f = 0; f < opts.nhydro; f++) {
					const auto up = U[f][k.index() + j];
					const auto u0 = U[f][k];
					const auto um = U[f][k.index() - j];
					double c0;
					if (f < sx_i || f >= sx_i + NDIM) {
						c0 = eps * (std::abs(up) + std::abs(um) + 2.0 * std::abs(u0));
					} else {
						double vp = 0.0;
						double v0 = 0.0;
						double vm = 0.0;
						for (int dim = 0; dim < NDIM; dim++) {
							vm += std::max(0.0, U[egas_i][k.index() - j] / U[rho_i][k.index() - 1]);
							v0 += std::max(0.0, U[egas_i][k] / U[rho_i][k]);
							vp += std::max(0.0, U[egas_i][k.index() + j] / U[rho_i][k.index() + 1]);
						}
						c0 = eps * (std::sqrt(vp) + std::sqrt(vm) + 2.0 * std::sqrt(v0));
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

void hydro_grid::initialize() {
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
				if (dim == 0) {
					if (x > 0.5) {
						x -= 1.0;
					}
				} else {
					x -= 0.5;
				}
				r += std::pow(x, 2);
			}
			r = std::sqrt(r);
			const auto n = 1.5;
			const auto alpha = 0.005;
			const auto theta = lane_emden(r / alpha, dx / alpha / 2.0, n);
			auto &rho = U[rho_i][i];
			rho = std::max(1.0e-10, std::pow(theta, n));
			const auto vx = 1.0;
			const auto vy = 0.00;
			U[sx_i][i] = vx * rho;
			U[sy_i][i] = vy * rho;
			const auto K = 4.0 * M_PI * opts.G * alpha * alpha / (n + 1);
			U[egas_i][i] = std::max(1.0e-20, K * std::pow(theta, 1.0 + n) / (opts.gamma - 1.0)); // * std::pow(unit * unit, opts.gamma);
			U[tau_i][i] = std::pow(U[egas_i][i], 1.0 / opts.gamma);
			U[egas_i][i] += 0.5 * (vx * vx + vy * vy) * rho;
		} else {
			printf("unknown problem %s\n", opts.problem.c_str());
			abort();
		}

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
	for (multi_iterator i(box.pad(-opts.hbw)); !i.end(); i++) {
		double max_egas = 0.0;
		double ekin = 0.0;
		for (multi_iterator d(multi_range(0).pad(1)); !d.end(); d++) {
			max_egas = std::max(max_egas, U[egas_i][i]);
		}
		for (int dim = 0; dim < NDIM; dim++) {
			ekin += 0.5 * std::pow(U[sx_i + dim][i], 2) / U[rho_i][i];
		}
		const double eint = U[egas_i][i] - ekin;
		if (eint > U[egas_i][i] * 0.1) {
			U[tau_i][i] = std::pow(eint, 1.0 / opts.gamma);
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
	if (!range.empty()) {
		boxes.push_back(range);
		for (const auto &amr : amr_boxes) {
			tmp.resize(0);
			for (const auto &b : boxes) {
				const auto amrbox = amr.pad(opts.window);
				if (amrbox.intersection(b).volume()) {
					auto tmp2 = b.subtract(amrbox);
					tmp.insert(tmp.end(), tmp2.begin(), tmp2.end());
				} else {
					tmp.push_back(b);
				}
			}
			boxes = std::move(tmp);
		}
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

std::vector<std::vector<double>> hydro_grid::pack_output() const {
	std::vector<std::vector<double>> data;
	auto bbox = box.pad(-opts.hbw);
//	auto bbox = box;
	data.resize(opts.nhydro + 2);
	if ( NDIM > 1) {
		std::swap(bbox.min[0], bbox.min[NDIM - 1]);
		std::swap(bbox.max[0], bbox.max[NDIM - 1]);
	}
	for (int f = 0; f <= opts.nhydro + 1; f++) {
		const auto &u = f == opts.nhydro ? phi : (f == opts.nhydro + 1 ? error : U[f]);
		for (multi_iterator i(bbox); !i.end(); i++) {
			multi_index j = i.index();
			if ( NDIM > 1) {
				std::swap(j[0], j[NDIM - 1]);
			}
			data[f].push_back(u[j]);
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
//
//std::vector<double> hydro_grid::pack_coarse_flux() {
//	const double a = compute_flux(0);
//	const auto inv = 1.0 / (1 << (NDIM - 1));
//	std::vector<double> data;
//	int size = 0;
//	for (int dim = 0; dim < NDIM; dim++) {
//		auto cbox = box.pad(-opts.gbw).half();
//		cbox.max[dim]++;
//		size += cbox.volume();
//	}
//	printf( "In\n");
//	data.reserve(size);
//	for (int dim = 0; dim < NDIM; dim++) {
//		auto cbox = box.pad(-opts.gbw).half();
//		cbox.max[dim]++;
//		for (int f = 0; f < opts.nhydro; f++) {
//			for (multi_iterator i(cbox); !i.end(); i++) {
//				double this_flux = 0.0;
//				auto this_box = multi_range(i.index()).double_();
//				this_box.max[dim]--;
//				for (multi_iterator j(this_box); !j.end(); j++) {
//					this_flux += F[dim][f][j] * inv;
//				}
//				data.push_back(this_flux);
//			}
//		}
//	}
//	printf( "out\n");
//	data.push_back(a);
//	return data;
//}
//
//double hydro_grid::unpack_fine_flux(const std::vector<double> &data, const multi_range &bbox) {
//	int k = 0;
//	double a;
//	for (int dim = 0; dim < NDIM; dim++) {
//		auto this_box = bbox;
//		this_box.max[dim]++;
//		for (int f = 0; f < opts.nhydro; f++) {
//			for (multi_iterator i(this_box); !i.end(); i++) {
//				assert(k < data.size());
////				if (f == rho_i) {
////					printf("%e %e\n", F[dim][f][i], data[k]);
////				}
//				F[dim][f][i] = data[k];
//				k++;
//			}
//		}
//	}
//	a = data[k];
//	k++;
//	assert(k == data.size());
//	return a;
//}
//
void hydro_grid::unpack_coarse_correction(const std::vector<double> &data, const multi_range &bbox, double dt) {
	int k = 0;
	const auto lambda = dt / dx;
	for (int dim = 0; dim < NDIM; dim++) {
		auto this_box = bbox;
		this_box.max[dim]++;
		for (int f = 0; f < opts.nhydro; f++) {
			for (multi_iterator i(this_box); !i.end(); i++) {
				assert(k < data.size());
				const auto flux = data[k] / dx - lambda * F0[dim][f][i];
				auto im = i.index();
				im[dim]--;
				U[f][i] += flux;
				U[f][im] -= flux;
				k++;
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
	reset_coarse_flux_registers();
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
	names.push_back("D");
	names.push_back("E");
	names.push_back("tau");
	for (int dim = 0; dim < NDIM; dim++) {
		names.push_back(std::string("S") + char('x' + dim));
	}
	names.push_back("phi");
	names.push_back("error");
	return names;
}

