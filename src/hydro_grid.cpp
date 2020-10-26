/*
 * grid.cpp
 *
 *  Created on: Oct 14, 2020
 *      Author: dmarce1
 */

#include <astrotiger/hydro_grid.hpp>
#include <astrotiger/hydro_flux.hpp>
#include <astrotiger/multi_array.hpp>
#include <vector>

hydro_grid::hydro_grid() {
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
	}
	return amax;
}

void hydro_grid::substep_update(int rk, double dt) {
	const double onembeta = 1.0 - opts.beta[rk];
	const double lambda = dt / dx;
	if (rk == 0) {
		for (int f = 0; f < opts.nhydro; f++) {
			U0[f] = U[f];
		}
		reset_flux_registers();
	}
	for (multi_iterator i(box.pad(-opts.hbw)); !i.end(); i++) {
		for (int f = 0; f < opts.nhydro; f++) {
			auto &u = U[f][i];
			u = U0[f][i];
			for (int dim = 0; dim < NDIM; dim++) {
				multi_index ip1 = i;
				ip1[dim]++;
				u -= (F[dim][f][ip1] - F[dim][f][i]) * lambda;
			}
		}
	}
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

void hydro_grid::resize(double dx_, multi_range box_) {
	dx = dx_;
	box = box_;
	const int sz = opts.nhydro * box.volume();
	double *u0 = new double[sz];
	double *u1 = new double[sz];
	int j = 0;
	U0.resize(opts.nhydro);
	U.resize(opts.nhydro);
	phi.resize(box_.pad(1));
	for (int i = 0; i < opts.nhydro; i++) {
		U0[i].resize(box);
		U[i].resize(box);
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
		for (int dim = 0; dim < NDIM; dim++) {
			all_zero = all_zero && (j[dim] == 0);
		}
		if (all_zero) {
			continue;
		}
		for (multi_iterator k(box.pad(-opts.hbw)); !k.end(); k++) {
			if (!R[k]) {
				for (int f = 0; f < opts.nhydro; f++) {
					const auto up = U[f][k.index() + j];
					const auto u0 = U[f][k];
					const auto um = U[f][k.index() - j];
					const auto num = std::abs(up + um - 2.0 * u0);
					const auto den = std::abs(up - u0) + std::abs(u0 - um) + 0.01 * (std::abs(up) + std::abs(um) + 2.0 * std::abs(u0));
					if (den > 0.0) {
						if (num / den > opts.refine_slope) {
							R[k] = true;
						}
					}
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
			for (int dim = 0; dim < NDIM; dim++) {
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
			const double rnd = 2.0 * (double) (rand() + 0.5) / RAND_MAX - 1.0;
			if (std::abs(coord(i.index()[dim]) - 0.5) > 0.25) {
				U[rho_i][i] = 1.0;
				U[sx_i][i] = +0.5 + rnd * 0.001;
			} else {
				U[sx_i][i] = -0.5;
				U[rho_i][i] = 2.0 + 2.0 * rnd * 0.001;
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
			if (r < 0.1) {
				U[rho_i][i] = 1.0;
			} else {
				U[rho_i][i] = 1.0e-20;
			}
			U[egas_i][i] = 1e-20;
		} else {
			printf("unknown problem %s\n", opts.problem.c_str());
			abort();
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
		for (const auto &amr : amr_boxes) {
			tmp.resize(0);
			for (const auto &b : boxes) {
				const auto amrbox = amr.pad(1);
				if (!amrbox.intersection(b).empty()) {
					auto tmp2 = b.subtract(amrbox);
					tmp.insert(tmp.end(), tmp2.begin(), tmp2.end());
				} else {
					tmp.push_back(b);
				}
			}
			boxes = std::move(tmp);
		}
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

std::vector<std::vector<double>> hydro_grid::pack_output() const {
	std::vector<std::vector<double>> data;
	auto bbox = box.pad(-opts.hbw);
	data.resize(opts.nhydro + 1);
	if ( NDIM > 1) {
		std::swap(bbox.min[0], bbox.min[NDIM - 1]);
		std::swap(bbox.max[0], bbox.max[NDIM - 1]);
	}
	for (int f = 0; f <= opts.nhydro; f++) {
		const auto &u = f == opts.nhydro ? phi : U[f];
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

void hydro_grid::unpack_coarse_flux(const std::vector<double> &data, const multi_range &bbox, double dt) {
	int k = 0;
	const auto lambda = dt / dx;
	for (int dim = 0; dim < NDIM; dim++) {
		auto this_box = bbox;
		this_box.max[dim]++;
		for (int f = 0; f < opts.nhydro; f++) {
			for (multi_iterator i(this_box); !i.end(); i++) {
				assert(k < data.size());
				const auto flux = data[k] / dx - lambda * F[dim][f][i];
				//			printf("%e %e\n", data[k] / dx, lambda * F[dim][f][i]);
				auto im = i.index();
				im[dim]--;
				U[f][i] += flux;
				U[f][im] -= flux;
				k++;
			}
		}
	}
}

std::vector<double> hydro_grid::pack_coarse_flux() {
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

statistics hydro_grid::get_statistics(const std::vector<multi_range>& child_ranges) const {
	statistics stats;
	stats.u.resize(opts.nhydro,0.0);
	for( multi_iterator i(box.pad(-opts.hbw)); !i.end(); i++) {
		bool refined = false;
		for( const auto& c : child_ranges) {
			if( c.contains(i)) {
				refined = true;
				break;
			}
		}
		if(!refined) {
			for( int f = 0; f < opts.nhydro; f++) {
				stats.u[f] += std::pow(dx,2) * U[f][i];
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
	return names;
}

