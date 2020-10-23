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
	const auto beta = rk == 0 ? 1.0 : 0.5;
	std::vector<multi_array<double>> UR(opts.nhydro);
	std::vector<multi_array<double>> UL(opts.nhydro);
	const auto urlbox = box.pad(2 - opts.hbw);
	const auto grad_box = box.pad(-opts.hbw + 1);
	for (int i = 0; i < opts.nhydro; i++) {
		UR[i].resize(urlbox);
		UL[i].resize(urlbox);
	}
	for (int dim = 0; dim < NDIM; dim++) {
		for (int i = 0; i < opts.nhydro; i++) {
			auto this_box = fbox[dim];
			this_box.min[dim]--;
			for (multi_iterator j(this_box); !j.end(); j++) {
				const auto du = 0.5 * U[i].minmod_gradient(dim, j);
				UR[i][j] = U[i][j] - du;
				auto jp1 = j.index();
				jp1[dim]++;
				UL[i][jp1] = U[i][j] + du;
			}
		}
		std::vector<double> ur(opts.nhydro), ul(opts.nhydro), flux(opts.nhydro);
		for (multi_iterator j(fbox[dim]); !j.end(); j++) {
			for (int f = 0; f < opts.nhydro; f++) {
				ur[f] = UR[f][j];
				ul[f] = UL[f][j];
			}
			amax = std::max(amax, hydro_flux(flux, ul, ur, dim));
			for (int f = 0; f < opts.nhydro; f++) {
				F[dim][f][j] = (1.0 - beta) * F[dim][f][j] + beta * flux[f];
			}
		}
	}
	return amax;
}

void hydro_grid::substep_update(int rk, double dt) {
	const double beta = rk == 0 ? 1.0 : 0.5;
	const double onembeta = 1.0 - beta;
	const double lambda = dt / dx;
	if (rk == 0) {
		for (int f = 0; f < opts.nhydro; f++) {
			U0[f] = U[f];
		}
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
	if (rk == opts.nrk - 1) {
		reset_coarse_flux_registers();
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
	for (multi_iterator i(box.pad(-opts.hbw)); !i.end(); i++) {
		double max_grad_rho = 0.0;
		double max_grad_egas = 0.0;
		for (int dim = 0; dim < NDIM; dim++) {
			max_grad_rho = std::max(max_grad_rho, std::abs(U[rho_i].smooth_gradient(dim, i)));
			max_grad_egas = std::max(max_grad_egas, std::abs(U[egas_i].smooth_gradient(dim, i)));
		}
		bool res = max_grad_rho / U[rho_i][i] > opts.refine_slope;
		res = res || max_grad_egas / U[egas_i][i] > opts.refine_slope;
		R[i] = res;
	}
}

void hydro_grid::initialize() {
	for (multi_iterator i(box); !i.end(); i++) {
		if (opts.problem == "sod") {
			for (int dim = 0; dim < NDIM; dim++) {
				U[sx_i + dim][i] = 0.0;
			}
			double xsum = 0.0;
			for (int dim = 0; dim < 1; dim++) {
				xsum += coord(i[dim]) - 0.5;
			}
			if (xsum > 0.0) {
				U[rho_i][i] = 1.0;
				U[egas_i][i] = 1.0;
			} else {
				U[rho_i][i] = 0.1;
				U[egas_i][i] = 0.125;
			}
		} else if (opts.problem == "blast") {
			double r = 0.0;
			U[rho_i][i] = 1.0;
			for (int dim = 0; dim < NDIM; dim++) {
				U[sx_i + dim][i] = 0.0;
				r += std::pow(coord(i[dim]) - 0.5, 2);
			}
			r = std::sqrt(r);
			if (r < 0.1) {
				U[egas_i][i] = 1.0;
			} else {
				U[egas_i][i] = 1.0e-3;
			}
		} else {
			printf("unknown problem %s\n", opts.problem.c_str());
			abort();
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
				if (this_vol > max_vol || (efficiency < opts.efficiency && this_vol >= min_vol)) {
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
				const auto amrbox = amr.pad(opts.window);
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
//	printf("Packing %s\n", bbox.to_string().c_str());
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
	if ( NDIM > 1) {
		std::swap(bbox.min[0], bbox.min[NDIM - 1]);
		std::swap(bbox.max[0], bbox.max[NDIM - 1]);
	}
	for (multi_iterator i(bbox); !i.end(); i++) {
		multi_index j = i.index();
		if ( NDIM > 1) {
			std::swap(j[0], j[NDIM - 1]);
		}
		data.push_back(U[f][j]);
	}
	return data;
}

std::vector<double> hydro_grid::pack_prolong(multi_range bbox, double w) const {
	std::vector<double> data;
	for (int f = 0; f < opts.nhydro; f++) {
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

std::vector<std::string> hydro_grid::field_names() {
	std::vector<std::string> names;
	names.push_back("D");
	names.push_back("E");
	for (int dim = 0; dim < NDIM; dim++) {
		names.push_back(std::string("S") + char('x' + dim));
	}
	return names;
}

