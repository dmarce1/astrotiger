/*
 * grid.cpp
 *
 *  Created on: Oct 14, 2020
 *      Author: dmarce1
 */

#include <astrotiger/hydro_grid.hpp>

hydro_grid::hydro_grid() {
}

void hydro_grid::resize(double dx_, multi_range box_) {
	dx = dx_;
	box = box_;
	const int sz = opts.nhydro * box.volume();
	double *u0 = new double[sz];
	double *u1 = new double[sz];
	int j = 0;
	U0.resize(opts.nhydro);
	U1.resize(opts.nhydro);
	for (int i = 0; i < opts.nhydro; i++) {
		U0[i].resize(box);
		U1[i].resize(box);
	}
	R.resize(box);
}

hydro_grid::~hydro_grid() {
}

void hydro_grid::compute_refinement_criteria() {
	const auto ibox = box.pad(-opts.bw[hydro_i]);
	for (multi_iterator i(box); !i.end(); i++) {
		R[i] = false;
	}
	for (multi_iterator i(ibox); !i.end(); i++) {
		double max_grad_rho = 0.0;
		for (int dim = 0; dim < NDIM; dim++) {
			vect<index_type> ip = i;
			vect<index_type> im = i;
			ip[dim]++;
			im[dim]--;
			max_grad_rho = std::max(max_grad_rho, std::abs(U0[rho_i][ip] - U0[rho_i][im]));
		}
		bool res = max_grad_rho / U0[rho_i][i] > opts.refine_slope;
		if (res) {
			R[i] = true;
			auto window = range<index_type>(i).pad(opts.window);
			for (multi_iterator j(window); !j.end(); j++) {
				R[j] = true;
			}
		}
	}
}

void hydro_grid::initialize() {
	for (multi_iterator i(box); !i.end(); i++) {
		if (opts.problem == "sod") {
			for (int dim = 0; dim < NDIM; dim++) {
				U0[sx_i + dim][i] = 0.0;
			}
			double xsum = 0.0;
			for (int dim = 0; dim < NDIM; dim++) {
				xsum += coord(i[dim]) - 0.5;
			}
			if (xsum < 0.0) {
				U0[rho_i][i] = 1.0;
				U0[egas_i][i] = 1.0;
			} else {
				U0[rho_i][i] = 0.1;
				U0[egas_i][i] = 0.125;
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

std::vector<multi_range> hydro_grid::refined_ranges() const {

	std::vector<multi_range> boxes;
	std::vector<multi_range> tmp;
	std::vector<multi_range> finished;
	multi_range range;
	const auto ibox = box.pad(-opts.bw[hydro_i]);
	range.min = box.max;
	range.max = box.min;
	for (multi_iterator i(ibox); !i.end(); i++) {
		if (R[i]) {
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
		do {
			tmp.resize(0);
			for (const auto &b : boxes) {
				if (b.volume() > max_vol) {
					auto new_boxes = b.split();
					range.min = new_boxes.first.max;
					range.max = new_boxes.first.min;
					for (multi_iterator i(new_boxes.first); !i.end(); i++) {
						if (R[i]) {
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
						if (R[i]) {
							for (int dim = 0; dim < NDIM; dim++) {
								range.min[dim] = std::min(range.min[dim], i[dim]);
								range.max[dim] = std::max(range.max[dim], i[dim] + 1);
							}
						}
					}
					new_boxes.second = range;
					tmp.push_back(new_boxes.first);
					tmp.push_back(new_boxes.second);
				} else {
					finished.push_back(b);
				}
			}
			change = tmp.size();
			boxes = std::move(tmp);
		} while (change);
	}
	return std::move(finished);
}

