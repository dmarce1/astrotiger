#include <astrotiger/options.hpp>
#include <astrotiger/particles.hpp>
#include <astrotiger/rand.hpp>
#include <array>

multi_array<double> particles::cloud_in_cell(double dt) const {
	multi_array<double> src(box.pad(1));
	multi_array<double> rho(box.pad(1));
	multi_array<double> drho_dt(box.pad(1));
	const auto dvinv = std::pow(dx, -NDIM);
	for (multi_iterator i(box.pad(1)); !i.end(); i++) {
		rho[i] = 0.0;
		drho_dt[i] = 0.0;
	}
	for (const auto &p : parts) {
		multi_index i;
		vect<double> q;
		for (int dim = 0; dim < NDIM; dim++) {
			const double tmp = (p.x[dim] - 0.5 * dx + opts.max_bw * dx) / dx - opts.max_bw;
			i[dim] = int((p.x[dim] - 0.5 * dx + opts.max_bw * dx) / dx) - opts.max_bw;
			;
			q[dim] = 1.0 - tmp + i[dim];
		}
		auto this_box = multi_range(i);
		for (int dim = 0; dim < NDIM; dim++) {
			this_box.max[dim]++;
		}
		for (multi_iterator j(this_box); !j.end(); j++) {
			if (box.contains(j)) {
				double wt = 1.0;
				double sgn;
				for (int dim = 0; dim < NDIM; dim++) {
					if (j[dim] == i[dim]) {
						wt *= q[dim];
						sgn = -1.0;
					} else {
						wt *= 1.0 - q[dim];
						sgn = +1.0;
					}
					drho_dt[j] += p.m * dvinv * sgn * p.v[dim] / dx;
				}
				assert(wt >= 0.0);
				assert(wt <= 1.0);
				rho[j] += p.m * dvinv * wt;
			}
		}
	}
	for (multi_iterator i(box.pad(1)); !i.end(); i++) {
		rho[i] += dt * drho_dt[i];
	}
	return rho;
}

multi_array<int> particles::particle_count() const {
	multi_array<int> count(box);
	for (multi_iterator i(box); !i.end(); i++) {
		count[i] = 0;
	}
	for (const auto &p : parts) {
		multi_index index;
		for (int dim = 0; dim < NDIM; dim++) {
			index[dim] = p.x[dim] / dx;
		}
		count[index]++;
	}
	return count;
}

std::vector<particle> particles::get_child_parts() {
	std::vector<particle> cparts;
	for (const auto &b : child_boxes) {
		int i = 0;
		while (i < parts.size()) {
			if (b.contains(parts[i].x)) {
				cparts.push_back(parts[i]);
				const int sz = parts.size() - 1;
				parts[i] = parts[sz];
				parts.resize(sz);
			} else {
				i++;
			}
		}
	}
	return cparts;
}

void particles::set_child_boxes(const std::vector<multi_range> &boxes) {
	child_boxes.resize(0);
	for (const auto &b : boxes) {
		range<double> rbox;
		for (int dim = 0; dim < NDIM; dim++) {
			rbox.min[dim] = b.min[dim] * dx;
			rbox.max[dim] = b.max[dim] * dx;
		}
		child_boxes.push_back(rbox);
	}
}

void particles::initialize() {
	constexpr int N = 10000;
	constexpr double a = 0.01;
	constexpr double c0 = 3.0 * 1.0 / (4.0 * M_PI * a * a * a);
	if (opts.problem == "part_test") {
		for (int i = 0; i < N; i++) {
			particle p;
			double r;
			double prob;
			double rnd;
			vect<double> x;
			do {
				r = 0.0;
				for (int dim = 0; dim < NDIM; dim++) {
					x[dim] = rand1();
					r += std::pow(x[dim] - 0.5, 2);
				}
				r = std::sqrt(r);
				prob = c0 * std::pow(1.0 + r * r / a / a, -2.5);
			} while (prob < rand1());
			p.x = x;
			const auto sigma = std::sqrt(1.0 / 6.0 / std::sqrt(1.0 + r * r / a / a));
			for (int dim = 0; dim < NDIM; dim++) {
				p.v[dim] = (2.0 * rand1() - 1.0) * sigma;
			}
			p.m = 1.0 / N;
			p.rung = 0;
			parts.push_back(p);
		}
	}
}

void particles::kick(int kick_level, int this_level, const std::vector<double> &dt0, const std::vector<double> &dt1,
		const std::array<multi_array<double>, NDIM> &g) {
	for (auto &p : parts) {
		if (p.rung >= kick_level) {
			multi_index i;
			vect<double> this_g;
			for (int dim = 0; dim < NDIM; dim++) {
				i[dim] = int(((p.x[dim] - 0.5 * dx) + opts.max_bw * dx) / dx) - opts.max_bw;
			}
			for (int dim = 0; dim < NDIM; dim++) {
				this_g[dim] = 0.0;
				multi_range this_box(i);
				for (int dim2 = 0; dim2 < NDIM; dim2++) {
					this_box.max[dim2]++;
				}
				for (multi_iterator j(this_box); !j.end(); j++) {
					double w = 1.0;
					for (int dim2 = 0; dim2 < NDIM; dim2++) {
						const auto tmp = ((p.x[dim2] - 0.5 * dx) + opts.max_bw * dx) / dx - opts.max_bw;
						if (j.index()[dim2] != i[dim2]) {
							w *= tmp - i[dim2];
						} else {
							w *= 1.0 - tmp + i[dim2];
						}
					}
					//		printf("%e %e\n", g[dim][j], w);
					assert(w >= 0.0);
					assert(w <= 1.0);
					this_g[dim] += w * g[dim][j];
				}

			}
//			printf("%e %e\n", this_g[0], this_g[1]);
			p.v += this_g * (dt0[p.rung] * 0.5);
			if (p.rung > this_level) {
				if (this_level >= kick_level) {
					p.rung = this_level;
				}
			} else if (p.rung < this_level) {
				p.rung = this_level;
			}
			p.v += this_g * (dt1[p.rung] * 0.5);
		}
	}
}

double particles::max_velocity() const {
	std::array<double, NDIM> max;
	for (int dim = 0; dim < NDIM; dim++) {
		max[dim] = 0.0;
	}
	for (const auto &p : parts) {
		for (int dim = 0; dim < NDIM; dim++) {
			max[dim] = std::max(max[dim], std::abs((double) p.v[dim]));
		}
	}
	for (int dim = 1; dim < NDIM; dim++) {
		max[0] = std::max(max[0], max[dim]);
	}
	return max[0];
}

std::vector<std::vector<double>> particles::pack_output() const {
	std::vector<std::vector<double>> data(NDIM + 1);
	for (const auto &part : parts) {
		for (int dim = 0; dim < NDIM; dim++) {
			data[dim].push_back(part.v[dim]);
		}
		data[NDIM].push_back(part.rung);
	}
	return data;
}

std::vector<std::vector<double>> particles::pack_coords() const {
	std::vector<std::vector<double>> coords(NDIM);
	for (const auto &part : parts) {
		for (int dim = 0; dim < NDIM; dim++) {
			coords[dim].push_back(part.x[dim]);
		}
	}
	return coords;
}

std::vector<std::string> particles::field_names() {
	std::vector<std::string> names;
	for (int dim = 0; dim < NDIM; dim++) {
		names.push_back(std::string("v_") + char('x' + char(dim)));
	}
	names.push_back("rung");
	return names;
}

std::vector<particle> particles::get_particles() {
	return parts;
}

void particles::add_parts(const std::vector<particle> &new_parts) {
#ifndef NDEBUG
	for (auto &p : new_parts) {
		assert(rbox.contains(p.x));
	}
#endif
	parts.insert(parts.end(), new_parts.begin(), new_parts.end());
}

std::vector<particle> particles::get_particles(const multi_range &bbox) {
	std::vector<particle> rparts;
	range<double> this_rbox;
	for (int dim = 0; dim < NDIM; dim++) {
		this_rbox.min[dim] = bbox.min[dim] * dx;
		this_rbox.max[dim] = bbox.max[dim] * dx;
	}
	int i = 0;
	while (i < parts.size()) {
		assert(rbox.contains(parts[i].x));
		if (this_rbox.contains(parts[i].x)) {
			rparts.push_back(parts[i]);
			const int sz = parts.size() - 1;
			parts[i] = parts[sz];
			parts.resize(sz);
		} else {
			i++;
		}
	}
	return rparts;
}

std::vector<particle> particles::drift(double dt) {
	std::vector<particle> escaped;
	int i = 0;
	while (i < parts.size()) {
		auto &part = parts[i];
		for (int dim = 0; dim < NDIM; dim++) {
			part.x[dim] += part.v[dim] * dt;
		}
		bool leave_box = !rbox.contains(part.x);
		if (!leave_box) {
			for (const auto &cbox : child_boxes) {
				leave_box = cbox.contains(part.x);
				if (leave_box) {
					break;
				}
			}
		}
		if (leave_box) {
			escaped.push_back(part);
			const auto sz = parts.size() - 1;
			parts[i] = parts[sz];
			parts.resize(sz);
		} else {
			i++;
		}
	}
	return escaped;
}

void particles::resize(double dx_, const multi_range &box_) {
	dx = dx_;
	box = box_;
	for (int dim = 0; dim < NDIM; dim++) {
		rbox.min[dim] = box.min[dim] * dx;
		rbox.max[dim] = box.max[dim] * dx;
	}
}
