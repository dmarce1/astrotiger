#include <astrotiger/options.hpp>
#include <astrotiger/particles.hpp>
#include <array>

multi_array<double> particles::cloud_in_cell() const {
	multi_array<double> cic;
	const auto dvinv = std::pow(dx, NDIM);
	for (multi_iterator i(box); !i.end(); i++) {
		cic[i] = 0.0;
	}
	for (const auto &p : parts) {
		multi_index i;
		vect<double> q;
		for (int dim = 0; dim < NDIM; dim++) {
			i[dim] = p.x[dim] / dx;
			q[dim] = 1.0 - p.x[dim] / dx + i[dim];
		}
		for (multi_iterator j(multi_range(i).pad(1)); !j.end(); j++) {
			if (box.contains(j)) {
				double wt = 1.0;
				for (int dim = 0; dim < NDIM; dim++) {
					wt *= q[dim];
				}
				cic[j] += p.m * dvinv * wt;
			}
		}
	}
	return cic;
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
	if (opts.problem == "part_test") {
		for (multi_iterator i(box); !i.end(); i++) {
			vect<double> x;
			double r = 0.0;
			for (int dim = 0; dim < NDIM; dim++) {
				x[dim] = i.index()[dim] * dx + 0.5 * dx;
				r += (x[dim] - 0.5) * (x[dim] - 0.5);
			}
			r = std::sqrt(r);
			if (r > 0.25 && r < 0.3) {
				for (int iter = 0; iter < rand() % 5; iter++) {
					particle p;
					double r1 = 0.0;
					for (int dim = 0; dim < NDIM; dim++) {
						p.x[dim] = i.index()[dim] * dx + dx * rand() / RAND_MAX;
						r1 += (p.x[dim] - 0.5) * (p.x[dim] - 0.5);
					}
					p.m = 1.0;
					r1 = std::sqrt(r1);
					p.v[0] = -(p.x[1]-0.5) / std::pow(r1,1.5);
					p.v[1] = (p.x[0]-0.5) / std::pow(r1,1.5);
					p.rung = 0;
					parts.push_back(p);
				}
			}
		}
	}
}

void particles::kick(int kick_level, int this_level, const std::vector<double> &dt0, const std::vector<double> &dt1) {
	for (auto &p : parts) {
		const auto x = (p.x - vect<double>(0.5));
		const auto r = abs(x);
		const auto f = -x / (r * r * r);
		if (p.rung >= kick_level) {
			p.v += f * (p.m * dt0[p.rung] * 0.5);
			if (p.rung > this_level) {
				if (this_level >= kick_level) {
					p.rung = this_level;
				}
			} else if (p.rung < this_level) {
				p.rung = this_level;
			}
			p.v += f * (p.m * dt1[p.rung] * 0.5);
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
