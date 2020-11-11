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

void particles::add_parts(std::vector<particle> &new_parts) {
	for (auto &p : new_parts) {
		for (int dim = 0; dim < NDIM; dim++) {
			while (p.x[dim] >= 1.0) {
				p.x[dim] -= 1.0;
			}
			while (p.x[dim] < 0.0) {
				p.x[dim] += 1.0;
			}
		}
	}
	parts.insert(parts.end(), new_parts.begin(), new_parts.end());
}

void particles::initialize() {
	if (opts.problem == "part_test") {
		for (int i = 0; i < 100; i++) {
			particle p;
			for (int dim = 0; dim < NDIM; dim++) {
				p.x[dim] = (double) rand() / RAND_MAX;
				p.v[dim] = 2.0 * (double) rand() / RAND_MAX - 1.0;
			}
			p.rung = 0;
			parts.push_back(p);
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
			max[dim] = std::max(max[dim], (double) p.v[dim]);
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

std::vector<particle> particles::drift(double dt) {
	std::vector<particle> escaped;
	int i = 0;
	while (i < parts.size()) {
		auto &part = parts[i];
		for (int dim = 0; dim < NDIM; dim++) {
			part.x[dim] += part.v[dim] * dt;
		}
		if (!rbox.contains(part.x)) {
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
