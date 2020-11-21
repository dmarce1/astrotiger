#include <astrotiger/options.hpp>
#include <astrotiger/particles.hpp>
#include <astrotiger/rand.hpp>
#include <astrotiger/cosmos.hpp>
#include <astrotiger/fileio.hpp>
#include <astrotiger/tree.hpp>
#include <array>

std::vector<double> particles::pack_cic_restrict() const {
	std::vector<double> data;
	const auto restrict_box = box.pad(2).half();
	multi_array<double> rho_restrict = rho.restrict_(restrict_box);
	for (multi_iterator i(restrict_box); !i.end(); i++) {
		data.push_back(rho_restrict[i]);
	}
	return data;
}

void particles::zero_cic(const multi_range &bbox) {
	for (multi_iterator i(bbox); !i.end(); i++) {
		rho[i] = 0.0;
	}
}

multi_array<double> particles::get_cic() const {
	return rho;
}

void particles::unpack_cic(const std::vector<double> &data, const multi_range &bbox) {
	int k = 0;
	for (multi_iterator i(bbox); !i.end(); i++) {
		assert(k < data.size());
//		assert( data[k] >= 0.0);
		rho[i] += data[k];
		k++;
	}
	assert(k == data.size());
}

std::vector<double> particles::pack_cic(const multi_range &bbox) const {
	std::vector<double> data;
	for (multi_iterator i(bbox); !i.end(); i++) {
		data.push_back(rho[i]);
	}
	return data;
}

void particles::unpack_cic_prolong(const std::vector<double> &data, const multi_range &bbox) {
	int k;
	for (multi_iterator i(bbox); !i.end(); i++) {
		assert(k < data.size());
		for (multi_iterator j(multi_range(i).double_()); !j.end(); j++) {
			rho[j] += data[k];
		}
		k++;
	}
	assert(k == data.size());
}

void particles::compute_cloud_in_cell(double dt) {
	const auto a = cosmos_a();
	const auto adot = cosmos_adot();
//	printf( "%i %s\n", parts.size(), box.to_string().c_str());
	const auto rhobox = box.pad(2);
	rho.resize(rhobox);
	multi_array<double> drho_dt(box.pad(2));
	const auto dvinv = std::pow(dx, -NDIM);
	for (multi_iterator i(rhobox); !i.end(); i++) {
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
		double wt_sum = 0.0;
		for (multi_iterator j(this_box); !j.end(); j++) {
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
			}
			assert(rbox.contains(p.x));
			if (rhobox.contains(j)) {
				double sgn;
				for (int dim = 0; dim < NDIM; dim++) {
					if (j[dim] == i[dim]) {
						sgn = -1.0;
					} else {
						sgn = +1.0;
					}
					drho_dt[j] += p.m * dvinv * sgn * p.v[dim] / (a * a * dx);
				}
				assert(wt >= 0.0);
				assert(wt <= 1.0);
				rho[j] += p.m * dvinv * wt;
			}
			wt_sum += wt;
		}
		if (std::abs(1.0 - wt_sum) > 1.0e-14) {
			printf("%e\n", wt_sum - 1.0);
			abort();
		}
	}

	for (multi_iterator i(box.pad(1)); !i.end(); i++) {
		rho[i] += dt * drho_dt[i];
	}
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
	if (opts.problem == "cosmos") {
		parts = fileio_get_particles();
		for (auto &p : parts) {
			p.m *= (opts.omega_m - opts.omega_b) / opts.omega_m;
		}
	}
}

void particles::kick(int kick_level, int this_level, const std::vector<double> &dt0, const std::vector<double> &dt1,
		const std::array<multi_array<double>, NDIM> &g) {
	const auto a = cosmos_a();
	const auto adot = cosmos_adot();
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
			auto dt = dt0[p.rung];
			p.v += this_g * 0.5 * dt * a;
			if (p.rung > this_level) {
				if (this_level >= kick_level) {
					p.rung = this_level;
				}
			} else if (p.rung < this_level) {
				p.rung = this_level;
			}
			dt = dt1[p.rung];
			p.v += this_g * 0.5 * dt * a;
		}
	}
}

energy_statistics particles::get_energy_statistics(const multi_array<double> &phi) const {
	energy_statistics e;
	e.epot = 0.0;
	e.ekin = 0.0;
	const auto a = cosmos_a();
	for (auto &p : parts) {
		multi_index i;
		vect<double> this_g;
		for (int dim = 0; dim < NDIM; dim++) {
			i[dim] = int(((p.x[dim] - 0.5 * dx) + opts.max_bw * dx) / dx) - opts.max_bw;
		}
		double this_phi = 0.0;
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
			assert(w >= 0.0);
			assert(w <= 1.0);
			this_phi += w * phi[j];
		}
		e.epot += 0.5 * p.m * this_phi;
		e.ekin += 0.5 * p.m * p.v.dot(p.v) / (a * a);
	}
	return e;
}

double particles::max_velocity() const {
	const auto a = cosmos_a();
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
	return max[0] / a;
}

std::vector<std::vector<double>> particles::pack_output() const {
	std::vector<std::vector<double>> data(NDIM + 2);
	for (const auto &part : parts) {
		for (int dim = 0; dim < NDIM; dim++) {
			data[dim].push_back(part.v[dim]);
		}
		data[NDIM].push_back(part.rung);
		data[NDIM + 1].push_back(part.m);
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
	names.push_back("m");
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
	const auto a = cosmos_a();
	std::vector<particle> escaped;
	int i = 0;
	while (i < parts.size()) {
		auto &part = parts[i];
		for (int dim = 0; dim < NDIM; dim++) {
			part.x[dim] += part.v[dim] * dt / (a * a);
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
