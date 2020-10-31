#include <astrotiger/gravity.hpp>
#include <astrotiger/options.hpp>

void gravity::resize(double dx_, const multi_range &box_) {
	dx = dx_;
	const auto cbox = box_.half().pad(opts.gbw);
	box = box_.pad(opts.gbw);
	phi0.resize(box);
	phi1.resize(box);
	phi.resize(box);
	X.resize(box);
	R.resize(box);
	active.resize(box);
	resid.resize(box);
	amr.resize(box);
	phi_c.resize(cbox);
	refined.resize(box);
	for (int dim = 0; dim < NDIM; dim++) {
		fbox[dim] = box_;
		fbox[dim].max[dim]++;
		flux[dim].resize(fbox[dim]);
	}

}

void gravity::set_refined(const std::vector<multi_range> &boxes) {
	for (multi_iterator i(box.pad(-1)); !i.end(); i++) {
		refined[i] = false;
	}
	for (multi_iterator i(box); !i.end(); i++) {
		for (const auto &b : boxes) {
			refined[i] = refined[i] || b.contains(i);
		}
	}
}

void gravity::set_amr_zones(const std::vector<multi_range> &boxes, const std::vector<multi_range> &boxes2, const std::vector<double> &data) {
	const auto cbox = box.pad(-opts.gbw).half().pad(opts.gbw);
	int k = 0;
	for (multi_iterator i(cbox); !i.end(); i++) {
		assert(k < data.size());
		phi_c[i] = data[k];
		k++;
	}

	for (multi_iterator i(box); !i.end(); i++) {
		if (!box.pad(-opts.gbw).contains(i)) {
			active[i] = false;
		}
	}
	for (const auto &b : boxes) {
		const auto inter = b.intersection(box.pad(-1));
		if (inter.volume()) {
			for (multi_iterator i(inter); !i.end(); i++) {
				active[i] = true;
			}
		}
	}

	assert(k == data.size());
	for (multi_iterator i(box); !i.end(); i++) {
		amr[i] = false;
		for (const auto &b : boxes2) {
			amr[i] = b.contains(i);
			if (amr[i]) {
				break;
			}
		}
	}
	compute_amr_bounds(true);
}
void gravity::compute_amr_bounds(bool plus_interior) {
	const auto to_dir = [](multi_index i) {
		vect<index_type> j(0);
		for (int dim = 0; dim < NDIM; dim++) {
			if (i[dim] % 2 == 0) {
				j[dim]--;
			} else {
				j[dim]++;
			}
		}
		return j;
	};
	if (plus_interior) {
		for (multi_iterator I(box); !I.end(); I++) {
			const auto i = I.index();
			const auto d = to_dir(i);
			const auto ic = i / 2;
			const auto icp = ic + d;
			const auto icm = ic - d;
			X[i] = (-(3.0 / 32.0) * phi_c[icm] + (15.0 / 16.0) * phi_c[ic] + (5.0 / 32.0) * phi_c[icp]);
			//		X[i]= phi_c[ic];
		}
	}
	const auto ibox = box.pad(-1);
	for (multi_iterator I(box); !I.end(); I++) {
		const auto i = I.index();
		if (ibox.contains(i)) {
			continue;
		}
		if (amr[i]) {
			const auto d = to_dir(i);
			const auto ip = i + d;
			const auto im = i - d;
			const auto ic = i / 2;
			const auto icp = ic + d;
			const auto icm = ic - d;
			X[i] = (-(3.0 / 32.0) * phi_c[icm] + (15.0 / 16.0) * phi_c[ic] + (5.0 / 32.0) * phi_c[icp]);
			//		X[i] = phi_c[ic];
		}
	}
}

void gravity::initialize_fine(const multi_array<double> &rho, double mtot) {
	for (multi_iterator i(box.pad(-1)); !i.end(); i++) {
		active[i] = true;
	}
	const auto rho0 = mtot;
	for (multi_iterator i(box.pad(-opts.gbw)); !i.end(); i++) {
		R[i] = 4.0 * M_PI * opts.G * (rho[i] - (opts.problem == "sphere" ? 0.0 : rho0));
	}
//	X = phi;
}

void gravity::initialize_coarse(double w) {
	if (!has_phi1) {
		for (multi_iterator i(box); !i.end(); i++) {
			phi[i] = 0.0;
		}
	} else {
		if (!has_phi0) {
			for (multi_iterator i(box); !i.end(); i++) {
				phi[i] = phi1[i];
//				printf( "--- %e\n", phi[i]);
			}
		} else {
			for (multi_iterator i(box); !i.end(); i++) {
				phi[i] = (1.0 - w) * phi0[i] + w * phi1[i];
			}
		}
	}
	for (multi_iterator i(box); !i.end(); i++) {
		X[i] = 0.0;
	}
	for (multi_iterator i(box); !i.end(); i++) {
		active[i] = false;
	}
}

double gravity::coord(index_type i) const {
	return (i + 0.5) * dx;
}

void gravity::unpack_coarse_source(const boundary &bnds) {
	for (int i = 0; i < bnds.data.size(); i++) {
		int k = 0;
		for (multi_iterator j(bnds.boxes[i]); !j.end(); j++) {
			assert(k < bnds.data[i].size());
			R[j.index() * 2] = 4.0 * M_PI * opts.G * bnds.data[i][k];
			k++;
		}
		assert(k == bnds.data[i].size());
	}
}

void gravity::set_outflow_boundaries() {
	for (multi_iterator i(box); !i.end(); i++) {
		if (!box.pad(-opts.gbw).contains(i)) {
			double sum = 0.0;
			double msum = 0.0;
			for (multi_iterator j(box); !j.end(); j++) {
				double r2 = 0.0;
				for (int dim = 0; dim < NDIM; dim++) {
					r2 += std::pow(coord(i.index()[dim]) - coord(j.index()[dim]), 2);
				}
				if (r2 != 0.0) {
					const auto M = std::pow(dx, NDIM) * R[j]/(4.0 * M_PI * opts.G);
					msum += M;
					const auto r = std::sqrt(r2);
					if ( NDIM == 1) {
						sum += opts.G * 4.0 * M_PI * std::abs(r) * M;
					} else if ( NDIM == 2) {
						sum += 2.0 * opts.G * std::log(r) * M;
					} else if ( NDIM == 3) {
						sum -= opts.G * M / r;
					}
				}
			}
			X[i] = sum;
			//		printf( "%e\n", sum);
			active[i] = false;
		} else {
			active[i] = true;
		}
	}
}

void gravity::set_avg_zero() {
	double avg = 0.0;
	for (multi_iterator i(box.pad(-opts.gbw)); !i.end(); i++) {
		avg += phi[i];
	}
	avg /= box.pad(-opts.gbw).volume();
	for (multi_iterator i(box); !i.end(); i++) {
		phi[i] -= avg;
	}
}

std::vector<double> gravity::pack(const multi_range &bbox, int type) const {
	assert(box.contains(bbox));
	std::vector<double> data;
	data.reserve(bbox.volume());
	const auto &Y = type == PACK_X ? X : R;
	for (multi_iterator i(bbox); !i.end(); i++) {
		data.push_back(Y[i]);
	}
	return data;
}

std::vector<double> gravity::pack_amr(const multi_range &bbox, double w) const {
	assert(box.contains(bbox));
	std::vector<double> data;
	data.reserve(bbox.volume());
	for (multi_iterator i(bbox); !i.end(); i++) {
		data.push_back((1.0 - w) * phi0[i] + w * phi1[i]);
	}
	return data;
}

multi_array<double> gravity::get_phi() const {
	multi_array<double> rphi;
	rphi.resize(box);
	for (multi_iterator i(box); !i.end(); i++) {
		rphi[i] = phi[i];
		//	rphi[i] = resid[i];
	}
	return rphi;
}

void gravity::unpack(const std::vector<double> &data, const multi_range &bbox, int type) {
	assert(box.contains(bbox));
	int k = 0;
	auto &Y = type == PACK_X ? X : R;
	for (multi_iterator i(bbox); !i.end(); i++) {
		assert(k < data.size());
		Y[i] = data[k];
		k++;
	}
}

void gravity::finish_fine() {
	has_phi0 = has_phi1;
	has_phi1 = true;
	phi = X;
	phi0 = phi1;
	phi1 = phi;
}

void gravity::relax(bool init_zero) {
	multi_array<double> X0;
	if (init_zero) {
		for (multi_iterator i(box); !i.end(); i++) {
			X[i] = 0.0;
		}
	}
	const auto ibox = box.pad(-1);
	X0 = X;
	auto *x0 = X0.data();
	auto *x1 = X.data();
	const auto s = X.get_strides();
	for (multi_iterator i(ibox); !i.end(); i++) {
		if (active[i]) {
			const auto j = X.index(i);
			x1[j] = 0.0;
			for (int dim = 0; dim < NDIM; dim++) {
				x1[j] += (0.5 / NDIM) * (x0[j + s[dim]] + x0[j - s[dim]]);
			}
			x1[j] = (x1[j] - (0.5 / NDIM) * (R[i] * dx * dx));
			if (i.index()[0] == -2) {
				printf("%e\n", x1[j - s[0]]);
			}
		}
	}
}

gravity_return gravity::get_restrict(double rho0) {
	gravity_return rc;
	multi_array<std::uint8_t> active_c;
	auto rbox = box.pad(-1).half();
	active_c.resize(rbox);
	for (multi_iterator i(rbox); !i.end(); i++) {
		const auto cbox = multi_range(i.index()).double_();
		bool this_active = true;
		for (multi_iterator ci(cbox); !ci.end(); ci++) {
			this_active = this_active && active[ci];
		}
		active_c[i] = this_active;
	}
	for (multi_iterator i(rbox); !i.end(); i++) {
		rc.active.push_back(active_c[i]);
	}
	for (multi_iterator i(box); !i.end(); i++) {
		resid[i] = 0.0;
	}
	const auto *x = X.data();
	const auto s = X.get_strides();
	double rmax = 0.0;
	rc.resid = 0.0;
	rc.mass = 0.0;
	for (multi_iterator i(box.pad(-1)); !i.end(); i++) {
		if (active[i]) {
			const auto j = X.index(i);
			for (int dim = 0; dim < NDIM; dim++) {
				resid[i] -= (x[j + s[dim]] + x[j - s[dim]]) / (dx * dx);
			}
			resid[i] += (2.0 * NDIM) * x[j] / (dx * dx);
			resid[i] += R[i];
			rc.resid += std::abs(resid[i]) * std::pow(dx, NDIM);
			rc.mass += R[i] * std::pow(dx, NDIM) / (4.0 * M_PI * opts.G );
		}
//		phi[i] = resid[i];
	}
	const auto resid_c = resid.restrict_(rbox);
	for (multi_iterator i(rbox); !i.end(); i++) {
		if (active_c[i]) {
			rc.R.push_back(resid_c[i]);
		}
	}
	rc.box = rbox;
	return rc;
}

void gravity::apply_restrict(const gravity_return &data) {
	int k = 0;
	for (multi_iterator i(data.box); !i.end(); i++) {
		assert(k < data.active.size());
		active[i] = data.active[k];
		k++;
	}
	assert(k == data.active.size());
	k = 0;
	for (multi_iterator i(data.box); !i.end(); i++) {
		if (active[i]) {
			assert(k < data.R.size());
			R[i] = data.R[k];
			k++;
		} else {
			R[i] = 0.0;
		}
	}
	assert(k == data.R.size());
}

std::vector<double> gravity::get_prolong(const multi_range &bbox) const {
	std::vector<double> data;
	for (multi_iterator i(bbox); !i.end(); i++) {
		data.push_back(X[i]);
	}
	return data;
}

void gravity::apply_prolong(const std::vector<double> &data) {
	int k = 0;
	const auto cbox = box.pad(-opts.gbw).half().pad(opts.gbw / 2);
	for (multi_iterator i(cbox); !i.end(); i++) {
		assert(k < data.size());
		const auto val = data[k];
		for (multi_iterator j(multi_range(i.index()).double_()); !j.end(); j++) {
			X[j] += val;
		}
		k++;
	}
	assert(k == data.size());
}

