#include <astrotiger/gravity.hpp>
#include <astrotiger/options.hpp>

void gravity::resize(double dx_, const multi_range &box_) {
	dx = dx_;
	const auto cbox = box_.half().pad(opts.gbw);
	box = box_.pad(opts.gbw);
	phi.resize(box);
	X.resize(box);
	R.resize(box);
	active.resize(box);
	amr.resize(box);
	phi_c.resize(cbox);
	refined.resize(box);
	for (int dim = 0; dim < NDIM; dim++) {
		fbox[dim] = box_;
		fbox[dim].max[dim]++;
		flux[dim].resize(fbox[dim]);
	}
	for (int dim = 0; dim < NDIM; dim++) {
		dir[dim] = 0;
		dir[dim][dim] = 1;
	}
	for (multi_iterator i(box); !i.end(); i++) {
		phi[i] = 0.0;
		refined[i] = false;
	}
}

void gravity::compute_flux() {
	for (int dim = 0; dim < NDIM; dim++) {
		for (multi_iterator i(fbox[dim]); !i.end(); i++) {
			const auto im = i.index() - dir[dim];
			if (!refined[im] && !refined[i]) {
				flux[dim][i] = (X[i] - X[im]) / dx;
			}
		}
	}
}

void gravity::set_refined(const std::vector<multi_range> &boxes) {
	for (multi_iterator i(box); !i.end(); i++) {
		refined[i] = false;
	}
	if (boxes.size()) {
		for (const auto &b : boxes) {
			const auto cbox = b.half();
			for (multi_iterator i(box); !i.end(); i++) {
				refined[i] = refined[i] || cbox.contains(i);
			}
		}
		for (multi_iterator i(box); !i.end(); i++) {
			active[i] = active[i] && !refined[i];
		}
	}
}

void gravity::set_amr_zones(const std::vector<multi_range> &boxes, const std::vector<double> &data) {
	const auto cbox = box.pad(-opts.gbw).half().pad(opts.gbw);
	int k = 0;
	for (multi_iterator i(cbox); !i.end(); i++) {
		assert(k < data.size());
		phi_c[i] = data[k];
		k++;
	}
	for (multi_iterator i(box); !i.end(); i++) {
		amr[i] = false;
		for (const auto &b : boxes) {
			amr[i] = b.contains(i);
			if (amr[i]) {
				break;
			}
		}
	}
}
void gravity::compute_amr_bounds(bool plus_interior) {
	const auto cbox = box.pad(-opts.gbw).half().pad(opts.gbw);
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
	const auto ibox = box.pad(-opts.gbw);
	if (plus_interior) {
		for (multi_iterator I(ibox); !I.end(); I++) {
			const auto i = I.index();
			const auto ibox = box.pad(-opts.gbw);
			const auto d = to_dir(i);
			const auto ic = i / 2;
			const auto icp = ic + d;
			const auto icm = ic - d;
			X[i] = (-(3.0 / 32.0) * phi_c[icm] + (15.0 / 16.0) * phi_c[ic] + (5.0 / 32.0) * phi_c[icp]);
		}
	}
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
			if (box.contains(ip) && !amr[ip]) {
				X[i] = (-(1.0 / 14.0) * phi_c[icm] + (5.0 / 6.0) * phi_c[ic] + (5.0 / 21.0) * X[ip]);
			} else {
				X[i] = (-(3.0 / 32.0) * phi_c[icm] + (15.0 / 16.0) * phi_c[ic] + (5.0 / 32.0) * phi_c[icp]);
			}
		}
	}
}

void gravity::initialize_fine(const multi_array<double> &rho, double mtot) {
	for (multi_iterator i(box.pad(-1)); !i.end(); i++) {
		active[i] = true;
	}
	const auto rho0 = mtot;
	for (multi_iterator i(box); !i.end(); i++) {
		R[i] = 4.0 * M_PI * opts.G * (rho[i] - rho0);
	}
	X = phi;
}

void gravity::initialize_coarse(double w) {
	for (multi_iterator i(box); !i.end(); i++) {
		X[i] = 0.0;
	}
	for (multi_iterator i(box); !i.end(); i++) {
		active[i] = false;
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
	const auto &Y = type == PACK_PHI ? phi : X;
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
		data.push_back(phi[i]);
	}
	return data;
}

multi_array<double> gravity::get_phi() const {
	multi_array<double> rphi;
	auto bbox = box.pad(-opts.gbw + 1);
	rphi.resize(bbox);
	for (multi_iterator i(bbox); !i.end(); i++) {
		rphi[i] = phi[i];
		rphi[i] = X[i];
	}
	return rphi;
}

void gravity::unpack(const std::vector<double> &data, const multi_range &bbox) {
	assert(box.contains(bbox));
	int k = 0;
	for (multi_iterator i(bbox); !i.end(); i++) {
		assert(k < data.size());
		X[i] = data[k];
		k++;
	}
}

void gravity::finish_fine() {
	phi = X;
}

void gravity::relax(bool init_zero) {
	multi_array<double> X0;
	compute_amr_bounds(false);
	compute_flux();
	if (init_zero) {
		for (multi_iterator i(box); !i.end(); i++) {
			X[i] = 0.0;
		}
	}
	const auto ibox = box.pad(-opts.gbw);
	for (multi_iterator i(ibox); !i.end(); i++) {
		if (active[i]) {
			double resid = R[i];
			for (int dim = 0; dim < NDIM; dim++) {
				resid -= (flux[dim][i.index() + dir[dim]] - flux[dim][i]) / dx;
			}
			X[i] -= (resid / (2 * NDIM)) * dx * dx;
		}
	}
}

std::vector<double> gravity::get_flux_restrict() const {
	std::vector<double> data;
	for (int dim = 0; dim < NDIM; dim++) {
		auto cfbox = box.pad(-opts.gbw).half();
		cfbox.max[dim]++;
		multi_array<double> this_flux;
		this_flux.resize(cfbox);
		for (multi_iterator i(cfbox); !i.end(); i++) {
			this_flux[i] = 0.0;
		}
		for (multi_iterator i(fbox[dim]); !i.end(); i++) {
			if (i.index()[dim] % 2 == 0) {
				const auto ic = i.index() / 2;
				const auto f = (phi[i] - phi[i.index() - dir[dim]]) / dx;
				this_flux[ic] += f / (1 << (NDIM - 1));
			}
		}
		for (multi_iterator i(cfbox); !i.end(); i++) {
			data.push_back(this_flux[i]);
		}
	}
	return std::move(data);
}

void gravity::apply_flux_restrict(const std::vector<double> &data, const multi_range &bbox) {
	int k = 0;
	for (int dim = 0; dim < NDIM; dim++) {
		auto this_box = bbox;
		this_box.max[dim]++;
		for (multi_iterator i(this_box); !i.end(); i++) {
			assert(k < data.size());
			flux[dim][i] = data[k];
			k++;
		}
	}
	assert(k == data.size());
}

gravity_return gravity::get_restrict(double rho0) {
	compute_amr_bounds(false);
	compute_flux();
	gravity_return rc;
	multi_array<std::uint8_t> active_c;
	auto rbox = box.pad(-opts.gbw).half();
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
	multi_array<double> resid(box);
	for (multi_iterator i(box); !i.end(); i++) {
		resid[i] = 0.0;
	}
	const auto *x = X.data();
	const auto s = X.get_strides();
	double rmax = 0.0;
	for (multi_iterator i(box.pad(-opts.gbw)); !i.end(); i++) {
		if (active[i]) {
			double res = R[i];
			for (int dim = 0; dim < NDIM; dim++) {
				res -= (flux[dim][i.index() + dir[dim]] - flux[dim][i]) / dx;
			}
			resid[i] = res;
			rmax = std::max(rmax, std::abs(resid[i] / (4.0 * M_PI * opts.G) / rho0));
		}
	}
	const auto resid_c = resid.restrict_(rbox);
	for (multi_iterator i(rbox); !i.end(); i++) {
		if (active_c[i]) {
			rc.R.push_back(resid_c[i]);
		}
	}
	rc.resid = rmax;
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
	const auto cbox = box.pad(-opts.gbw).half();
	for (multi_iterator i(cbox); !i.end(); i++) {
		assert(k < data.size());
		const auto val = data[k];
		for (multi_iterator j(multi_range(i.index()).double_()); !j.end(); j++) {
//			X[j] += val;
		}
		k++;
	}
	assert(k == data.size());
}

