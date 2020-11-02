#include <astrotiger/gravity.hpp>
#include <astrotiger/options.hpp>

constexpr double omega = 1.9;

void gravity::resize(double dx_, const multi_range &box_) {
	inner_box = box_;
	dx = dx_;
	box = box_.pad(opts.gbw);
	phi0.resize(box);
	phi.resize(box);
	X.resize(box);
	R.resize(box);
	active.resize(box);
#ifdef OUTPUT_RESID
	resid.resize(box);
#endif

}

void gravity::allocate(const multi_range &box_) {
	box = box_;
	X.resize(box);
	R.resize(box);
	active.resize(box);
#ifdef OUTPUT_RESID
	resid.resize(box);
#endif

}

void gravity::set_amr_zones(const std::vector<multi_range> &boxes, const std::vector<double> &data) {
	const auto cbox = inner_box.half().pad(opts.gbw / 2 + opts.gbw % 2 + 1);
	int k = 0;
	multi_array<double> phi_c(cbox);
	for (multi_iterator i(cbox); !i.end(); i++) {
		assert(k < data.size());
		phi_c[i] = data[k];
		k++;
	}
	assert(k==data.size());

	for (multi_iterator i(box); !i.end(); i++) {
		if (!inner_box.contains(i)) {
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
//	for (multi_iterator i(box); !i.end(); i++) {
//		amr[i] = false;
//		for (const auto &b : boxes2) {
//			amr[i] = b.contains(i);
//			if (amr[i]) {
//				break;
//			}
//		}
//	}
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
	for (multi_iterator I(box); !I.end(); I++) {
		const auto i = I.index();
		const auto d = to_dir(i);
		const auto ic = i / 2;
		const auto icp = ic + d;
		const auto icm = ic - d;
		X[i] = (-(3.0 / 32.0) * phi_c[icm] + (15.0 / 16.0) * phi_c[ic] + (5.0 / 32.0) * phi_c[icp]);
	}
}

void gravity::initialize_fine(const multi_array<double> &rho, double mtot) {
	for (multi_iterator i(box.pad(-1)); !i.end(); i++) {
		active[i] = true;
	}
	const auto rho0 = mtot;
	for (multi_iterator i(inner_box); !i.end(); i++) {
		R[i] = 4.0 * M_PI * opts.G * (rho[i] - (opts.problem == "sphere" ? 0.0 : rho0));
	}
//	X = phi;
}

void gravity::initialize_coarse(double w) {
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
			const auto l = j.index() * 2;
			for (multi_iterator l(multi_range(j.index()).double_()); !l.end(); l++) {
				if (box.contains(l)) {
					R[l] = 4.0 * M_PI * opts.G * bnds.data[i][k];
				}
			}
			k++;
		}
		assert(k == bnds.data[i].size());
	}
}

void gravity::set_outflow_boundaries() {
	for (multi_iterator i(box); !i.end(); i++) {
		if (!box.pad(-1).contains(i)) {
			double sum = 0.0;
			double msum = 0.0;
			for (multi_iterator j(box); !j.end(); j++) {
				double r2 = 0.0;
				for (int dim = 0; dim < NDIM; dim++) {
					r2 += std::pow(coord(i.index()[dim]) - coord(j.index()[dim]), 2);
				}
				if (r2 != 0.0) {
					const auto M = std::pow(dx, NDIM) * R[j] / (4.0 * M_PI * opts.G);
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
	for (multi_iterator i(inner_box); !i.end(); i++) {
		avg += phi[i];
	}
	avg /= inner_box.volume();
	for (multi_iterator i(box); !i.end(); i++) {
		phi[i] -= avg;
	}
}

std::vector<double> gravity::pack(const multi_range &bbox) const {
	assert(box.contains(bbox));
	std::vector<double> data;
	data.reserve(bbox.volume());
	for (multi_iterator i(bbox); !i.end(); i++) {
		data.push_back(X[i]);
	}
	return data;
}

std::vector<double> gravity::pack_phi(const multi_range &bbox) const {
	assert(box.contains(bbox));
	std::vector<double> data;
	data.reserve(bbox.volume());
	for (multi_iterator i(bbox); !i.end(); i++) {
		data.push_back(phi[i]);
	}
	return data;
}

void gravity::to_array(multi_array<double> &a, const multi_range &bbox, double w) const {
	for (multi_iterator i(bbox); !i.end(); i++) {
		if (box.contains(i)) {
			a[i] = w == 0.0 ? phi[i] : w * phi[i] + (1.0 - w) * phi0[i];
		}
	}
}

multi_array<double> gravity::get_phi() const {
	multi_array<double> rphi;
	rphi.resize(box);
	for (multi_iterator i(box); !i.end(); i++) {
#ifdef OUTPUT_RESID
			rphi[i] = resid[i];
#else
		rphi[i] = phi[i];
#endif
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
	phi0 = phi;
	phi = X;
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
			int sum = 0;
			for (int dim = 0; dim < NDIM; dim++) {
				sum += i.index()[dim];
			}
			if (sum % 2 == red % 2) {
				double r = 0.0;
				for (int dim = 0; dim < NDIM; dim++) {
					r += (0.5 / NDIM) * (x1[j + s[dim]] + x1[j - s[dim]]);
				}
				r += -x1[j] - R[i] * dx * dx / (2.0 * NDIM);
				x1[j] += omega * r;
			}
//			if (i.index()[0] == -2) {
//				printf("%e\n", x1[j - s[0]]);
//			}
		}
	}
	red++;
}

gravity_return gravity::get_restrict(double rho0) {
	gravity_return rc;
	multi_array<std::uint8_t> active_c;
	const auto bigbox = inner_box.pad(opts.gbw - 1);
	auto rbox = bigbox.half();
	active_c.resize(rbox);
	for (multi_iterator i(rbox); !i.end(); i++) {
		const auto id = i.index();
		if (box.pad(-1).contains(id * 2)) {
			const auto cbox = multi_range(id).double_();
			bool this_active = true;
			for (multi_iterator ci(cbox); !ci.end(); ci++) {
				this_active = this_active && active[ci];
			}
			active_c[i] = this_active;
		} else {
			active_c[i] = false;
		}
	}
	for (multi_iterator i(rbox); !i.end(); i++) {
		rc.active.push_back(active_c[i]);
	}
#ifndef OUTPUT_RESID
	multi_array<double> resid;
#endif
	resid.resize(inner_box.pad(opts.gbw));
	for (multi_iterator i(bigbox); !i.end(); i++) {
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
			rc.mass += R[i] * std::pow(dx, NDIM) / (4.0 * M_PI * opts.G);
		}
//		phi[i] = resid[i];
	}
	const auto resid_c = resid.restrict_(rbox);
	for (multi_iterator i(rbox); !i.end(); i++) {
		if (box.contains(i)) {
			rc.R.push_back(resid_c[i]);
		} else {
			rc.R.push_back(0.0);
		}
	}
	rc.box = rbox;
	return rc;
}

void gravity::apply_restrict(const gravity_return &data) {
	int k = 0;
	for (multi_iterator i(data.box); !i.end(); i++) {
		assert(k < data.active.size());
		if (box.contains(i)) {
			active[i] = data.active[k];
		}
		k++;
	}
	assert(k == data.active.size());
	k = 0;
	for (multi_iterator i(data.box); !i.end(); i++) {
		assert(k < data.R.size());
		if (box.contains(i)) {
			R[i] = data.R[k];
		}
		k++;
	}
	assert(k == data.R.size());
}

std::vector<double> gravity::get_prolong(const multi_range &bbox) const {
	std::vector<double> data;
	for (multi_iterator i(bbox); !i.end(); i++) {
		if (box.contains(i)) {
			data.push_back(X[i]);
		} else {
			data.push_back(0.0);
		}
	}
	return data;
}

void gravity::apply_prolong(const std::vector<double> &data) {
	int k = 0;
	const auto cbox = inner_box.half().pad(opts.gbw - 1);
	for (multi_iterator i(cbox); !i.end(); i++) {
		assert(k < data.size());
		const auto val = data[k];
		for (multi_iterator j(multi_range(i.index()).double_()); !j.end(); j++) {
			if (box.pad(-1).contains(j)) {
				X[j] += val;
			}
		}
		k++;
	}
	assert(k == data.size());
}

