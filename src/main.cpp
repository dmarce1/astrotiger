#include <astrotiger/tree.hpp>
#include <astrotiger/levels.hpp>
#include <astrotiger/output.hpp>

static tree_client root;

std::vector<double> dx;
std::vector<double> dt;
std::vector<double> last_dt;
std::vector<double> tm;
std::vector<int> super_step;

void solve_gravity(int l, double t, double mtot) {

	const double toler = 5.0e-4;
	int pass = 0;
	double r;
	int oi = 101;
	if (opts.particles) {
		root.compute_cic(std::vector<double>(), t, l).get();
	}
	printf("Solving gravity on level %i %e\n", l, mtot);
	bool kill = false;
	do {
		auto tmp = root.gravity_solve(pass, l, std::vector<double>(), t, mtot).get();
		r = tmp.resid;
		printf("%i %e\n", pass, r);
		if (pass > 500) {
			std::string fname = "X." + std::to_string(oi++) + ".silo";
			output_silo(fname);
			printf("Gravity solver failed to converged\n");
			//		abort();
			if (pass > 1000) {
				kill = true;
				break;
			}
		}
		pass++;
	} while (r > toler);
	auto tmp = root.gravity_solve(GRAVITY_FINAL_PASS, l, std::vector<double>(), t, mtot).get();
	std::string fname = "X." + std::to_string(oi++) + ".silo";
//	output_silo(fname);
	r = tmp.resid / mtot;
	if (kill) {
		abort();
	}
	//printf("%i %e\n", pass, tmp.resid);
}

bool master(int level, int coarse_level, double tmax) {
	if (level > opts.max_level) {
		return false;
	}
	int oi = 0;
	double nstep = -1;
	const int refine_freq = 1;
	statistics stats;
	if (opts.self_gravity) {
		stats = root.get_statistics(opts.max_level, tm[level]).get();
	}
	const int max_refined = stats.min_level;
	do {
//		if (level > 0 && nstep == -1) {
//			levels_set_child_families(level - 1);
//		}
		bool refine = (nstep == -1);
		if (refine) {
//			printf("Refining on level %i\n", level);
		}
		if (nstep != -1) {
			coarse_level = level;
		}
//		printf("%i\n", coarse_level);
//		std::string fname = "X." + std::to_string(oi++) + ".silo";
//		output_silo(fname);

//		printf("Hydro pre-step\n");
		last_dt[level] = dt[level];
		dt[level] = levels_hydro_initialize(level, refine);
		if (refine) {
			levels_show();
		}
		if (refine) {
			for (int l = level; l <= opts.max_level; l++) {
				levels_set_child_families(l);
			}
		}
		dt[level] = std::max(dt[level], levels_fine_fluxes(level));
		if (dt[level] == 0.0) {
			tm[level] = tm[level - 1];
			master(level + 1, coarse_level, tm[level]);
			return false;
		}
//		levels_show();
		dt[level] = opts.cfl * dx[level] / dt[level];
		nstep = std::ceil((tmax - tm[level]) / dt[level]);
		dt[level] = (tmax - tm[level]) / nstep;
		printf("Advancing level %i from %e to %e\n", level, tm[level], tm[level] + dt[level]);
		if (opts.self_gravity) {
			auto tmp = root.get_statistics(std::min(level, max_refined), tm[level]).get();
			assert(tmp.u.size());
			const auto mtot = tmp.u[rho_i];
//			printf("max_refined = %i level = %i\n", max_refined, level, mtot);
			solve_gravity(level, tm[level], mtot);
		}
		if (level == opts.max_level) {
			root.kick(coarse_level, tm[level], last_dt, dt).get();
		}
		levels_hydro_substep(level, 0, dt[level], false);
		if (opts.self_gravity) {
			const auto mtot = root.get_statistics(std::min(level, max_refined), tm[level] + dt[level]).get().u[rho_i];
//			printf("max_refined = %i level = %i\n", max_refined, level, mtot);
			solve_gravity(level, tm[level] + dt[level], mtot);
		}
		levels_hydro_substep(level, 1, dt[level], nstep == 1.0);
		tm[level] += dt[level];
		const bool has_next_level = master(level + 1, coarse_level, tm[level]);
		///	printf( "-\n");
		if (opts.particles && !has_next_level) {
			//		printf( "%i\n", coarse_level);
			root.drift(dt[level]).get();
			root.finish_drift(std::vector<particle>()).get();
		}
	} while (nstep != 1.0);
	levels_apply_coarse_correction(level);
	super_step[level]++;
	return true;
}

int hpx_main(int argc, char *argv[]) {
//	srand(1000);
	options opts;
	opts.process_options(argc, argv);
	levels_init();
	multi_range box;
	for (int i = 0; i < NDIM; i++) {
		box.min[i] = 0;
		box.max[i] = opts.max_box;
	}
	std::vector<sibling> sibs;
	root = tree::allocate(0, box).get();
	multi_range periodic;
	for (int dim = 0; dim < NDIM; dim++) {
		periodic.min[dim] = -1;
		periodic.max[dim] = 2;
	}
	for (multi_iterator i(periodic); !i.end(); i++) {
		bool zero = true;
		for (int dim = 0; dim < NDIM; dim++) {
			zero = zero && (i[dim] == 0);
		}
		if (!zero) {
			bool use_periodic = true;
			for (int dim = 0; dim < NDIM; dim++) {
				if ((i.index()[dim] == -1 && opts.bnd[2 * dim] != PERIODIC) || (i.index()[dim] == +1 && opts.bnd[2 * dim + 1] != PERIODIC)) {
					use_periodic = false;
				}
			}
			if (use_periodic) {
				const auto shift = i.index() * opts.max_box;
				sibling sib;
				sib.client = root;
				sib.shift = shift;
				sibs.push_back(sib);
			}
		}
	}
	dx.resize(opts.max_level + 1);
	dt.resize(opts.max_level + 1);
	last_dt.resize(opts.max_level + 1);
	tm.resize(opts.max_level + 1);
	super_step.resize(opts.max_level + 1);
	for (int l = 0; l <= opts.max_level; l++) {
		printf("Forming level %i\n", l);
		super_step[l] = 0;
		dx[l] = 1.0 / (opts.max_box * (1 << l));
		tm[l] = 0.0;
		auto fut = root.initialize(l);
		auto amax = fut.get();
		printf("Initialized. Setting families\n");
		if (l >= 1) {
			levels_set_child_families(l - 1);
		} else {
			root.set_family(tree_client(), root, sibs).get();
		}
//		if (amax != 0.0) {
//			dt[l] = opts.cfl * dx[l] / amax;
//		} else {
		dt[l] = 0.0;
//		}
		printf("Done\n");
	}
	root.restrict_all().get();
//	master(0, 1.0e-100);
	levels_show();
	if (opts.self_gravity) {
		const auto tmp = root.get_statistics(opts.max_level, 0.0).get();
		const auto max_refined = tmp.min_level;
		for (int l = 0; l <= opts.max_level; l++) {
			const auto tmp = root.get_statistics(std::min(l, max_refined), 0.0).get();
			assert(tmp.u.size());
			const auto mtot = tmp.u[rho_i];
			solve_gravity(l, 0.0, mtot);
		}
	}
	output_silo("X.0.silo");
	int i = 0;
	const auto dt = 0.01;
	levels_show();
	for (double t = 0.0; t < opts.tmax; t += dt) {
		i++;
		master(0, 0, std::min(t + dt, opts.tmax));
		std::string fname = "X." + std::to_string(i) + ".silo";
		output_silo(fname);
	}
	std::string fname = "X." + std::to_string(i) + ".silo";
	output_silo(fname);
	return hpx::finalize();
}

#ifndef HPX_LITE
int main(int argc, char *argv[]) {

	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };

	hpx::init(argc, argv, cfg);
	printf("exiting..\n");
}
#endif
