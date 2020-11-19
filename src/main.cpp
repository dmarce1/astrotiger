#include <astrotiger/tree.hpp>
#include <astrotiger/levels.hpp>
#include <astrotiger/output.hpp>
#include <astrotiger/cosmos.hpp>

static tree_client root;

std::vector<double> dx;
std::vector<double> dt;
std::vector<double> last_dt;
std::vector<double> tm;
std::vector<int> super_step;

double solve_gravity(int l, double t, double mtot) {

	const double toler = 5.0e-4;
	int pass = 0;
	double r;
	int oi = 101;
	if (opts.particles) {
		root.compute_cic(std::vector<double>(), t, l).get();
	}
//	printf("Solving gravity on level %i %e\n", l, mtot);
	bool kill = false;
	int iters = 1;
	int n = 0;
	double last_r;
	r = 0.0;
	do {
		auto tmp = root.gravity_solve(pass, l, std::vector<double>(), t, mtot, iters).get();
		last_r = r;
		r = tmp.resid;
//		printf("%3i %3i %e\n", pass, iters, r);
		if (r > 0.5 * last_r && pass > 0) {
			n++;
			iters = std::pow(std::sqrt(2), n);
		}
		if (pass > 500) {
			std::string fname = "X." + std::to_string(oi++) + ".silo";
			output_silo(fname);
//			printf("Gravity solver failed to converged\n");
			//		abort();
			if (pass > 1000) {
				kill = true;
				break;
			}
		}
		pass++;
	} while (r > toler);
	auto tmp = root.gravity_solve(GRAVITY_FINAL_PASS, l, std::vector<double>(), t, mtot, iters).get();
	std::string fname = "X." + std::to_string(oi++) + ".silo";
//	output_silo(fname);
	r = tmp.resid / mtot;
	if (kill) {
		abort();
	}
	return tmp.vmax;
//	printf("%3i %3i %e\n", pass, iters, tmp.resid);
}

bool master(int level, int coarse_level, double tmax, bool already_refined = false) {
	if (level == 0) {
		levels_show();
	}
	cosmos_advance(tm[level]);
	double coarse_a0 = cosmos_a();
	if (level > opts.max_level) {
//		printf( "**Applying coarse correction to level %i\n", level - 1);
		levels_apply_coarse_correction(level - 1, coarse_a0, coarse_a0);
		return false;
	}
//	int oi = 0;
	double nstep = -1;
	const int refine_freq = 1;
	statistics stats;
	if (opts.self_gravity) {
		stats = root.get_statistics(opts.max_level, tm[level]).get();
	}
	const int max_refined = stats.min_level;
	int this_step = 0;
	int oi = 100;
	static auto e = root.get_energy_statistics().get();
	static auto last_e = e;
	cosmos_advance(0.0);
	static double last_a = cosmos_a();
	static double a = cosmos_a();
	static double epec = 0.0;
	do {
		if (level == 0) {
			cosmos_advance(tm[level]);
			last_e = e;
			last_a = a;
			e = root.get_energy_statistics().get();
			a = cosmos_a();
			epec += 0.5 * (e.ekin + last_e.ekin) * (a - last_a);
			const auto etot = a * e.ekin + a * e.epot + epec;
			static const auto etot0 = etot;
			printf("%e %e %e %e %e %e %e\n", tm[level], a * e.ekin, a * e.epot, epec, etot, a, (etot-etot0)/e.ekin);
		}
		const auto a = cosmos_a();
		bool refine = ((nstep == -1) && !already_refined) || (this_step % 2 == 0 && this_step > 0 && level == 0);
		if (nstep != -1) {
			coarse_level = level;
		}
		last_dt[level] = dt[level];
		cosmos_advance(tm[level]);
		if (refine) {
			for (int l = level; l < opts.max_level; l++) {
//				printf("Refining level %i\n", l);
				levels_hydro_initialize(l, refine);
//				levels_show();
				levels_set_child_families(l);
				levels_get_hydro_boundaries(l + 1, tm[l + 1]);
//				levels_energy_update(l + 1);
			}
		} else {
			levels_hydro_initialize(level, false);
		}
		double amax;
		amax = levels_fine_fluxes(level);
		if (amax == 0.0) {
			tm[level] = tm[level - 1];
			master(level + 1, coarse_level, tm[level]);
			double coarse_a1 = cosmos_a();
//			printf( "--Applying coarse correction to level %i\n", level - 1);
			levels_apply_coarse_correction(level - 1, coarse_a0, coarse_a1);
			return false;
		}
//		levels_show();
		double mtot1, mtot2;
		if (opts.self_gravity) {
			auto tmp = root.get_statistics(std::min(level, max_refined), tm[level]).get();
			assert(tmp.u.size());
			const auto mtot = tmp.u[rho_i];
			mtot1 = mtot;
//			printf("max_refined = %i level = %i\n", max_refined, level, mtot);
			const auto gmax = solve_gravity(level, tm[level], mtot);
//			printf( "gmax  = %e\n", gmax);
			amax = std::max(amax, gmax);
//			amax = std::max(amax, opts.cfl * a * dx[level] * cosmos_adot() / cosmos_a());
		}
		dt[level] = opts.cfl * a * dx[level] / amax;
		nstep = std::ceil((tmax - tm[level]) / dt[level]);
		dt[level] = ((tmax - tm[level]) / nstep);
		dt[level] = std::min(dt[level], tmax - tm[level]);
		const auto a1 = cosmos_a();
		const auto H1 = cosmos_adot();
		cosmos_advance(tm[level] + dt[level]);
		const auto a2 = cosmos_a();
		const auto H2 = cosmos_adot();
		cosmos_advance(tm[level]);
//		printf("Advancing level %i from %e to %e scale factor %e to %e %e %e\n", level, tm[level], tm[level] + dt[level], a1, a2, H1, H2);
		if (opts.particles) {
			if (level == levels_max_level()) {
				root.kick(coarse_level, tm[level], last_dt, dt).get();
			}
		}
		levels_hydro_substep(level, 0, dt[level], false, a1, a2);
		cosmos_advance(tm[level] + dt[level]);
		if (opts.self_gravity) {
			const auto mtot = root.get_statistics(std::min(level, max_refined), tm[level] + dt[level]).get().u[rho_i];
//			printf("max_refined = %i level = %i\n", max_refined, level, mtot);
			mtot2 = mtot;
			solve_gravity(level, tm[level] + dt[level], mtot);
		}
		levels_hydro_substep(level, 1, dt[level], nstep == 1.0, a1, a2);
		tm[level] += dt[level];
		const bool has_next_level = master(level + 1, coarse_level, tm[level], refine);
		///	printf( "-\n");
		//	std::string fname = "X." + std::to_string(oi++) + ".silo";
		//	output_silo(fname);
		if (opts.particles && !has_next_level) {
			//e		printf("%e %e %e %i %e %e %e %e %e %e\n", amax, tm[level], dt[level], coarse_level, a1, a2, H1, H2, mtot1, mtot2);
			cosmos_advance(tm[level] + 0.5 * dt[level]);
			root.drift(dt[level]).get();
			root.finish_drift(std::vector<particle>()).get();
		}
		this_step++;
	} while (nstep != 1.0);
	cosmos_advance(tm[level]);
	double coarse_a1 = cosmos_a();
	if (level > 0) {
//		printf( "Applying coarse correction to level %i\n", level - 1);
		levels_apply_coarse_correction(level - 1, coarse_a0, coarse_a1);
	}
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
	double dt = opts.tmax / 100.0;
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
