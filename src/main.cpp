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

double solve_gravity(int l, double t, double mtot, int max_refined) {

	const double toler = 1.0e-3;
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
			iters = std::pow(2, n);
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
	if (l <= max_refined) {
		const auto dif = -root.get_average_phi(l).get();
		root.set_average_phi(l, dif).get();
	}
	return tmp.vmax;
//	printf("%3i %3i %e\n", pass, iters, tmp.resid);
}

bool master(int level, int coarse_level, double tmax, bool already_refined = false) {
	if (level == 0) {
		levels_show();
	}
//	cosmos_set_z(opts.z, opts.H0);
//	cosmos_advance(tm[level]);
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
	double mtot;
	if (opts.self_gravity) {
		stats = root.get_statistics(opts.max_level, tm[level]).get();
		mtot = stats.u[rho_i];
	} else {
		mtot = 0.0;
	}
	const int max_refined = stats.min_level;
	int this_step = 0;
	int oi = 100;
	static auto e = root.get_energy_statistics(mtot).get();
	static auto last_e = e;
	static bool first = true;
	if (first) {
		cosmos_advance(0.0);
		first = false;
	}
	static double last_a = cosmos_a();
	static double a = cosmos_a();
	static double epec = 0.0;
	static int base_step = 0;
	do {
		cosmos_advance(tm[level]);
		levels_energy_update(level);
		//	const auto a = cosmos_a();
		bool refine = ((nstep == -1) && !already_refined) || (this_step % 2 == 0 && this_step > 0 && level == 0);
//		refine = false;
//		printf( "%i\n", refine);
		if (nstep != -1) {
			coarse_level = level;
		}
		last_dt[level] = dt[level];
		double amax;
		if (refine) {
			for (int l = level; l <= opts.max_level; l++) {
//				printf("Refining level %i\n", l);
				double this_amax = levels_hydro_initialize(l, refine, nstep == -1);
				if (l == level) {
					amax = this_amax;
				}
				//			levels_show();
				levels_set_child_families(l);
				if (l < opts.max_level) {
					levels_get_hydro_boundaries(l + 1, tm[l + 1]);
				}
				auto tmp = root.get_statistics(std::min(level, max_refined), tm[level]).get();
				const auto mtot = tmp.u[rho_i];
				solve_gravity(l, tm[l], mtot, max_refined);
			}
		} else {
			amax = levels_hydro_initialize(level, false, nstep == -1);
		}
//		printf("First step %i\n", nstep == -1);
		levels_get_hydro_boundaries(level, tm[level]);
		amax = std::max(amax, levels_compute_fluxes(level, 0));
		if (level < opts.max_level) {
			levels_get_hydro_boundaries(level + 1, tm[level]);
			amax = std::max(amax, levels_compute_fluxes(level + 1, 0));
		}
		amax = levels_fine_fluxes(level);
		auto vmax = root.max_part_velocity(level).get();
		auto cmax = 10 * opts.cfl * dx[level] * std::abs(cosmos_adot());
		//	printf("%e\n", a / std::abs(cosmos_adot()));
		//	printf("%e\n", amax);
		double mtot1, mtot2;
		double gmax;
		if (opts.self_gravity) {
			auto tmp = root.get_statistics(std::min(level, max_refined), tm[level]).get();
			assert(tmp.u.size());
			const auto mtot = tmp.u[rho_i];
			mtot1 = mtot;
//			printf("max_refined = %i level = %i\n", max_refined, level, mtot);
			gmax = solve_gravity(level, tm[level], mtot, max_refined);
//			printf( "amax = %e gmax  = %e\n", amax, gmax);
			gmax;
		}
		auto hmax = amax;
		amax = std::max(amax, gmax);
		amax = std::max(amax, cmax);
		amax = std::max(amax, vmax);
		if (level > levels_max_level()) {
			tm[level] = tm[level - 1];
//			master(level + 1, coarse_level, tm[level]);
			double coarse_a1 = cosmos_a();
//			printf( "--Applying coarse correction to level %i\n", level - 1);
			levels_apply_coarse_correction(level - 1, coarse_a0, coarse_a1);
			return false;
		}
//		levels_show();
		dt[level] = opts.cfl * a * dx[level] / amax;
//		printf( "%e %e %e", tmax, tm[level], dt[level]);
		nstep = std::ceil((tmax - tm[level]) / dt[level]);
		dt[level] = ((tmax - tm[level]) / nstep);
		dt[level] = std::min(dt[level], tmax - tm[level]);
		if (level == 0) {
			last_e = e;
			last_a = a;
			stats = root.get_statistics(max_refined, tm[level]).get();
			mtot = stats.u[rho_i];
			e = root.get_energy_statistics(mtot).get();
			a = cosmos_a();
			epec += 0.5 * (e.ekin + last_e.ekin) * (a - last_a);
			const auto etot = a * e.ekin + a * e.epot + epec;
			static const auto etot0 = etot;
			const auto adot = cosmos_adot();
			const auto adotdot = cosmos_adotdot();
			printf(
					"\n step      time         dt            ekin          epot          epec          etot         a          adot        adotdot       econ      mtot");
			printf("\n%5i  %e %e %e %e %e %e %e %e %e %e %e\n", base_step + this_step, tm[level], dt[level], a * e.ekin, a * e.epot, epec, etot, a, adot,
					adotdot, (etot - etot0) / (e.ekin + epec), mtot);
		}
		char limit;
		if (amax == hmax) {
			limit = 'H';
		} else if (amax == cmax) {
			limit = 'C';
		} else if (amax == vmax) {
			limit = 'V';
		} else if (amax == gmax) {
			limit = 'G';
		}
		if (opts.particles) {
			printf("\t Hydro level %i; ", level);
			if (level == levels_max_level()) {
				printf(" Kicking level %i+; Drifting all levels; ", coarse_level);
			} else {
				printf("                                        ");
			}
			printf("time from %e to %e %c\n", tm[level], tm[level] + dt[level], limit);
		} else {
			printf("\t Advancing level %i from %e to %e.\n", level, tm[level], tm[level] + dt[level]);
		}
		if (opts.particles) {
			if (level == levels_max_level()) {
				root.kick(coarse_level, tm[level], last_dt, dt).get();
			}
		}
		levels_hydro_substep(level, 0, dt[level], false, a, a);
		const auto a0 = cosmos_a();
		cosmos_advance(tm[level] + dt[level]);
		const auto a1 = cosmos_a();
		if (opts.self_gravity) {
			const auto mtot = root.get_statistics(std::min(level, max_refined), tm[level] + dt[level]).get().u[rho_i];
//			printf("max_refined = %i level = %i\n", max_refined, level, mtot);
			mtot2 = mtot;
			solve_gravity(level, tm[level] + dt[level], mtot, max_refined);
		}
		levels_hydro_substep(level, 1, dt[level], nstep == 1.0, a, a);
		tm[level] += dt[level];
		const bool has_next_level = master(level + 1, coarse_level, tm[level], refine);
		///	printf( "-\n");
		//	std::string fname = "X." + std::to_string(oi++) + ".silo";
		//	output_silo(fname);
		if (opts.particles && !has_next_level) {
			//e		printf("%e %e %e %i %e %e %e %e %e %e\n", amax, tm[level], dt[level], coarse_level, a1, a2, H1, H2, mtot1, mtot2);
			root.drift(dt[level], a0, a1).get();
			root.finish_drift(std::vector<particle>()).get();
		}
		this_step++;
		cosmos_advance(tm[level]);
	} while (nstep != 1.0);
	double coarse_a1 = cosmos_a();
	if (level > 0) {
		levels_apply_coarse_correction(level - 1, coarse_a0, coarse_a1);
	}
	super_step[level]++;
	if (level == 0) {
		base_step += this_step;
	}
	return true;
}

void chemistry_test();

int hpx_main(int argc, char *argv[]) {
//	srand(1000);
//	chemistry_test();
//	return hpx::finalize();

	options opts;

	opts.process_options(argc, argv);

//	for (double t = 0.0; t < 125; t += 0.1 * std::abs(cosmos_a() / cosmos_adot())) {
//		cosmos_advance(t);
//		printf("%e %e %e %e\n", t, cosmos_a(), cosmos_adot(), cosmos_adotdot());
//		if (cosmos_a() < 1.0 / 129.0 && t > 0.0) {
//			break;
//		}
//	}
//	return hpx::finalize();
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
			solve_gravity(l, 0.0, mtot, max_refined);
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
