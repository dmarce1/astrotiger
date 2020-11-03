#include <astrotiger/tree.hpp>
#include <astrotiger/levels.hpp>
#include <astrotiger/output.hpp>

static tree_client root;

std::vector<double> dx;
std::vector<double> dt;
std::vector<double> tm;
std::vector<int> super_step;

statistics solve_gravity() {
	const double toler = 5.0e-4;
	statistics stats;
	stats = root.get_statistics().get();
	const auto mtot = stats.u[rho_i];
	printf( "Total Mass = %e\n", mtot);
	int oi = 0;
	for (int l = 0; l <= opts.max_level; l++) {
		int pass = 0;
		double r;
		printf("Solving gravity on level %i\n", l);
		//	for (int i = 0; i < 10; i++) {
		do {
			auto tmp = root.gravity_solve(pass, l, std::vector<double>(), 0.0, mtot).get();
			r = tmp.resid /mtot;
			printf("%i %e %e %e\n", pass, r, tmp.resid, tmp.mass);
//			if( r < 5e-03 && l == opts.max_level) {
//				break;
//			}
			pass++;
			if (pass > 250) {
			//	break;
			}
	//		if( l == 6 )
	//		output_silo(std::string("X.") + std::to_string(oi++) + ".silo");

		} while (r > toler);
		//	}
		auto tmp = root.gravity_solve(GRAVITY_FINAL_PASS, l, std::vector<double>(), 0.0, mtot).get();
		r = tmp.resid / mtot;
		printf("%i %e\n", pass, r);
//		output_silo(std::string("X.") + std::to_string(oi++) + ".silo");
//		r = root.gravity_solve(0, l, std::vector<double>(), 0.0).get().resid;
//		printf( "%e\n", r);
//		r = root.gravity_solve(GRAVITY_FINAL_PASS, l, std::vector<double>(), 0.0).get().resid;
//		printf( "%e\n", r);
	}
	auto rms = root.compute_error().get();
	rms = std::sqrt(rms);
	printf( "Error = %e\n",rms);
	return stats;
}

void master(int level, double tmax) {
	if (level > opts.max_level) {
		return;
	}

	double nstep = -1;
	const int refine_freq = opts.window / opts.cfl;
	do {
		printf("Advancing level %i from %e to %e\n", level, tm[level], tm[level] + dt[level]);
		if (level > 0 && nstep == -1) {
			levels_set_child_families(level - 1);
		}
		bool refine = (nstep == -1) && (super_step[level] % refine_freq == 0);
		if (refine) {
			printf("Refining on level %i\n", level);
		}
		dt[level] = levels_hydro_initialize(level, refine);
		if (dt[level] == 0.0) {
			tm[level] = tm[level - 1];
			master(level + 1, tm[level]);
			return;
		}
		levels_show();
		dt[level] = opts.cfl * dx[level] / dt[level];
		nstep = std::ceil((tmax - tm[level]) / dt[level]);
		dt[level] = (tmax - tm[level]) / nstep;
		levels_hydro_substep(level, 0, dt[level]);
//		printf("...\n");
		for (int rk = 1; rk < opts.nrk; rk++) {
			levels_hydro_substep(level, rk, dt[level]);
		}
		tm[level] += dt[level];
		master(level + 1, tm[level]);
	} while (nstep != 1.0);
	super_step[level]++;

}

int hpx_main(int argc, char *argv[]) {
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
	printf("%i\n", sibs.size());
	dx.resize(opts.max_level + 1);
	dt.resize(opts.max_level + 1);
	tm.resize(opts.max_level + 1);
	super_step.resize(opts.max_level + 1);
	for (int l = 0; l <= opts.max_level; l++) {
		super_step[l] = 0;
		dx[l] = 1.0 / (opts.max_box * (1 << l));
		tm[l] = 0.0;
		auto fut = root.initialize(l);
		auto amax = fut.get();
		if (l >= 1) {
			levels_set_child_families(l - 1);
		}
		if (amax != 0.0) {
			dt[l] = opts.cfl * dx[l] / amax;
		} else {
			dt[l] = 0.0;
		}
	}
	root.set_family(tree_client(), root, sibs).get();
	root.restrict_all().get();
//	master(0, 1.0e-100);
	levels_show();
	solve_gravity();
	output_silo("X.0.silo");
//	int i = 0;
//	const auto dt = 0.01;
//	levels_show();
//	for (double t = 0.0; t < opts.tmax; t += dt) {
//		i++;
//		master(0, std::min(t + dt, opts.tmax));
//		std::string fname = "X." + std::to_string(i) + ".silo";
//		output_silo(fname);
//	}
//	std::string fname = "X." + std::to_string(i) + ".silo";
//	output_silo(fname);
	return hpx::finalize();
}

#ifndef HPX_LITE
int main(int argc, char *argv[]) {

	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };

	hpx::init(argc, argv, cfg);
	printf("exiting..\n");
}
#endif
