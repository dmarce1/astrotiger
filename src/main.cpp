#include <astrotiger/tree.hpp>
#include <astrotiger/levels.hpp>
#include <astrotiger/output.hpp>

static tree_client root;

std::vector<double> dx;
std::vector<double> dt;
std::vector<double> tm;
std::vector<int> super_step;

void solve_gravity(int l, double dt) {
	const double toler = 5.0e-4;
	int pass = 0;
	double r;
	printf("Solving gravity on level %i\n", l);
	auto tmp2 = root.set_gravity_source(l, dt).get();
	auto src = tmp2.first;
	static double this_mass;
	if (tmp2.second) {
		multi_range box;
		for (int i = 0; i < NDIM; i++) {
			box.min[i] = 0;
			box.max[i] = opts.max_box;
		}
		this_mass = 0.0;
		double dx = 1.0 / opts.max_box;
		for (multi_iterator i(box); !i.end(); i++) {
			this_mass += src[i] * std::pow(dx, NDIM);
	//		printf( "%e\n", src[i]);
		}
	//	printf( "!!! %e\n", this_mass);
	}
	do {
		auto tmp = root.gravity_solve(pass, l, std::vector<double>(), this_mass).get();
		r = tmp.resid;
//		printf("%i %e\n", pass, r);
		if (pass > 200) {
			printf("Gravity solver failed to converged\n");
			abort();
		}
		pass++;
	} while (r > toler);
	auto tmp = root.gravity_solve(GRAVITY_FINAL_PASS, l, std::vector<double>(), this_mass).get();
//	std::string fname = "X." + std::to_string(oi++) + ".silo";
//	output_silo(fname);
	r = tmp.resid / this_mass;
//	printf("%i %e\n", pass, tmp.resid);
}

void master(int level, double tmax) {
	if (level > opts.max_level) {
		return;
	}
	int oi = 0;
	double nstep = -1;
	const int refine_freq = opts.window / opts.cfl;
	statistics stats;
	if (opts.self_gravity) {
		stats = root.get_statistics(opts.max_level).get();
	}
	do {
		if (level > 0 && nstep == -1) {
			levels_set_child_families(level - 1);
		}
		bool refine = (nstep == -1) && (super_step[level] % refine_freq == 0);
		if (refine) {
			printf("Refining on level %i\n", level);
		}
//		std::string fname = "X." + std::to_string(oi++) + ".silo";
//		output_silo(fname);

		printf("Hydro pre-step\n");
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
		printf("Advancing level %i from %e to %e\n", level, tm[level], tm[level] + dt[level]);
		if (opts.self_gravity) {
			//		assert(tmp.u.size());
			solve_gravity(level, 0.5 * dt[level]);
		}
//		printf("Rk1\n");
		levels_hydro_substep(level, 0, dt[level]);
//		printf("Rk2\n");
		levels_hydro_substep(level, 1, dt[level]);
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
	dx.resize(opts.max_level + 1);
	dt.resize(opts.max_level + 1);
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
		if (amax != 0.0) {
			dt[l] = opts.cfl * dx[l] / amax;
		} else {
			dt[l] = 0.0;
		}
		printf("Done\n");
	}
	root.restrict_all().get();
//	master(0, 1.0e-100);
	levels_show();
	if (opts.self_gravity) {
		const auto tmp = root.get_statistics(opts.max_level).get();
		for (int l = 0; l <= opts.max_level; l++) {
			//	assert(tmp.u.size());
			solve_gravity(l, 0.0);
		}
	}
	output_silo("X.0.silo");
	int i = 0;
	const auto dt = 0.01;
	levels_show();
	for (double t = 0.0; t < opts.tmax; t += dt) {
		i++;
		master(0, std::min(t + dt, opts.tmax));
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
