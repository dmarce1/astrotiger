#include <astrotiger/tree.hpp>
#include <astrotiger/levels.hpp>
#include <astrotiger/output.hpp>

static tree_client root;

std::vector<double> dx;
std::vector<double> dt;
std::vector<double> tm;

void master(int level, double tmax) {
	if (level > opts.max_level) {
		return;
	}

	double nstep = -1;
	do {
		if (level > 0 && nstep == -1) {
			levels_set_child_families(level - 1);
		}
		bool refine = nstep == -1;
		dt[level] = levels_hydro_initialize(level, refine);
		dt[level] = opts.cfl * dx[level] / dt[level];
		nstep = std::ceil((tmax - tm[level]) / dt[level]);
		dt[level] = (tmax - tm[level]) / nstep;
		printf("Advancing level %i from %e to %e\n", level, tm[level], tm[level] + dt[level]);
		levels_hydro_substep(level, 0, dt[level]);
		printf("...\n");
		levels_hydro_substep(level, 1, dt[level]);
		tm[level] += dt[level];
		master(level + 1, tm[level]);
	} while (nstep != 1.0);

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
			const auto shift = i.index() * opts.max_box;
			sibling sib;
			sib.client = root;
			sib.shift = shift;
			sibs.push_back(sib);
		}
	}
	dx.resize(opts.max_level + 1);
	dt.resize(opts.max_level + 1);
	tm.resize(opts.max_level + 1);
	for (int l = 0; l <= opts.max_level; l++) {
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
	output_silo("X.0.silo");
	master(0, 0.25);
	output_silo("X.1.silo");
	return hpx::finalize();
}

#ifndef HPX_LITE
int main(int argc, char *argv[]) {

	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };

	hpx::init(argc, argv, cfg);
	printf("exiting..\n");
}
#endif
