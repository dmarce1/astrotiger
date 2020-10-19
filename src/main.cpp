#include <astrotiger/tree.hpp>
#include <astrotiger/levels.hpp>

static tree_client root;

double master(int level) {
	if (level > 0) {
		levels_set_child_families(level - 1);
	}
	if (level < opts.max_level) {
		master(level + 1);
	}
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
	for (int l = 0; l < opts.max_level; l++) {
		root.initialize(l).get();
	}
	root.set_family(tree_client(), root, sibs).get();
	master(0);
	return hpx::finalize();
}

#ifndef HPX_LITE
int main(int argc, char *argv[]) {

	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };

	hpx::init(argc, argv, cfg);
	printf("exiting..\n");
}
#endif
