#include <astrotiger/tree.hpp>

int hpx_main(int argc, char *argv[]) {
	options opts;
	opts.process_options(argc, argv);

	multi_range box;
	for( int i = 0; i < NDIM; i++) {
		box.min[i] = 0;
		box.max[i] = opts.max_box;
	}
	auto root = tree::allocate(0,box).get();

	root.initialize().get();

	return hpx::finalize();
}

#ifndef HPX_LITE
int main(int argc, char *argv[]) {

	std::vector < std::string > cfg = { "hpx.commandline.allow_unknown=1" };

	hpx::init(argc, argv, cfg);
	printf("exiting.2..\n");
}
#endif
