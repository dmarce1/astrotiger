/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#include "Definitions.hpp"

#include "Bits.hpp"
#include "FixedPrecision.hpp"
#include "KDTree.hpp"

#include <hpx/hpx_init.hpp>

int hpx_main(int argc, char *argv[]) {
	auto test = KDTree<3>::create();
	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
	hpx::init_params init_params;
	std::string sizeString = "hpx.stacks.small_size=" + std::to_string(size_t(1) << size_t(18)) + "\n";
	init_params.cfg.push_back("hpx.commandline.allow_unknown=1");
	init_params.cfg.push_back(sizeString);
	return hpx::init(argc, argv, init_params);
}
