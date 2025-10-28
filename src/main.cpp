#include <hpx/hpx_init.hpp>
#include "gas_flux.hpp"
#include "rad_flux.hpp"
#include "interval.hpp"
#include "multi_array.hpp"

#define NDIM 3

static constexpr PhysicalConstants<double> pc { };

int hpx_main(int argc, char *argv[]) {
	constexpr std::array<int, 3> dims { { 4, 4, 8 } };
	MultiArray<double, 3, dims> arr{};
	constexpr Interval<int, 3> I(1, 8);
	constexpr auto rng = I.literal();
	arr.subarray<rng>();
	auto test = arr.transpose<0, 2>();

	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
	enableFPE();
	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };
	cfg.push_back("hpx.stacks.small_size=1048576");
	hpx::init_params init_params;
	init_params.cfg = std::move(cfg);
	auto rc = hpx::init(argc, argv, init_params);
	return rc;
}
