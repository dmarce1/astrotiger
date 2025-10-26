#include <hpx/hpx_init.hpp>
#include "gas_flux.hpp"
#include "rad_flux.hpp"

#define NDIM 3

static constexpr PhysicalConstants<double> pc { };

int hpx_main(int argc, char *argv[]) {
	constexpr EquationOfState<double> eos(1.0);
	GasConserved<double, NDIM> Ul;
	GasConserved<double, NDIM> Ur;
	riemannHLLC(Ul, Ur, eos, 1);

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
