#include <hpx/hpx_init.hpp>
#include "radiation.hpp"

#define NDIM 3

static constexpr PhysicalConstants<double> pc { };

int hpx_main(int argc, char *argv[]) {
	RadConserved<double, NDIM> vars;
	constexpr auto c = pc.c.value();
	constexpr auto δ = identity<double, NDIM + 1>();
	double E = 1.0;
	double F = c * E / sqrt(1.6);
	for (int d = 0; d < NDIM; d++) {
		vars.setEnergy(E);
		vars.setFlux(0, F);
		vars.setFlux(1, .01);
		vars.setFlux(2, -.2);
		auto const [λ, R] = vars.eigenSystem(d);
		auto Λ = δ;
		for(int k = 0; k < NDIM + 1; k++) {
			Λ(k, k) = λ[k];
		}
		auto const J = vars.jacobian(d);
//		std::cout << "Λ = \n" << Λ;
//		std::cout << "J = \n" << J;
//		std::cout << "R = \n" << R;
//		std::cout << "J * R = \n" << J * R;
//		std::cout << "R * Λ = \n" << R * Λ;
		std::cout << "J * R - R * Λ = \n" << (J * R - R * Λ);
	}

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
