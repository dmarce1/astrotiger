#include "Analysis.hpp"

#include <hpx/hpx_init.hpp>
#include "Options.hpp"
#include "MultiIndex.hpp"
#include "HyperGrid.hpp"
#include "EulerState.hpp"
#include "Real.hpp"
#include "RungeKutta.hpp"

#include "Quadrature.hpp"

#include <numbers>

int hpx_main(int argc, char *argv[]) {
	enableFPE();
	printf("Reading options...\n");
	processOptions(argc, argv);

//	constexpr int D = 3;
//	constexpr int M = 4;
//	constexpr int N = 4;
//
//	std::array<double, binco(M + D - 1, D)> x { };
//	std::fill(x.begin(), x.end(), 0);
//	for (int i = 0; i < binco(N + D - 1, D); i++) {
//		//	x[i] = std::legendre(2, quadX[(i)%N]) + std::legendre(1, quadX[(i/N/N)%N]);
//		x[i] = 2.0 * (rand() / double(RAND_MAX)) - 1.0;
//	}
//	auto y = modalSynthesis<double, D, M, N>(x);
//	auto z = modalAnalysis<double, D, M, N>(y);
//	double err = 0.0;
//	for (int i = 0; i < binco(M + D - 1, D); i++) {
//		err += sqr(z[i] - x[i]);
//		printf("n=%i of %i:  %e %e %e\n", i, power(N, D), z[i], x[i], z[i] - x[i]);
//	}
//	printf("err = %e\n", err);

	using T = Real;
	constexpr int P = 3;
	constexpr int D = 2;
	constexpr int N = 128;
	using RK = typename RungeKutta<T, P>::type;
	using S = EulerState<T, D>;
	HyperGrid<S, N, P, RK> grid;
	grid.initialize(initSodShockTube<T, D>);
	grid.enforceBoundaryConditions();
	T t = T(0);
	T tmax = T(.15);
	T dt;
	RK const rk;
	int iter = 0;
	while (t < tmax) {
		grid.output("X", iter, t);
		std::cout << "i = " << std::to_string(iter);
		std::cout << "  t = " << std::to_string(t);
		dt = grid.beginStep();
		std::cout << "  dt = " << dt << std::endl;
		for (int s = 0; s < rk.stageCount(); s++) {
			grid.subStep(dt, s);
			grid.enforceBoundaryConditions();
		}
		grid.endStep();
		grid.enforceBoundaryConditions();
		iter++;
		t += dt;
	}
	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
#ifndef NDEBUG
	installFpeHandler();
#endif
	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };
	cfg.push_back("hpx.stacks.small_size=524288");
	hpx::init_params init_params;
	init_params.cfg = std::move(cfg);
	auto rc = hpx::init(argc, argv, init_params);
	return rc;
}
