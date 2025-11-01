#include <hpx/hpx_init.hpp>
#include "gauss_legendre.hpp"
#include "dgTransforms.hpp"
#include "gas_flux.hpp"
#include "rad_flux.hpp"
#include "interval.hpp"
#include "multi_array.hpp"

#define NDIM 3

static constexpr PhysicalConstants<double> pc { };

int hpx_main(int argc, char *argv[]) {
	enableFPE();
	constexpr int order = 3;
	constexpr int dimensionCount = 3;
	LegendreBasis<double, order, dimensionCount> basis;
	std::array<double, pow(order, dimensionCount)> f;
	constexpr auto q = gaussLegendrePoints<double, order>();
	for (int degree = 0; degree < order; degree++) {
		for (int l = 0; l <= degree; l++) {
			for (int m = 0; m + l <= degree; m++) {
				int const n = degree - (m + l);
				if (n >= 0) {
					auto func = [l, m, n](double x, double y, double z) {
						return std::legendre(l, x) * std::legendre(m, y) * std::legendre(n, z);
					};
					int b = 0;
					for (int i = 0; i < order; i++) {
						for (int j = 0; j < order; j++) {
							for (int a = 0; a < order; a++) {
								f[b++] = func(q.x[i], q.x[j], q.x[a]);
							}
						}
					}
					printf("%i %i %i: ", l, m, n);
					auto F = basis.analyze(f);
					for (size_t p = 0; p < F.size(); p++) {
						printf("%i ", (int) std::round(F[p]));
					}
					printf("\n");
					printf("%i %i %i: ", l, m, n);
					auto const g = basis.synthesize(F);
					for (size_t p = 0; p < g.size(); p++) {
						printf("%e %e\n ", f[p], g[p]);
					}
					printf("\n");
				}
			}
		}
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
