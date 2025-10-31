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
	std::array<double, pow(order, dimensionCount)> f;
	constexpr auto q = gaussLegendrePoints<double, order>();
	for (int l = 0; l < order; l++) {
		for (int m = 0; m <= l; m++) {
			for (int n = 0; n <= m; n++) {
				auto func = [l, m, n](double x, double y, double z) {
					return std::legendre(l - m, z) * std::legendre(m - n, y) * std::legendre(n, x);
				};
				int p = 0;
				for (int i = 0; i < order; i++) {
					for (int j = 0; j < order; j++) {
						for (int k = 0; k < order; k++) {
							f[p++] = func(q.x[i], q.x[j], q.x[k]);
						}
					}
				}
				printf("%i %i %i: ", n, m - n, l - m);
				auto const F = analyzeLegendre<double, order, dimensionCount>(std::move(f));
				for (size_t p = 0; p < F.size(); p++) {
					printf("%i \t", (int) std::round(F[p]));
				}
				printf("\n");
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
