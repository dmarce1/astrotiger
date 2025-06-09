#include <hpx/hpx_init.hpp>

#include "HyperGrid.hpp"
#include "Options.hpp"
#include "Util.hpp"
#include "Triangular.hpp"

void testq();

template<int dimensionCount, int modeCount>
auto triangularMap() {
	constexpr int size = binomialCoefficient(modeCount + dimensionCount - 1, dimensionCount);
	std::array<std::array<int, dimensionCount>, size> map;
//	for (int flatIndex = 0; flatIndex < size; flatIndex++) {
//		for (int dimension = dimensionCount - 1; dimension >= 0; dimension--) {
//			auto thisIndex = triangleIndex[dimension];
//			thisIndex -= triangleIndex[dimension + 1];
//			flatIndex += binomialCoefficient(thisIndex + dimensionCount - 1, dimensionCount);
//		}
//	}
	return map;
}

template<typename Type, int modeCount, int dimensionCount>
std::array<Type, binomialCoefficient(modeCount + dimensionCount - 1, dimensionCount)> nodalToModal(std::array<Type, power(modeCount, dimensionCount)> x) {
	static auto const quadrature = getLegendreQuadraturePoints(modeCount);
	int blockCount = power(modeCount, dimensionCount - 1);
	std::array<Type, power(modeCount, dimensionCount)> y;
	std::array<Type, binomialCoefficient(modeCount + dimensionCount - 1, dimensionCount)> z;
	for (int dimension = dimensionCount - 1; dimension >= 0; dimension--) {
		int outputIndex = 0;
		int lowSize = binomialCoefficient(modeCount + dimensionCount - dimension - 2, dimensionCount - dimension - 1);
		for (int hiIndex = 0; hiIndex < blockCount; hiIndex++) {
			for (int modeIndex = 0; modeIndex < modeCount; modeIndex++) {
				int lowMaximum = binomialCoefficient(modeCount - modeIndex + dimensionCount - dimension - 2, dimensionCount - dimension - 1);
				for (int lowIndex = 0; lowIndex < lowMaximum; lowIndex++) {
					y[outputIndex] = Type(0);
					for (int nodeIndex = 0; nodeIndex < modeCount; nodeIndex++) {
						y[outputIndex] +=
								quadrature[nodeIndex].w * std::legendre(modeIndex, quadrature[nodeIndex].x) * (Type(modeIndex) + Type(0.5)) * x[hiIndex * modeCount * lowSize + nodeIndex * lowSize + lowIndex];
					}
					outputIndex++;
				}
			}
		}
		blockCount /= modeCount;
		std::swap(x, y);
	}
	std::copy(x.begin(), x.begin() + z.size(), z.begin());
	return z;
}

int hpx_main(int argc, char *argv[]) {
	enableFPE();
	processOptions(argc, argv);
	using namespace Constants;
	constexpr int D = 3;
	constexpr int M = 3;
	std::array<double, power(M, D)> x { };
	std::fill(x.begin(), x.end(), 0);
	static auto const Q = getLegendreQuadraturePoints(M);
	for (int i = 0; i < power(M, D); i++) {
		x[i] = std::legendre(1, Q[(i / M) % M].x);
	}
	auto y = nodalToModal<double, M, D>(x);
	for (int i = 0; i < binomialCoefficient(M + D - 1, D); i++) {
		printf("n=%i  y=%e\n", i, y[i]);
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
