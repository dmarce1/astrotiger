#pragma once

#include <algorithm>
#include <limits>
#include <numbers>

#include "math.hpp"

template<typename T>
struct QuadraturePoint {
	T x;
	T w;
};

template<typename T, int order, int N>
constexpr auto gaussLegendrePoint() {
	using std::abs;
	using std::atan;
	using std::cos;
	using std::sin;
	using std::numeric_limits;
	static constexpr T ε = 4 * numeric_limits < T > ::epsilon();
	static constexpr T zero = T(0);
	static constexpr T one = T(1);
	static constexpr T two = T(2);
	static constexpr T half = one / two;
	static constexpr T π = std::numbers::pi_v<T>;
	constexpr auto legendre = [](T x) {
		T Pn, Pnp1, Pnm1;
		T dPndx, dPnp1dx, dPnm1dx;
		Pn = Pnm1 = one;
		dPndx = dPnm1dx = zero;
		for (int n = 0; n < order; n++) {
			Pnp1 = (T(2 * n + 1) * x * Pn - T(n) * Pnm1) / T(n + 1);
			dPnp1dx = (T(2 * n + 1) * (Pn + x * dPndx) - T(n) * dPnm1dx) / T(n + 1);
			Pnm1 = Pn;
			dPnm1dx = dPndx;
			Pn = Pnp1;
			dPndx = dPnp1dx;
		}
		return std::pair(Pn, dPndx);
	};
	T x, w;
	if constexpr (order == ((N << 1) | 1)) {
		x = zero;
	} else {
		T θ, θ1;
		θ1 = π * (one - half * (T(2 * N + 1)) / T(order));
		do {
			θ = θ1;
			x = cos(θ);
			auto const [Pn, dPndx] = legendre(x);
			θ1 = θ + Pn / (sin(θ) * dPndx);
		} while (abs(θ - θ1) > ε);
	}
	auto const [Pn, dPndx] = legendre(x);
	w = two / (sqr(dPndx) * (one - sqr(x)));
	return QuadraturePoint<T> { x, w };
}

template<typename T, int order, int N = 0>
constexpr auto gaussLegendrePoints() {
	QuadraturePoint<std::array<T, order>> rc;
	if constexpr (N < order) {
		rc = gaussLegendrePoints<T, order, N + 1>();
		auto const pt = gaussLegendrePoint<T, order, N>();
		rc.x[N] = pt.x;
		rc.w[N] = pt.w;
	}
	return rc;
}

#define SHOW1( a ) printf("%s = %i\n", #a, a)
#define SHOW2( a, b) printf("%s = %i %s = %i\n", #a, a, #b, b)
#define SHOW3( a, b, c ) printf("%s = %i %s = %i %s = %i\n", #a, a, #b, b, #c, c)
#define SHOW5( a, b, c, d , e ) printf("%s = %i %s = %i %s = %i %s = %i %s = %i\n", #a, a, #b, b, #c, c, #d, d, #e, e)

//#define SHOW1( a )
//#define SHOW2( a, b, c )
//#define SHOW3( a, b, c )
//#define SHOW4( a, b, c, d )

template<typename Type, int order, int dimensionCount>
struct LegendreBasis {
	static constexpr int physicalSize = pow(order, dimensionCount);
	static constexpr int spectralSize = binco(order + dimensionCount - 1, dimensionCount);
	static constexpr auto qPts = gaussLegendrePoints<Type, order>();
private:
	static constexpr auto zero = Type(0);
	static constexpr auto one = Type(1);
	static constexpr auto two = Type(2);
	static constexpr auto half = one / two;
	static constexpr int totalDegree(int index) {
		if (index == 0) {
			return 0;
		}
		return totalDegree(index / order) + (index % order);
	}
	static constexpr int nextIndex(int idx) {
		int deg = totalDegree(idx);
		int beg = idx + 1;
		while (deg < order) {
			for (idx = beg; idx != physicalSize; idx++) {
				if (totalDegree(idx) == deg) {
					return idx;
				}
			}
			deg++;
			beg = 0;
		}
		return physicalSize;
	}
	template<int dimension, int totalDegree, bool doSynthesis, int ... blockIndexes>
	static void compute(Type *source, std::integer_sequence<int, blockIndexes...>) {
		constexpr auto batchSize = pow<dimensionCount - dimension - 1>(order);
		constexpr auto blockSize = batchSize * order;
		(compute<dimension, blockIndexes, totalDegree + blockIndexes, doSynthesis>(source + blockIndexes * blockSize), ...);
	}
	template<int dimension, int totalDegree, bool doSynthesis>
	static void compute(Type *source) {
		compute<dimension, totalDegree, doSynthesis>(source, std::make_integer_sequence<int, order> { });
	}
	template<int dimension, int block, int totalDegree, bool doSynthesis>
	static constexpr auto compute(Type *source) {
		if constexpr (dimension == dimensionCount) {
			return;
		} else {
			if constexpr (doSynthesis) {
				compute<dimension + 1, totalDegree, true>(source);
			}
			constexpr auto batchSize = pow<dimensionCount - dimension - 1>(order);
			auto const endOrder = order - totalDegree;
			auto const copySize = (doSynthesis ? order : endOrder) * batchSize;
			std::array<Type, order * batchSize> destin { };
			if (endOrder > 0) {
				for (int k = 0; k < order; k++) {
					Type Pn = one, Pnm1 = zero, Pnp1;
					for (int n = 0; n < endOrder; n++) {
						for (int index = 0; index < batchSize; index++) {
							if constexpr (doSynthesis) {
								destin[batchSize * k + index] += Pn * source[batchSize * n + index];
							} else {
								destin[batchSize * n + index] += Pn * source[batchSize * k + index] * (half * Type(2 * n + 1)) * qPts.w[k];
							}
						}
						if (n + 1 < order) {
							Pnp1 = (Type(2 * n + 1) * qPts.x[k] * Pn - Type(n) * Pnm1) / Type(n + 1);
							Pnm1 = Pn;
							Pn = Pnp1;
						}
					}
				}
				std::copy_n(destin.begin(), copySize, source);
			}
			if constexpr (!doSynthesis) {
				compute<dimension + 1, totalDegree, false>(source);
			}
		}
	}
public:
	static constexpr auto analyze(std::array<Type, physicalSize> inVector) {
		compute<0, 0, 0, false>(inVector.data());
		std::array<Type, spectralSize> outVector { };
		int inIndex = 0;
		for (int outIndex = 0; outIndex < spectralSize; outIndex++) {
			outVector[outIndex] = inVector[inIndex];
			inIndex = nextIndex(inIndex);
		}
		return outVector;
	}
	static constexpr auto synthesize(std::array<Type, spectralSize> inVector) {
		std::array<Type, physicalSize> outVector { };
		int outIndex = 0;
		for (int inIndex = 0; inIndex < spectralSize; inIndex++) {
			outVector[outIndex] = inVector[inIndex];
			outIndex = nextIndex(outIndex);
		}
		compute<0, 0, 0, true>(outVector.data());
		return outVector;
	}
};

