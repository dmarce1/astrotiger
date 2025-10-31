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

namespace detail {

//inline constexpr size_t basisSize(size_t P, size_t D, BasisT basisT) {
//	if (basisT == BasisT::totalDegree) {
//		return binco(P + D - 1, D);
//	} else {
//		return pow(P, D);
//	}
//}

}

template<typename Type, int order, int dimensionCount>
constexpr auto analyzeLegendre(auto inVector) {
	constexpr auto q = gaussLegendrePoints<Type, order>();
	constexpr int inSize = pow(order, dimensionCount);
	constexpr int outSize = binco(order + dimensionCount - 1, dimensionCount);
	static_assert(inVector.size() == inSize);
	int blockCount = 1;
	int batchSize = inSize / order;
	int blockSize = batchSize * order;
	std::array<Type, inSize> tmpStorage { };
	std::array<int, dimensionCount> indices { };
	for (int dimension = 0; dimension < dimensionCount; dimension++) {
		indices.fill(0);
		Type *sourcePtr = inVector.data();
		for (int block = 0; block < blockCount; block++) {
			int const endOrder = order - std::accumulate(indices.begin(), indices.end(), 0);
			if (endOrder > 0) {
				auto destin = std::span(tmpStorage.data(), blockSize);
				auto source = std::span(sourcePtr, blockSize);
				std::fill_n(destin.begin(), blockSize, Type(0));
				for (int k = 0; k < order; k++) {
					Type Pnp1;
					Type Pn = 1;
					Type Pnm1 = 0;
					for (int n = 0; n < endOrder; n++) {
						for (int index = 0; index < batchSize; index++) {
							destin[batchSize * n + index] += ((2 * Type(n) + 1) / 2) * Pn * source[batchSize * k + index] * q.w[k];
						}
						if (n + 1 < order) {
							Pnp1 = ((2 * Type(n) + 1) * q.x[k] * Pn - Type(n) * Pnm1) / Type(n + 1);
							Pnm1 = Pn;
							Pn = Pnp1;
						}
					}
				}
				std::copy_n(destin.begin(), endOrder * batchSize, source.begin());
			}
			sourcePtr += blockSize;
			if (block + 1 < blockCount) {
				int dimIndex = 0;
				while (++indices[dimIndex] == order) {
					indices[dimIndex++] = 0;
				}
			}
		}
		blockCount *= order;
		blockSize = batchSize;
		batchSize /= order;
	}
	std::array<Type, outSize> outVector { };
	indices.fill(0);
	for (int inIndex = 0; inIndex < inSize; inIndex++) {
		int outIndex = 0, degree = 0;
		for (int dimension = 0; (dimension < dimensionCount) && (degree < order); dimension++) {
			degree += indices[dimension];
			outIndex += binco(degree + dimension, dimension + 1);
		}
		if (degree < order) {
			outVector[outIndex] = inVector[inIndex];
		}
		if (inIndex + 1 < inSize) {
			int dimIndex = dimensionCount - 1;
			while (++indices[dimIndex] == order) {
				indices[dimIndex--] = 0;
			}
		}
	}
	return outVector;
}
//
//for (size_t k = 0; k < P; k++) {
//	T Pnm1 = T(0);
//	T Pn = T(1);
//	T Pnp1;
//	size_t const Mlo = order - k;
//	SHOW(Mlo);
//		for (size_t n = 0; n < order; n++) {
////						F[dsti + ilo] += (n + T(0.5)) * Pn * f[srci + ilo] * q.w[k];
//		}
////					Pnp1 = (T(2 * n + 1) * q.x[k] * Pn - T(n) * Pnm1) / T(n + 1);
////					Pnm1 = Pn;
////					Pn = Pnp1;
////					dsti += Mlo;
//	srci += sourceStride;
//}
