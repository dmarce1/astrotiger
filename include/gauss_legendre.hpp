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

template<typename T, int P, int N>
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
		for (int n = 0; n < P; n++) {
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
	if constexpr (P == ((N << 1) | 1)) {
		x = zero;
	} else {
		T θ, θ1;
		θ1 = π * (one - half * (T(2 * N + 1)) / T(P));
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

template<typename T, int P, int N = 0>
constexpr auto gaussLegendrePoints() {
	QuadraturePoint<std::array<T, P>> rc;
	if constexpr (N < P) {
		rc = gaussLegendrePoints<T, P, N + 1>();
		auto const pt = gaussLegendrePoint<T, P, N>();
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

template<typename T, int P, int D>
constexpr auto analyzeLegendre(auto f) {
	constexpr auto q = gaussLegendrePoints<T, P>();
	constexpr int inSize = pow(P, D);
	constexpr int outSize = binco(P + D - 1, D);
	static_assert(f.size() == inSize);
	int Nhi = 1;
	int PNlo = inSize;
	int Nlo = inSize / P;
	std::array<T, inSize> g { };
	std::array<int, D> idx { };
	for (int d = 0; d < D; d++) {
		idx.fill(0);
		T *srcPointer = f.data();
		for (int nhi = 0; nhi < Nhi; nhi++) {
			int const Pend = P - std::accumulate(idx.begin(), idx.end(), 0);
			if (Pend > 0) {
				auto dst = std::span(g.data(), PNlo);
				auto src = std::span(srcPointer, PNlo);
				std::fill_n(dst.begin(), PNlo, T(0));
				for (int k = 0; k < P; k++) {
					T Pnp1;
					T Pn = T(1);
					T Pnm1 = T(0);
					for (int n = 0; n < Pend; n++) {
						for (int nlo = 0; nlo < Nlo; nlo++) {
							dst[Nlo * n + nlo] += (T(2 * n + 1) / T(2)) * Pn * src[Nlo * k + nlo] * q.w[k];
						}
						if (n + 1 < P) {
							Pnp1 = (T(2 * n + 1) * q.x[k] * Pn - T(n) * Pnm1) / T(n + 1);
							Pnm1 = Pn;
							Pn = Pnp1;
						}
					}
				}
				std::copy_n(dst.begin(), Pend * Nlo, src.begin());
			}
			srcPointer += PNlo;
			if (nhi + 1 < Nhi) {
				int q = 0;
				while (++idx[q] == P) {
					idx[q++] = 0;
				}
			}
		}
		Nhi *= P;
		PNlo /= P;
		Nlo /= P;
	}
	std::array<T, outSize> F { };
	idx.fill(0);
	for (int j = 0; j < inSize; j++) {
		int k = 0, deg = 0;
		for (int dim = 0; (dim < D) && (deg < P); dim++) {
			deg += idx[dim];
			k += binco(deg + dim, dim + 1);
		}
		if (deg < P) {
			F[k] = f[j];
		}
		if (j + 1 < inSize) {
			int d = D - 1;
			while (++idx[d] == P) {
				idx[d--] = 0;
			}
		}
	}
	return F;
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
