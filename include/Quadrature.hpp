#pragma once

#include "Definitions.hpp"

#include <vector>

struct QuadraturePoint {
	Real x;
	Real w;
};

std::vector<QuadraturePoint> getLobattoQuadraturePoints(int);

template<typename T>
constexpr T legendre(int N, T x) {
	T Pnp1, Pnm1;
	T Pn = T(1);
	if (N >= 1) {
		Pnp1 = x;
		Pnm1 = Pn;
		Pn = Pnp1;
		for (int n = 1; n < N; n++) {
			Pnp1 = (T(2 * n + 1) * x * Pn - T(n) * Pnm1) / T(n + 1);
			Pnm1 = Pn;
			Pn = Pnp1;
		}
	}
	return Pn;
}

template<typename T, int N>
constexpr auto getLegendreQuadraturePoints() {
	std::array<T, N> x;
	std::array<T, N> w;
	T const Np1 = T(N + 1);
	T const N_ = T(N);
	for (int j = 0; j < N; j++) {
		x[j] = std::cos(std::numbers::pi_v < T > *(T(1) - (T(j) + (T(1) / T(2))) / N_));
		T xPrevious, Pnp1;
		do {
			xPrevious = x[j];
			Pnp1 = legendre(N + 1, x[j]);
			T const Pn = legendre(N, x[j]);
			T const x2 = squared(x[j]);
			T const num = (T(1) - x2) * Pn;
			T const den = Np1 * (Pnp1 - x[j] * Pn);
			x[j] += num / den;
		} while (T(x[j]) != T(xPrevious));
		T const dpdx = Np1 * (legendre(N + 1, x[j]) - x[j] * legendre(N, x[j])) / (squared(x[j]) - T(1));
		w[j] = T(2) / ((T(1) - squared(x[j])) * squared(dpdx));
	}
	return std::pair<std::array<T, N>, std::array<T, N>>({x, w});
}
