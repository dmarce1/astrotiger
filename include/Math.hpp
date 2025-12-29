/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include <cmath>
#include <concepts>
#include <cstdint>

#define MONO_INCREASES 0x1
#define MONO_DECREASES 0x2
#define MONO_EQUALS 0x4

inline constexpr auto sqr(auto x) {
	return x * x;
}

constexpr auto binco(std::integral auto n, std::integral auto k) {
	using std::min;
	using type = decltype(n);
	if (k > n) {
		return 0;
	}
	k = min(k, n - k);
	type c = 1;
	for (type j = 1; j <= k; j++) {
		c *= (n + 1 - j);
		c /= j;
	}
	return c;
}

template <std::integral auto N>
consteval int factorial() {
	if constexpr (N == 0) {
		return 1;
	} else {
		return N * factorial<N - 1>();
	}
}

constexpr auto factorial(std::integral auto n) {
	using type = decltype(n);
	type fac = 1;
	for (type i = 2; i <= n; i++) {
		fac *= i;
	}
	return fac;
}

constexpr auto icopysign(auto x, auto y) {
	x = (x > 0) ? x : -x;
	if (y >= 0) {
		return x;
	} else {
		return -x;
	}
}

constexpr auto ipow(auto x, std::integral auto p) -> decltype(x) {
	using type = decltype(x);
	constexpr type one = type(1);
	if (p < 0) {
		return one / ipow(x, -p);
	}
	type y = x;
	type z = one;
	while (p) {
		if (p & 1) {
			z *= y;
		}
		p >>= 1;
		y *= y;
	};
	return z;
}

constexpr int nonepow(std::integral auto p) {
	if ((p & 1) == 0) {
		return +1;
	} else {
		return -1;
	}
}

template <int P>
consteval int nonepow() {
	if constexpr ((P & 1) == 0) {
		return +1;
	} else {
		return -1;
	}
}


template <int λo_, int... λ>
consteval int monotonicity() {
	int mt = 0;
	int λo = λo_;
	(((mt |= ((λ > λo) ? MONO_INCREASES : ((λ < λo) ? MONO_DECREASES : MONO_EQUALS))), λo = λ), ...);
	return mt;
}

template <int... λ>
consteval bool strictlyIncreasing() {
	constexpr int mt = monotonicity<λ...>();
	return mt == MONO_INCREASES;
}

template <int... λ>
consteval bool nonDecreasing() {
	constexpr int mt = monotonicity<λ...>();
	return ((mt & MONO_DECREASES) == 0);
}

template <int... λ>
consteval bool nonIncreasing() {
	constexpr int mt = monotonicity<λ...>();
	return ((mt & MONO_INCREASES) == 0);
}

template <int... λ>
consteval bool strictlyDecreasing() {
	constexpr int mt = monotonicity<λ...>();
	return mt == MONO_DECREASES;
}

