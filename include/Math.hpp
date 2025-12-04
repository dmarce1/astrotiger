/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "Definitions.hpp"

#include <cassert>
#include <numeric>

template<typename T>
constexpr T iPow(T x, int p) {
	constexpr auto one = T(1);
	if (p < 0) {
		return one / iPow(x, -p);
	}
	T y = x;
	T z = one;
	while (p) {
		if (p & 1) {
			z *= y;
		}
		p >>= 1;
		y *= y;
	};
	return z;
}

constexpr uint64_t factorial(uint64_t n) {
	assert(n <= 20);
	uint64_t fac = 1;
	for (uint64_t i = 2; i <= n; i++) {
		fac *= i;
	}
	return fac;
}

constexpr uint64_t binCo(uint64_t n, uint64_t k) {
	using std::min;
	if (k > n) {
        return 0;
    }
    k = min(k, n - k);
    uint64_t c = 1;
    for (uint64_t j = 1; j <= k; j++) {
        c *= (n + 1 - j);
        c /= j;
    }
    return c;
}
