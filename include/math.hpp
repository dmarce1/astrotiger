#pragma once

#include <type_traits>

#include "rational.hpp"

inline constexpr intmax_t binomialCoefficient(intmax_t n, intmax_t k) {
	constexpr intmax_t one(1);
	intmax_t num = one;
	intmax_t den = one;
	for (intmax_t i = one; i <= k; i++) {
		num *= n + one - i;
		den *= i;
	}
	return num / den;
}

template<typename BaseType, typename ExponentType, std::enable_if_t<std::is_integral_v<ExponentType>, int> = 0>
inline constexpr BaseType pow(BaseType x, ExponentType n) {
	if (n >= ExponentType(0)) {
		BaseType xm = x;
		BaseType xn = BaseType(1);
		while (n) {
			if (n & ExponentType(1)) {
				xn *= xm;
			}
			n >>= ExponentType(1);
			if (n) {
				xm *= xm;
			}
		}
		return xn;
	} else {
		return BaseType(1) / pow(x, -n);
	}
}

template<auto exponent, typename BaseType = double, std::enable_if_t<std::is_integral_v<decltype(exponent)>, int> = 0>
inline constexpr BaseType pow(BaseType x) {
	return pow(x, exponent);
}

inline constexpr int nonepow(int n) {
	using std::abs;
	return (abs(n) & 1) ? -1 : +1;
}

inline constexpr intmax_t factorial(intmax_t n) {
	constexpr intmax_t zero(0);
	constexpr intmax_t one(1);
	if (n > zero) {
		return n * factorial(n - one);
	} else {
		return one;
	}
}

inline constexpr Rational factorialPower(Rational x, intmax_t n) {
	constexpr intmax_t one(1);
	if (n) {
		n--;
		return (x - Rational(n)) * factorialPower(x, n);
	} else {
		return one;
	}
}

template<typename T>
inline constexpr auto sqr(T r) {
	return r * r;
}

template<auto n, size_t r>
struct CeilDiv {
	static constexpr auto value = (((n + r - 1) / r));
};

template<auto n, size_t r>
struct CeilRound {
	static constexpr auto value = r * ceil_div(n, r);
};

