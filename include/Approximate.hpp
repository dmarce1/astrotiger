/*
 * Approximate.hpp
 *
 *  Created on: Dec 27, 2025
 *      Author: dmarce1
 */

#pragma once

#include "Math.hpp"

#include <ostream>

template<std::floating_point T, int W = 2>
struct Approximate {
	constexpr Approximate(Approximate const&) = default;
	constexpr Approximate(Approximate&&) = default;
	constexpr Approximate& operator=(Approximate const&) = default;
	constexpr Approximate& operator=(Approximate&&) = default;
	constexpr Approximate() :
			x(0.0), σ(0.0) {
	}
	constexpr Approximate(std::floating_point auto xo = 0.0) :
	x(xo), σ(eps * xo) {
	}

	constexpr Approximate(std::integral auto xo = 0.0) :
	x(xo), σ(eps * xo) {
	}

	constexpr Approximate& operator=(std::floating_point auto xo) {
		x = xo;
		σ = eps * xo;
		return*this;
	}

	explicit constexpr operator T() const {
		if (std::abs(x) < std::abs(σ)) {
			return T(0.0);
		} else {
			return x;
		}
	}
	constexpr Approximate& operator=(std::integral auto xo) {
		x = xo;
		σ = eps * xo;
		return*this;
	}
	constexpr Approximate& operator+=(Approximate const &other) {
		*this = *this + other;
		return *this;
	}
	constexpr Approximate& operator-=(Approximate const &other) {
		*this = *this - other;
		return *this;
	}
	constexpr Approximate& operator*=(Approximate const &other) {
		*this = *this * other;
		return *this;
	}
	constexpr Approximate& operator/=(Approximate const &other) {
		*this = *this / other;
		return *this;
	}
	constexpr Approximate operator+() const {
		return *this;
	}
	constexpr Approximate operator-() const {
		return Approximate(-x, σ);
	}
	constexpr Approximate operator+(Approximate const &B) const {
		Approximate const &A = *this;
		T const x = A.x + B.x;
		T const σ = sqrt(sqr(A.σ) + sqr(B.σ));
		return Approximate(x, σ);
	}
	constexpr Approximate operator-(Approximate const &B) const {
		Approximate const &A = *this;
		return A + (-B);
	}
	constexpr Approximate operator*(Approximate const &B) const {
		Approximate const &A = *this;
		T const x = A.x * B.x;
		T const σ = sqrt(sqr(B.x * A.σ) + sqr(A.x * B.σ));
		return Approximate(x, σ);
	}
	constexpr Approximate operator/(Approximate const &B) const {
		Approximate const &A = *this;
		T const x = A.x / B.x;
		T const σ = sqrt(sqr(A.σ / B.x) + sqr(A.x * B.σ / sqr(B.x)));
		return Approximate(x, σ);
	}
	constexpr bool operator==(Approximate const &C) const {
		auto const &B = *this;
		auto const A = B - C;
		return (std::abs(A.x) <= A.σ);
	}
	constexpr bool operator!=(Approximate const &C) const {
		auto const &B = *this;
		return !(B == C);
	}
	constexpr bool operator<(Approximate const &C) const {
		auto const &B = *this;
		auto const A = B - C;
		return (A.x < -A.σ);
	}
	constexpr bool operator<=(Approximate const &C) const {
		auto const &B = *this;
		auto const A = B - C;
		return (A.x <= A.σ);
	}
	constexpr bool operator>(Approximate const &C) const {
		auto const &B = *this;
		auto const A = B - C;
		return (A.x >= -A.σ);
	}
	constexpr bool operator>=(Approximate const &C) const {
		auto const &B = *this;
		auto const A = B - C;
		return (A.x > A.σ);
	}
	constexpr T absError() const {
		return σ;
	}
	constexpr T relError() const {
		return σ / std::abs(x);
	}
	friend constexpr auto abs(Approximate A) {
		A.x = std::abs(A.x);
		return A;
	}
	friend std::ostream& operator<<(std::ostream &os, Approximate const &A) {
		char *str;
		if (A.x && A.σ) {
			asprintf(&str, "%g±(%.2g)%%", A.x, (100 * A.σ / std::abs(A.x)));
		} else {
			asprintf(&str, "%g", A.x);
		}
		os << std::string(str);
		free(str);
		return os;
	}
	constexpr Approximate(T a, T b) :
			x(a), σ(b) {
	}
	static constexpr T eps = W * std::numeric_limits < T > ::epsilon();
	T x;
	T σ;
};

template<typename T>
struct IsApproximate {
	static constexpr bool value = false;
};

template<typename T, int W>
struct IsApproximate<Approximate<T, W>> {
	static constexpr bool value = true;
};


template<typename T>
concept ApproximateType = IsApproximate<T>::value;

