/*
 * Real.hpp
 *
 *  Created on: Dec 11, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_REAL_HPP_
#define INCLUDE_REAL_HPP_

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stacktrace>
#include <stdexcept>
#include <type_traits>

namespace Math {

#ifndef NDEBUG
#define CHECK true
#else
#define CHECK  false
#endif

struct Real {
	using Type = double;
	constexpr Real() {
		if constexpr (CHECK) {
			value = std::numeric_limits<Type>::signaling_NaN();
		}
	}
	constexpr explicit Real(double a) {
		value = Type(a);
	}
	Real& operator=(Real const &a) {
		value = a.value;
		return *this;
	}
	constexpr operator Type() const {
		return value;
	}
	Real operator+() const {
		debug_check(*this);
		return *this;
	}
	Real operator-() const {
		debug_check(*this);
		return Real(-value);
	}
	Real operator+(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return Real(value + a.value);
	}
	Real operator-(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return Real(value - a.value);
	}
	Real operator*(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return Real(value * a.value);
	}
	Real operator/(Real const &a) const {
		Real result;
		zero_check(a);
		debug_check(a);
		debug_check(*this);
		result.value = value / a.value;
		return result;
	}
	Real& operator+=(Real const &a) {
		debug_check(a);
		debug_check(*this);
		*this = *this + a;
		return *this;
	}
	Real& operator-=(Real const &a) {
		debug_check(a);
		debug_check(*this);
		*this = *this - a;
		return *this;
	}
	Real& operator*=(Real const &a) {
		debug_check(a);
		debug_check(*this);
		*this = *this * a;
		return *this;
	}
	Real& operator/=(Real const &a) {
		zero_check(a);
		debug_check(a);
		*this = *this / a;
		return *this;
	}
	bool operator==(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return value == a.value;
	}
	bool operator!=(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return value != a.value;
	}
	bool operator<=(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return value <= a.value;
	}
	bool operator>=(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return value >= a.value;
	}
	bool operator<(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return value < a.value;
	}
	bool operator>(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return value > a.value;
	}
	static Real zero() {
		Real z;
		z.value = Type(0);
		return z;
	}
	static Real tiny() {
		Real z;
		z.value = Type(std::numeric_limits<double>::min());
		return z;
	}
	friend Real abs(Real a) {
		debug_check(a);
		a.value = std::fabs(a.value);
		return a;
	}
	friend Real copysign(Real a, Real b) {
		debug_check(a);
		debug_check(b);
		a.value = std::copysign(a.value, b.value);
		return a;
	}
	friend Real expm1(Real a) {
		debug_check(a);
		a.value = std::expm1(a.value);
		return a;
	}
	friend Real exp(Real a) {
		debug_check(a);
		a.value = std::exp(a.value);
		return a;
	}
	friend Real log(Real a) {
		debug_check(a);
		a.value = std::log(a.value);
		return a;
	}
	friend Real cos(Real a) {
		debug_check(a);
		a.value = std::cos(a.value);
		return a;
	}
	friend Real sin(Real a) {
		debug_check(a);
		a.value = std::sin(a.value);
		return a;
	}
	friend Real acos(Real a) {
		debug_check(a);
		range_check( -Real(1), a, Real(1));
		a.value = std::acos(a.value);
		return a;
	}
	friend Real asin(Real a) {
		debug_check(a);
		range_check( -Real(1), a, Real(1));
		a.value = std::asin(a.value);
		return a;
	}
	friend auto lround(Real a) {
		debug_check(a);
		return std::lround(a.value);
	}
	friend Real tgamma(Real a) {
		nonneg_check(a);
		debug_check(a);
		a.value = std::tgamma(a.value);
		return a;
	}
	friend Real sqrt(Real a) {
		nonneg_check(a);
		debug_check(a);
		a.value = std::sqrt(a.value);
		return a;
	}
	friend Real pow(Real a, Real b) {
		nonneg_check(a);
		debug_check(a);
		debug_check(b);
		a.value = std::pow(a.value, b.value);
		return a;
	}
	friend Real pow(Real x, int n) {
		nonneg_check(x);
		debug_check(x);
		if (n < 0) {
			return Real(1) / pow(x, -n);
		} else {
			Real y = Real(1);
			Real xn = x;
			while (n) {
				if (n & 1) {
					y *= xn;
				}
				xn *= xn;
				n >>= 1;
			}
			return y;
		}
	}
	friend Real max(Real a, Real b) {
		debug_check(a);
		debug_check(b);
		a.value = std::max(a.value, b.value);
		return a;
	}
	friend Real min(Real a, Real b) {
		debug_check(a);
		debug_check(b);
		a.value = std::min(a.value, b.value);
		return a;
	}
	friend std::string to_string(Real r) {
		return std::to_string(r.value);
	}
private:
	Type value;
	static void nonneg_check(Real a) {
		if constexpr (CHECK) {
			if (a.value < 0.0) {
				std::string errorString = "FATAL ERROR: Illegal operation on negative number.\n";
				errorString += "Stack trace:\n";
				errorString += std::to_string(std::stacktrace::current());
				std::cout << errorString;
				abort();
//				throw std::invalid_argument(errorString);
			}
		}
	}
	static void zero_check(Real a) {
		if constexpr (CHECK) {
			if (a.value == 0.0) {
				std::string errorString = "FATAL ERROR: Divide by zero\n";
				errorString += "Stack trace:\n";
				errorString += std::to_string(std::stacktrace::current());
				std::cout << errorString;
				abort();
//				throw std::invalid_argument(errorString);
			}
		}
	}
	static void debug_check(Real a) {
		if constexpr (CHECK) {
			if (!std::isfinite(a.value)) {
				std::string errorString = "FATAL ERROR: Operation on NaN\n";
				errorString += "Stack trace:\n";
				errorString += std::to_string(std::stacktrace::current());
				std::cout << errorString;
				abort();
//				throw std::invalid_argument(errorString);
			}
		}
	}
	static void range_check(Real a, Real b, Real c) {
		if constexpr (CHECK) {
			if ((b < a) || (b > c)) {
				std::string errorString = "FATAL ERROR: Range violation\n";
				errorString += "Stack trace:\n";
				errorString += std::to_string(std::stacktrace::current());
				std::cout << errorString;
				abort();
//				throw std::invalid_argument(errorString);
			}
		}
	}
};

}

constexpr Math::Real operator""_Real(long double a) {
	return Math::Real(a);
}

#undef CHECK

#endif /* INCLUDE_REAL_HPP_ */
