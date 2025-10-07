/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_RATIONAL_HPP_
#define INCLUDE_RATIONAL_HPP_

#include <cmath>
#include <cstdint>
#include <fstream>
#include <numeric>
#include <type_traits>

struct Rational {
	using Type = intmax_t;
	constexpr Rational() :
			n_(0), d_(1) {
	}
	constexpr Rational(Type n) :
			n_(n), d_(1) {
	}
	constexpr Rational(Type n, Type d) :
			n_(n), d_(d) {
	}
	constexpr Rational& operator=(Rational i) {
		n_ = i.n_;
		d_ = i.d_;
		return *this;
	}
	constexpr Rational& operator+=(Rational i) {
		*this = *this + i;
		return *this;
	}
	constexpr Rational& operator-=(Rational i) {
		*this = *this - i;
		return *this;
	}
	constexpr Rational& operator*=(Rational i) {
		*this = *this * i;
		return *this;
	}
	constexpr Rational& operator/=(Rational i) {
		*this = *this / i;
		return *this;
	}
	constexpr Rational& operator=(Type i) {
		n_ = i;
		d_ = 1;
		return *this;
	}
	constexpr Rational& operator+=(Type i) {
		*this = *this + i;
		return *this;
	}
	constexpr Rational& operator-=(Type i) {
		*this = *this - i;
		return *this;
	}
	constexpr Rational& operator*=(Type i) {
		*this = *this * i;
		return *this;
	}
	constexpr Rational& operator/=(Type i) {
		*this = *this / i;
		return *this;
	}
	friend constexpr Rational operator+(Rational b) {
		return b;
	}
	friend constexpr Rational operator-(Rational a) {
		a.n_ = -a.n_;
		a.normalize();
		return a;
	}
	friend constexpr Rational operator+(Rational b, Rational c) {
		Rational a;
		a.n_ = b.d_ * c.n_ + b.n_ * c.d_;
		a.d_ = b.d_ * c.d_;
		a.normalize();
		return a;
	}
	friend constexpr Rational operator+(Rational a, Type b) {
		a.n_ += a.d_ * b;
		a.normalize();
		return a;
	}
	friend constexpr Rational operator+(Type a, Rational b) {
		b.n_ += b.d_ * a;
		b.normalize();
		return b;
	}
	friend constexpr Rational operator-(Rational a, Rational b) {
		return a + -b;
	}
	friend constexpr Rational operator-(Rational a, Type b) {
		return a + -b;
	}
	friend constexpr Rational operator-(Type a, Rational b) {
		return a + -b;
	}
	friend constexpr Rational operator*(Rational a, Rational b) {
		a.n_ *= b.n_;
		a.d_ *= b.d_;
		a.normalize();
		return a;
	}
	friend constexpr Rational operator*(Rational a, Type b) {
		a.n_ *= b;
		a.normalize();
		return a;
	}
	template<typename Real, std::enable_if<std::is_floating_point<Real>::value, int>::type = 0>
	friend constexpr Real operator+(Rational a, Real b) {
		return (a.n_ + a.d_ * b) / a.d_;
	}
	template<typename Real, std::enable_if<std::is_floating_point<Real>::value, int>::type = 0>
	friend constexpr Real operator+(Real a, Rational b) {
		return b + a;
	}
	template<typename Real, std::enable_if<std::is_floating_point<Real>::value, int>::type = 0>
	friend constexpr Real operator-(Rational a, Real b) {
		return a + (-b);
	}
	template<typename Real, std::enable_if<std::is_floating_point<Real>::value, int>::type = 0>
	friend constexpr Real operator-(Real a, Rational b) {
		return a + (-b);
	}
	template<typename Real, std::enable_if<std::is_floating_point<Real>::value, int>::type = 0>
	friend constexpr Real operator*(Rational a, Real b) {
		b *= a.n_;
		b /= a.d_;
		return b;
	}
	template<typename Real, std::enable_if<std::is_floating_point<Real>::value, int>::type = 0>
	friend constexpr Real operator*(Real a, Rational b) {
		return b * a;
	}
	template<typename Real, std::enable_if<std::is_floating_point<Real>::value, int>::type = 0>
	friend constexpr Real operator/(Rational a, Real b) {
		return a.n_ / (b * a.d_);
	}
	template<typename Real, std::enable_if<std::is_floating_point<Real>::value, int>::type = 0>
	friend constexpr Real operator/(Real a, Rational b) {
		a *= b.d_;
		a /= b.n_;
		return a;
	}
	friend constexpr Rational operator*(Type a, Rational b) {
		return b * a;
	}
	friend constexpr Rational operator/(Rational a, Rational b) {
		a.n_ *= b.d_;
		a.d_ *= b.n_;
		a.normalize();
		return a;
	}
	friend constexpr Rational operator/(Type b, Rational a) {
		std::swap(a.n_, a.d_);
		a.n_ *= b;
		a.normalize();
		return a;
	}
	friend constexpr Rational operator/(Rational a, Type b) {
		a.d_ *= b;
		a.normalize();
		return a;
	}
	friend constexpr bool operator==(Rational a, Rational b) {
		return (a.n_ * b.d_) == (a.d_ * b.n_);
	}
	friend constexpr bool operator!=(Rational a, Rational b) {
		return (a.n_ * b.d_) != (a.d_ * b.n_);
	}
	friend constexpr bool operator<(Rational a, Rational b) {
		return (a.n_ * b.d_) < (a.d_ * b.n_);
	}
	friend constexpr bool operator>(Rational a, Rational b) {
		return (a.n_ * b.d_) > (a.d_ * b.n_);
	}
	friend constexpr bool operator<=(Rational a, Rational b) {
		return (a.n_ * b.d_) <= (a.d_ * b.n_);
	}
	friend constexpr bool operator>=(Rational a, Rational b) {
		return (a.n_ * b.d_) >= (a.d_ * b.n_);
	}
	constexpr operator double() const {
		return double(n_) / double(d_);
	}
	friend std::ostream& operator<<(std::ostream &os, Rational f) {
		os << f.n_;
		if (f.d_ != 1) {
			os << "/" << f.d_;
		}
		return os;
	}
	constexpr void normalize() {
		if (n_ == 0) {
			d_ = 1;
		} else {
			auto const d = std::copysign(std::gcd(n_, d_), d_);
			n_ /= d;
			d_ /= d;
		}
	}
	Type n_ = 0;
	Type d_ = 1;
};

#endif /* INCLUDE_RATIONAL_HPP_ */
