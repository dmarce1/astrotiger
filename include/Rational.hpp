/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once


#include <cmath>
#include <cstdint>
#include <numeric>
#include <ostream>
#include <type_traits>

struct Rational {
	using Type = int64_t;
	constexpr Rational() :
			n_(0), d_(1) {
		normalize();
	}
	constexpr Rational(Type n) :
			n_(n), d_(1) {
		normalize();
	}
	constexpr Rational(Type n, Type d) :
			n_(n), d_(d) {
		normalize();
	}
	constexpr Rational& operator=(Rational i) {
		n_ = i.n_;
		d_ = i.d_;
		normalize();
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
//	constexpr Rational& operator+=(Type i) {
//		*this = *this + i;
//		return *this;
//	}
//	constexpr Rational& operator-=(Type i) {
//		*this = *this - i;
//		return *this;
//	}
//	constexpr Rational& operator*=(Type i) {
//		*this = *this * i;
//		return *this;
//	}
//	constexpr Rational& operator/=(Type i) {
//		*this = *this / i;
//		return *this;
//	}
	friend constexpr Rational operator+(Rational b) {
		return b;
	}
	friend constexpr Rational operator-(Rational a) {
		a.n_ = -a.n_;
		a.normalize();
		return a;
	}
	friend constexpr Rational operator+(Rational b, Rational c) {
	    Type const g = std::gcd(b.d_, c.d_);
	    Type const bd1 = b.d_ / g;
	    Type const cd1 = c.d_ / g;

	    Rational a;
	    a.n_ = b.n_ * cd1 + c.n_ * bd1;
	    a.d_ = bd1 * c.d_;            // == bd1 * (g*cd1) == lcm(b.d_, c.d_)
	    a.normalize();
	    return a;
	}
//	friend constexpr Rational operator+(Rational a, Type b) {
//		a.n_ += a.d_ * b;
//		a.normalize();
//		return a;
//	}
//	friend constexpr Rational operator+(Type a, Rational b) {
//		b.n_ += b.d_ * a;
//		b.normalize();
//		return b;
//	}
	friend constexpr Rational operator-(Rational a, Rational b) {
		return a + (-b);
	}
//	friend constexpr Rational operator-(Rational a, Type b) {
//		return a + (-b);
//	}
//	friend constexpr Rational operator-(Type a, Rational b) {
//		return a + (-b);
//	}
	friend constexpr Rational operator*(Rational a, Rational b) {
		Type const g1 = std::gcd(a.n_, b.d_);
		a.n_ /= g1;
		b.d_ /= g1;
		Type const g2 = std::gcd(b.n_, a.d_);
		b.n_ /= g2;
		a.d_ /= g2;
		a.n_ *= b.n_;
		a.d_ *= b.d_;
		a.normalize();
		return a;
	}
//	friend constexpr Rational operator*(Rational a, Type b) {
//		a.n_ *= b;
//		a.normalize();
//		return a;
//	}
//	template<typename Real, std::enable_if<std::is_floating_point<Real>::value, int>::type = 0>
//	friend constexpr Real operator+(Rational a, Real b) {
//		return (a.n_ + a.d_ * b) / a.d_;
//	}
//	template<typename Real, std::enable_if<std::is_floating_point<Real>::value, int>::type = 0>
//	friend constexpr Real operator+(Real a, Rational b) {
//		return b + a;
//	}
//	template<typename Real, std::enable_if<std::is_floating_point<Real>::value, int>::type = 0>
//	friend constexpr Real operator-(Rational a, Real b) {
//		return a + (-b);
//	}
//	template<typename Real, std::enable_if<std::is_floating_point<Real>::value, int>::type = 0>
//	friend constexpr Real operator-(Real a, Rational b) {
//		return a + (-b);
//	}
//	template<typename Real, std::enable_if<std::is_floating_point<Real>::value, int>::type = 0>
//	friend constexpr Real operator*(Rational a, Real b) {
//		b *= a.n_;
//		b /= a.d_;
//		return b;
//	}
//	template<typename Real, std::enable_if<std::is_floating_point<Real>::value, int>::type = 0>
//	friend constexpr Real operator*(Real a, Rational b) {
//		return b * a;
//	}
//	template<typename Real, std::enable_if<std::is_floating_point<Real>::value, int>::type = 0>
//	friend constexpr Real operator/(Rational a, Real b) {
//		return a.n_ / (b * a.d_);
//	}
//	template<typename Real, std::enable_if<std::is_floating_point<Real>::value, int>::type = 0>
//	friend constexpr Real operator/(Real a, Rational b) {
//		a *= b.d_;
//		a /= b.n_;
//		return a;
//	}
//	friend constexpr Rational operator*(Type a, Rational b) {
//		return b * a;
//	}
	friend constexpr Rational operator/(Rational a, Rational b) {
	    Type const g1 = std::gcd(a.n_, b.n_);
	    a.n_ /= g1;
	    b.n_ /= g1;

	    Type const g2 = std::gcd(b.d_, a.d_);
	    b.d_ /= g2;
	    a.d_ /= g2;

	    a.n_ *= b.d_;
	    a.d_ *= b.n_;

	    a.normalize();
	    return a;
	}
//	friend constexpr Rational operator/(Type b, Rational a) {
//		std::swap(a.n_, a.d_);
//		a.n_ *= b;
//		a.normalize();
//		return a;
//	}
//	friend constexpr Rational operator/(Rational a, Type b) {
//		a.d_ *= b;
//		a.normalize();
//		return a;
//	}
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
	constexpr void normalize() {
		if (n_ == 0) {
			d_ = 1;
		} else {
			auto d = icopysign(std::gcd(n_, d_), d_);
			n_ /= d;
			d_ /= d;
		}
	}
	constexpr Type denominator() const {
		return d_;
	}
	constexpr Type numerator() const {
		return n_;
	}
	friend constexpr Rational abs(Rational v) {
		v.n_ = (v.n_ >= 0) ? v.n_ : -v.n_;
		v.normalize();
		return v;
	}
	friend std::ostream& operator<<(std::ostream &os, Rational f) {
		os << f.n_;
		if (f.d_ != 1) {
			os << "/" << f.d_;
		}
		return os;
	}
	Type n_ = 0;
	Type d_ = 1;
};
