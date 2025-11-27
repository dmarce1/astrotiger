/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "Definitions.hpp"
#include "Integer.hpp"

#include <cmath>

template<int width, typename Real = double>
struct FixedPrecision {
	using real_type = Real;
	constexpr FixedPrecision(FixedPrecision const&) = default;
	constexpr FixedPrecision(FixedPrecision&&) = default;
	constexpr FixedPrecision& operator=(FixedPrecision const&) = default;
	constexpr FixedPrecision& operator=(FixedPrecision&&) = default;
	template<typename Arc>
	void serialize(Arc &&arc, unsigned) {
		arc & bits_;
	}
	constexpr FixedPrecision() :
			bits_(0) {
	}
	constexpr FixedPrecision(Real value) :
			bits_(value * real2int) {
	}
	constexpr FixedPrecision& operator=(Real value) {
		bits_ = value * real2int;
		return *this;
	}
	constexpr explicit operator Real() const {
		return bits_ * int2real;
	}
	constexpr FixedPrecision& operator+=(FixedPrecision const &other) {
		bits_ += other.bits_;
		return *this;
	}
	constexpr FixedPrecision& operator-=(FixedPrecision const &other) {
		bits_ -= other.bits_;
		return *this;
	}
	constexpr FixedPrecision& operator+=(Real value) {
		bits_ *= value;
		return *this;
	}
	constexpr FixedPrecision& operator/=(Real value) {
		bits_ *= one / value;
		return *this;
	}
	constexpr FixedPrecision operator+(FixedPrecision const &other) const {
		return FixedPrecision(bits_ + other.bits_);
	}
	constexpr FixedPrecision operator-(FixedPrecision const &other) const {
		return FixedPrecision(bits_ - other.bits_);
	}
	constexpr FixedPrecision operator*(Real value) const {
		return FixedPrecision(Type(std::round(bits_ * value)));
	}
	constexpr FixedPrecision operator/(Real value) const {
		return *this * (one / value);
	}
	constexpr bool operator==(FixedPrecision const &other) const {
		return bool(bits_ == other.bits_);
	}
	constexpr bool operator!=(FixedPrecision const &other) const {
		return bool(bits_ != other.bits_);
	}
	constexpr bool operator>(FixedPrecision const &other) const {
		return bool(bits_ > other.bits_);
	}
	constexpr bool operator<(FixedPrecision const &other) const {
		return bool(bits_ < other.bits_);
	}
	constexpr bool operator>=(FixedPrecision const &other) const {
		return bool(bits_ >= other.bits_);
	}
	constexpr bool operator<=(FixedPrecision const &other) const {
		return bool(bits_ <= other.bits_);
	}
	friend constexpr FixedPrecision nexttoward(FixedPrecision from, FixedPrecision to) {
		if (to > from) {
			from.bits_ += 1;
		} else if (to < from) {
			from.bits_ -= 1;
		}
		return from;
	}
	friend constexpr FixedPrecision operator*(Real valueA, FixedPrecision const &valueB) {
		return valueA * valueB;
	}
	friend constexpr FixedPrecision max(FixedPrecision const &valueA, FixedPrecision const &valueB) {
		return (valueA >= valueB) ? valueA : valueB;
	}
	friend constexpr FixedPrecision min(FixedPrecision const &valueA, FixedPrecision const &valueB) {
		return (valueA <= valueB) ? valueA : valueB;
	}
	friend constexpr FixedPrecision midpoint(FixedPrecision const &valueA, FixedPrecision const &valueB) {
		Type const &a = valueA.bits_;
		Type const &b = valueA.bits_;
		return FixedPrecision(((a >> 1) + (b >> 1)) + (a & b & 1));

	}
	static constexpr FixedPrecision min() {
		return FixedPrecision(Type(0));
	}
	static constexpr FixedPrecision max() {
		return FixedPrecision((Type(1) << Type(width)) - Type(1));
	}
	static constexpr FixedPrecision epsilon() {
		return FixedPrecision(1);
	}
private:
	static constexpr Real real2int = std::ldexp(Real(1), +width);
	static constexpr Real int2real = std::ldexp(Real(1), -width);
	static constexpr Real one = Real(1);
	using Type = UnsignedInteger<width>;
	Type bits_;
	constexpr FixedPrecision(Type b) :
			bits_(b) {
	}
};

template<typename >
struct numeric_limits;

template<int width, typename Real>
struct numeric_limits<FixedPrecision<width, Real>> {
	using value_type = FixedPrecision<width, Real>;
	static constexpr value_type max() {
		return value_type::max();
	}
	static constexpr value_type min() {
		return value_type::min();
	}
	static constexpr value_type epsilon() {
		return value_type::epsilon();
	}
};
