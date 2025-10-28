/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/
#pragma once

#include <array>

#include "util.hpp"

template<typename T, int D>
struct Interval { // @formatter:off
	CONSTEXPR_DEFAULT_MEMBERS(Interval); // @formatter:on
	using literal_type = std::pair<std::array<T, D>, std::array<T, D>>;
	constexpr Interval(literal_type const &a) :
			A_(a.first), B_(a.second) {
	}
	template<template<typename, int> typename Container>
	constexpr Interval(Container<T, D> const &a, Container<T, D> const &b) {
		for (int k = 0; k < D; k++) {
			A_[k] = a[k];
			B_[k] = b[k];
		}
	}
	template<template<typename, int> typename Container>
	constexpr Interval(Container<T, D> const &b) {
		A_.fill(T(0));
		for (int k = 0; k < D; k++) {
			B_[k] = b[k];
		}
	}
	constexpr Interval(T const &a, T const &b) {
		A_.fill(a);
		B_.fill(b);
	}
	constexpr Interval(T const &b) {
		A_.fill(T(0));
		B_.fill(b);
	}
	constexpr literal_type literal() const {
		return literal_type(A_, B_);
	}
	template<template<typename, int> typename Container>
	constexpr bool contains(Container<T, D> const &C) const {
		for (int k = 0; k < D; k++) {
			if (C[k] < A_[k]) {
				return false;
			}
			if (C[k] >= B_[k]) {
				return false;
			}
		}
		return true;
	}
	constexpr T span(int k) const {
		return B_[k] - A_[k];
	}
	constexpr std::array<T, D> span() const {
		std::array<T, D> s;
		for (int j = 0; j < D; j++) {
			s[j] = span(j);
		}
		return s;
	}
	constexpr size_t size() const {
		size_t sz = 1;
		for (int k = 0; k < D; k++) {
			auto const s = span(k);
			if (s <= 0) {
				return 0;
			}
			sz *= s;
		}
		return sz;
	}
	constexpr bool null() const {
		return size() <= 0;
	}
	constexpr bool intersects(Interval<T, D> const &I2) {
		return intersection(*this, I2).size() > 0;
	}
	friend constexpr Interval intersection(Interval<T, D> const &I1, Interval<T, D> const &I2) {
		using std::min;
		using std::max;
		Interval<T, D> I3;
		for (int j = 0; j < D; j++) {
			I3.A_[j] = max(I1.A_[j], I2.A_[j]);
			I3.B_[j] = min(I1.B_[j], I2.B_[j]);
		}
		return I3;
	}
	friend constexpr Interval bounding(Interval<T, D> const &I1, Interval<T, D> const &I2) {
		using std::min;
		using std::max;
		Interval<T, D> I3;
		for (int j = 0; j < D; j++) {
			I3.A_[j] = min(I1.A_[j], I2.A_[j]);
			I3.B_[j] = max(I1.B_[j], I2.B_[j]);
		}
		return I3;
	}
	template<template<typename, int> typename Container, std::enable_if_t<std::is_integral_v<T>, int> = 0>
	constexpr T flatten(Container<T, D> const &C) const {
		T i = C[0];
		for(int k = 1; k < D; k++) {
			i *= span(k);
			i += C[k];
		}
		return i;
	}
private:
	std::array<T, D> A_;
	std::array<T, D> B_;
};

