/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/
#pragma once


#include "util.hpp"
#include "vector.hpp"

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
	template<template<typename, int> typename Container>
	constexpr Interval& operator+=(Container<T, D> const &C) {
		*this = *this + C;
		return *this;
	}
	template<template<typename, int> typename Container>
	constexpr Interval& operator-=(Container<T, D> const &C) {
		*this = *this / C;
		return *this;
	}
	constexpr Interval& operator*=(T const &value) {
		*this = *this * value;
		return *this;
	}
	constexpr Interval& operator/=(T const &value) {
		*this = *this / value;
		return *this;
	}
	template<template<typename, int> typename Container>
	constexpr Interval operator+(Container<T, D> const &C) const {
		Interval i;
		for (int d = 0; d < D; d++) {
			i.A_[d] = A_ + C[d];
			i.B_[d] = B_ + C[d];
		}
		return i;
	}
	template<template<typename, int> typename Container>
	constexpr Interval operator-(Container<T, D> const &C) const {
		Interval i;
		for (int d = 0; d < D; d++) {
			i.A_[d] = A_ - C[d];
			i.B_[d] = B_ - C[d];
		}
		return i;
	}
	constexpr Interval operator*(T const &value) const {
		Interval i;
		for (int d = 0; d < D; d++) {
			i.A_[d] *= value;
			i.B_[d] *= value;
		}
		return i;
	}
	constexpr Interval operator/(T const &value) const {
		Interval i;
		for (int d = 0; d < D; d++) {
			i.A_[d] /= value;
			i.B_[d] /= value;
		}
		return i;
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
	constexpr bool contains(Interval const &other) const {
		for (int k = 0; k < D; k++) {
			if (other.A_[k] < A_[k]) {
				return false;
			}
			if (other.B_[k] > B_[k]) {
				return false;
			}
		}
		return true;
	}
	constexpr T span(int k) const {
		return B_[k] - A_[k];
	}
	constexpr Vector<T, D> span() const {
		Vector<T, D> s;
		for (int j = 0; j < D; j++) {
			s[j] = span(j);
		}
		return s;
	}
	constexpr std::pair<T, T> sub(int k) const {
		return std::pair<T, T>(A_[k], B_[k]);
	}
	constexpr T stride(int k) const {
		return size(k + 1, D);
	}
	constexpr Vector<T, D> stride() const {
		Vector<T, D> s;
		for (int j = 0; j < D; j++) {
			s[j] = stride(j);
		}
		return s;
	}
	constexpr size_t size(int b = 0, int e = D) const {
		size_t sz = 1;
		for (int k = b; k < e; k++) {
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
		for (int k = 1; k < D; k++) {
			i *= span(k);
			i += C[k];
		}
		return i;
	}
	constexpr Interval swapDims(int dim1, int dim2) {
		Interval I = *this;
		std::swap(I.A_[dim1], I.A_[dim2]);
		std::swap(I.B_[dim1], I.B_[dim2]);
		return I;
	}
	constexpr T begin(int d) const {
		return A_[d];
	}
	constexpr T end(int d) const {
		return B_[d];
	}
	constexpr Vector<T, D> begin() const {
		return A_;
	}
	constexpr Vector<T, D> end() const {
		return B_;
	}
	constexpr Interval& setBegin(int d, T const &v) {
		A_[d] = v;
		return *this;
	}
	constexpr Interval& setEnd(int d, T const &v) {
		B_[d] = v;
		return *this;
	}
	constexpr Interval& setRange(int d, T const &b, T const &e) {
		setBegin(d, b);
		setEnd(d, e);
		return *this;
	}
	friend constexpr Interval operator*(T const &value, Interval const &i) {
		return i * value;
	}
private:
	Vector<T, D> A_;
	Vector<T, D> B_;
};

