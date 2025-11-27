/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "vector.hpp"

template<typename Type>
struct LeftRightState;

template<typename Type>
struct IsLeftRightState {
	static constexpr bool value = false;
};

template<typename Type>
struct IsLeftRightState<LeftRightState<Type>> {
	static constexpr bool value = true;
};

template<typename Type>
struct LeftRightState {
	static constexpr bool isVectorType = IsVector<Type>::value;
	constexpr LeftRightState() = default;
	constexpr LeftRightState(LeftRightState const&) = default;
	constexpr LeftRightState(LeftRightState&&) = default;
	constexpr LeftRightState(Type const &init) {
		L = R = init;
	}
	constexpr LeftRightState(Type const &l, Type const &r) {
		L = l;
		R = r;
	}
	constexpr LeftRightState(Type &&l, Type &&r) :
			L(std::move(l)), R(std::move(r)) {
	}
	constexpr LeftRightState& operator=(LeftRightState const&) = default;
	constexpr LeftRightState& operator=(LeftRightState&&) = default;
	constexpr LeftRightState& operator=(Type const &init) {
		L = R = init;
		return *this;
	}
	constexpr LeftRightState& operator+=(LeftRightState const &other) {
		*this = *this + other;
		return *this;
	}
	constexpr LeftRightState& operator-=(LeftRightState const &other) {
		*this = *this - other;
		return *this;
	}
	constexpr LeftRightState& operator*=(LeftRightState const &other) {
		*this = *this * other;
		return *this;
	}
	constexpr LeftRightState& operator/=(LeftRightState const &other) {
		*this = *this / other;
		return *this;
	}
	constexpr auto operator+() const {
		return LeftRightState(+L, +R);
	}
	constexpr auto operator-() const {
		return LeftRightState(-L, -R);
	}
	constexpr auto operator+(LeftRightState<Type> const &other) const {
		auto const l = L + other.L;
		auto const r = R + other.R;
		return LeftRightState<std::remove_cv_t<decltype(l)>>(l, r);
	}
	constexpr auto operator-(LeftRightState<Type> const &other) const {
		auto const l = L - other.L;
		auto const r = R - other.R;
		return LeftRightState<std::remove_cv_t<decltype(l)>>(l, r);
	}
	template<typename OtherType, std::enable_if_t<!IsLeftRightState<OtherType>::value, int> = 0>
	constexpr auto operator+(OtherType const &other) const {
		return LeftRightState(L + other, R + other);
	}
	template<typename OtherType, std::enable_if_t<!IsLeftRightState<OtherType>::value, int> = 0>
	constexpr auto operator-(OtherType const &other) const {
		return LeftRightState(L - other, R - other);
	}
	template<typename OtherType, std::enable_if_t<!IsLeftRightState<OtherType>::value, int> = 0>
	friend constexpr auto operator+(OtherType const &A, IsLeftRightState<Type> const &B) {
		return B + A;
	}
	template<typename OtherType, std::enable_if_t<!IsLeftRightState<OtherType>::value, int> = 0>
	friend constexpr auto operator-(OtherType const &A, IsLeftRightState<Type> const &B) {
		return -B + A;
	}
	template<typename OtherType>
	constexpr auto operator*(LeftRightState<OtherType> const &other) const {
		auto const l = L * other.L;
		auto const r = R * other.R;
		return LeftRightState<std::remove_cv_t<decltype(l)>>(l, r);
	}
	template<typename OtherType, std::enable_if_t<!IsLeftRightState<OtherType>::value, int> = 0>
	constexpr auto operator*(OtherType const &other) const {
		auto const l = L * other;
		auto const r = R * other;
		return LeftRightState<std::remove_cv_t<decltype(l)>>(l, r);
	}
	template<typename TypeA, std::enable_if_t<!IsLeftRightState<TypeA>::value && !IsVector<TypeA>::value, int> = 0>
	friend constexpr auto operator*(TypeA const &A, LeftRightState<Type> const &B) {
		return B * A;
	}
	template<typename OtherType>
	constexpr auto operator/(LeftRightState<OtherType> const &other) const {
		auto const l = L / other.L;
		auto const r = R / other.R;
		return LeftRightState<std::remove_cv_t<decltype(l)>>(l, r);
	}
	template<typename TypeA, std::enable_if_t<!IsLeftRightState<TypeA>::value, int> = 0>
	friend constexpr auto operator/(TypeA const &A, LeftRightState<Type> const &B) {
		auto const l = A / B.L;
		auto const r = A / B.R;
		return LeftRightState<std::remove_cv_t<decltype(l)>>(l, r);
	}
	template<typename OtherType, std::enable_if_t<!IsLeftRightState<OtherType>::value, int> = 0>
	constexpr auto operator/(OtherType const &other) const {
		auto const l = L / other;
		auto const r = R / other;
		return LeftRightState<std::remove_cv_t<decltype(l)>>(l, r);
	}
	friend constexpr auto sqrt(LeftRightState<Type> const &lr) {
		auto const l = sqrt(lr.L);
		auto const r = sqrt(lr.R);
		return LeftRightState<std::remove_cv_t<decltype(l)>>(l, r);
	}
	template<typename OtherType, int dimensionCount>
	auto dot(LeftRightState<Vector<OtherType, dimensionCount>> const &other) const {
		LeftRightState<decltype(OtherType() * typename Type::value_type())> rc;
		rc.L = L.dot(other.L);
		rc.R = R.dot(other.R);
		return rc;
	}
	auto operator[](int n) const {
		LeftRightState<typename Type::value_type> rc;
		rc.L = L[n];
		rc.R = R[n];
		return rc;
	}
	constexpr Type operator()(bool f) const {
		return f ? L : R;
	}
	Type L;
	Type R;
};

template<typename T>
using LRS = LeftRightState<T>;


