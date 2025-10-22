/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include <array>
#include <initializer_list>
#include <type_traits>

#include "matrix_fwd.hpp"

template<typename Type, int Size>
struct Vector {
	using value_type = Type;
	constexpr Vector() = default;
	constexpr Vector(Vector const&) = default;
	constexpr Vector(Vector&&) = default;
	template<typename OtherType>
	constexpr Vector(Vector<OtherType, Size> const &other) {
		this->operator=(other);
	}
	constexpr Vector(std::array<Type, Size> const &initList) :
			ν_(initList) {
	}
	constexpr Vector(std::initializer_list<Type> initList) {
		std::copy(initList.begin(), initList.end(), ν_.begin());
	}
	constexpr Vector(Type const &init) {
		for (int i = 0; i < Size; i++) {
			ν_[i] = init;
		}
	}
	constexpr Vector& operator=(Vector const&) = default;
	constexpr Vector& operator=(Vector&&) = default;
	template<typename OtherType>
	constexpr Vector& operator=(Vector<OtherType, Size> const &other) {
		for (int n = 0; n < Size; n++) {
			ν_[n] = Type(other[n]);
		}
		return *this;
	}
	constexpr Type const& operator[](int i) const {
		return ν_[i];
	}
	constexpr Type& operator[](int i) {
		return ν_[i];
	}
	constexpr Vector& operator+=(Vector const &a) {
		*this = *this + a;
		return *this;
	}
	constexpr Vector& operator-=(Vector const &a) {
		*this = *this - a;
		return *this;
	}
	constexpr Vector& operator*=(Type const &a) {
		*this = *this * a;
		return *this;
	}
	constexpr Vector& operator/=(Type const &a) {
		*this = *this / a;
		return *this;
	}
	constexpr bool operator==(Vector const &a) const {
		return ν_ == a.ν_;
	}
	constexpr bool operator!=(Vector const &a) const {
		return !(*this == a);
	}
	friend constexpr Vector operator+(Vector const &a) {
		return a;
	}
	friend constexpr Vector operator-(Vector a) {
		for (int n = 0; n < Size; n++) {
			a[n] = -a[n];
		}
		return a;
	}
	friend constexpr Vector operator+(Vector a, Vector const &b) {
		for (int n = 0; n < Size; n++) {
			a[n] += b[n];
		}
		return a;
	}
	friend constexpr Vector operator-(Vector a, Vector const &b) {
		for (int n = 0; n < Size; n++) {
			a[n] -= b[n];
		}
		return a;
	}
	template<typename OtherType, typename = std::enable_if_t<!IsVector<OtherType>::value && !IsMatrix<OtherType>::value>>
	friend constexpr auto operator*(Vector const &b, OtherType const &c) {
		using ReturnType = decltype(Type() * OtherType());
		Vector<ReturnType, Size> a;
		for (int n = 0; n < Size; ++n) {
			a[n] = b[n] * c;
		}
		return a;
	}
	template<typename OtherType, typename = std::enable_if_t<!IsVector<OtherType>::value && !IsMatrix<OtherType>::value>>
	friend constexpr auto operator*(OtherType const &c, Vector const &b) {
		return b * c;
	}
	template<typename OtherType, typename = std::enable_if_t<!IsVector<OtherType>::value>>
	friend constexpr auto operator/(Vector const &b, OtherType const &c) {
		using ReturnType = decltype(Type() / OtherType());
		Vector<ReturnType, Size> a;
		for (int n = 0; n < Size; n++) {
			a[n] = b[n] / c;
		}
		return a;
	}
	template<typename OtherType, typename ReturnType = decltype(Type() * OtherType())>
	constexpr auto dot(Vector<OtherType, Size> const &b) const {
		Vector const &a = *this;
		ReturnType sum(0);
		for (int n = 0; n < Size; n++) {
			sum += a[n] * b[n];
		}
		return sum;
	}
	static constexpr size_t size() {
		return Size;
	}
	constexpr Type& front() {
		return ν_.front();
	}
	constexpr Type& back() {
		return ν_.back();
	}
	constexpr auto begin() {
		return ν_.begin();
	}
	constexpr auto end() {
		return ν_.end();
	}
	constexpr Type front() const {
		return ν_.front();
	}
	constexpr Type back() const {
		return ν_.back();
	}
	constexpr auto begin() const {
		return ν_.begin();
	}
	constexpr auto end() const {
		return ν_.cend();
	}
private:
	std::array<Type, Size> ν_;

};

template<typename Type, int Size>
Type abs(Vector<Type, Size> const &v) {
	return sqrt(v.dot(v));
}

template<typename Type, int Size>
Type normalize(Vector<Type, Size> const &v) {
	auto const v1 = abs(v);
	if (v1) {
		return v / v1;
	} else {
		return Type(0);
	}
}

template<typename Type, int Size>
static constexpr Vector<Type, Size> unitVector(int i) {
	Vector<Type, Size> ζ;
	for (int j = 0; j < Size; j++) {
		ζ[j] = Type(i == j);
	}
	return ζ;
}

template<typename Type, int Size>
std::ostream& operator<<(std::ostream& os, Vector<Type, Size> const& v) {
	os << '[';
	for (int i = 0; i < Size; i++) {
		if (i > 0) {
			os << ", ";
		}
		os << v[i];
	}
	os << ']';
	return os;
}
