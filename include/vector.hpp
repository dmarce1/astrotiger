/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include <array>
#include <initializer_list>
#include <type_traits>

template<typename, int>
struct Vector;

template<typename >
struct IsVector {
	static constexpr bool value = false;
};

template<typename Type, int size>
struct IsVector<Vector<Type, size>> {
	static constexpr bool value = true;
};

template<typename Type, int Size>
struct Vector {
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
	template<typename OtherType, typename = std::enable_if_t<!IsVector<OtherType>::value>>
	friend constexpr auto operator*(Vector const &b, OtherType const &c) {
		using ReturnType = decltype(Type() * OtherType());
		Vector<ReturnType, Size> a;
		for (int n = 0; n < Size; ++n) {
			a[n] = b[n] * c;
		}
		return a;
	}
	template<typename OtherType, typename = std::enable_if_t<!IsVector<OtherType>::value>>
	friend constexpr auto operator*(OtherType const &c, Vector const &b) {
		return b * c;
	}
	template<typename OtherType, typename = std::enable_if_t<!IsVector<OtherType>::value>>
	friend constexpr auto operator/(Vector const &b, OtherType const &c) {
		using ReturnType = decltype(Type() / OtherType());
		auto const ic = 1 / c;
		Vector<ReturnType, Size> a;
		for (int n = 0; n < Size; n++) {
			a[n] = b[n] * ic;
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
private:
	std::array<Type, Size> ν_;

};

