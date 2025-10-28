/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include <algorithm>

#include "math.hpp"

template<int maxIndex, int dimensionCount, typename Integer = int, std::enable_if_t<std::is_integral<Integer>::value, int> = 0>
struct Indices {
	constexpr Indices() = default;
	constexpr Indices(Indices const&) = default;
	constexpr Indices(Indices&&) = default;
	constexpr Indices(std::array<Integer, dimensionCount> const &init) :
			α_(init) {
	}
	constexpr Indices& operator=(Indices const&) = default;
	constexpr Indices& operator=(Indices&&) = default;
	constexpr Integer operator[](int i) const {
		return α_[i];
	}
	constexpr Integer& operator[](int i) {
		return α_[i];
	}
	constexpr Integer abs() const {
		return std::accumulate(α_.begin(), α_.end(), 0);
	}
	constexpr Integer max() const {
		return *std::max_element(α_.begin(), α_.end());
	}
	constexpr operator std::array<Integer, dimensionCount>() const {
		return α_;
	}
	constexpr operator Integer() const {
		Indices const &α = *this;
		Integer a = 0;
		for (int d = 0; d < dimensionCount; d++) {
			a = (maxIndex + 1) * a + α[d];
		}
		return a;
	}
	constexpr Indices& operator++() {
		using std::min;
		Indices &α = *this;
		for (int d = dimensionCount - 1; d >= 0; d--) {
			if (++α[d] <= maxIndex) {
				return α;
			}
			α[d] = 0;
		}
		α = end();
		return α;
	}
	constexpr Indices operator++(int) {
		Indices const old = *this;
		this->operator++();
		return old;
	}
	friend constexpr bool operator<=(Indices const &α, Integer b) {
		for (int d = 0; d < dimensionCount; d++) {
			if (α[d] > b) {
				return false;
			}
		}
		return true;
	}
	friend constexpr bool operator <(Indices const &α, Integer b) {
		return !(α >= b);
	}
	friend constexpr bool operator >(Indices const &α, Integer b) {
		return !(α <= b);
	}
	friend constexpr bool operator >=(Indices const &α, Integer b) {
		for (int d = 0; d < dimensionCount; d++) {
			if (α[d] < b) {
				return false;
			}
		}
		return true;
	}
	friend constexpr bool operator <=(Indices const &α, Indices const &β) {
		for (int d = 0; d < dimensionCount; d++) {
			if (α[d] > β[d]) {
				return false;
			}
		}
		return true;
	}
	friend constexpr bool operator <(Indices const &α, Indices const &β) {
		return !(α >= β);
	}
	friend constexpr bool operator >(Indices const &α, Indices const &β) {
		return !(α <= β);
	}
	friend constexpr bool operator >=(Indices const &α, Indices const &β) {
		return β <= α;
	}
	friend constexpr bool operator ==(Indices const &α, Indices const &β) {
		for (int d = 0; d < dimensionCount; d++) {
			if (α[d] != β[d]) {
				return false;
			}
		}
		return true;
	}
	friend constexpr bool operator!=(Indices const &α, Indices const &β) {
		return !(α == β);
	}
	friend constexpr Indices operator+(Indices const &α, Indices const &β) {
		Indices γ;
		for (int i = 0; i < dimensionCount; i++) {
			γ[i] = α[i] + β[i];
		}
		return γ;
	}
	friend constexpr Indices operator-(Indices const &α, Indices const &β) {
		Indices γ;
		for (int i = 0; i < dimensionCount; i++) {
			γ[i] = α[i] - β[i];
		}
		return γ;
	}
	friend constexpr Indices operator*(Indices const &α, Integer b) {
		Indices γ;
		for (size_t i = 0; i < dimensionCount; i++) {
			γ[i] = α[i] * b;
		}
		return γ;
	}
	template<typename Real, std::enable_if_t<std::is_floating_point<Real>::value, int> = 0>
	friend constexpr Real pow(std::array<Real, dimensionCount> const &x, Indices const &α) {
		Real z(1.0);
		for (int d = 0; d < dimensionCount; d++) {
			for (int n = 0; n < α[d]; n++) {
				z *= x[d];
			}
		}
		return z;
	}
	static constexpr size_t size() {
		return ::pow<dimensionCount>(maxIndex + 1);
	}
	static constexpr Indices zero() {
		Indices ζ;
		ζ.α_.fill(Integer(0));
		return ζ;
	}
	static constexpr Indices unit(size_t i) {
		Indices ζ = zero();
		ζ[i] = 1;
		return ζ;
	}
	template<auto i>
	static constexpr Indices unit() {
		Indices ζ = zero();
		ζ[i] = 1;
		return ζ;
	}
	static constexpr Indices end() {
		Indices ζ;
		ζ.α_.fill(maxIndex + 1);
		return ζ;
	}
	static constexpr Indices one() {
		Indices ζ;
		ζ.fill(+1);
		return ζ;
	}
private:
	std::array<Integer, dimensionCount> α_;
};

