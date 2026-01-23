/*
 * WeightedPermutation.hpp
 *
 *  Created on: Jan 8, 2026
 *      Author: dmarce1
 */

#pragma once

#include "Permutation.hpp"
#include "SymmetricGroup.hpp"

#include <type_traits>

template <typename T, int N>
struct WeightedPermutation {
	static constexpr int length() {
		return N;
	}
	constexpr WeightedPermutation() :
		p_{}, c_(1) {
	}
	constexpr WeightedPermutation(Permutation<N> const &p) :
		p_(p), c_(1) {
	}
	constexpr WeightedPermutation(std::pair<T, int> const &p) :
		p_(Sn[p.second]), c_(p.first) {
	}
	template <typename Q>
	struct CoefficientSelector {
		using type = Q;
	};

	template <typename Q>
		requires ConvertsToLiteral<Q>::value
	struct CoefficientSelector<Q> {
		using type = typename Q::LiteralType;
	};

	using CoefficientType = typename CoefficientSelector<T>::type;
	using LiteralType = std::pair<CoefficientType, typename Permutation<N>::LiteralType>;
	constexpr LiteralType literal() const {
		LiteralType rc;
		rc.second = p_.literal();
		if constexpr (ConvertsToLiteral<T>::value) {
			rc.first = c_.lteral();
		} else {
			rc.first = c_;
		}
		return rc;
	}
	constexpr operator LiteralType() const {
		return literal();
	}
	constexpr WeightedPermutation(LiteralType const &lit) :
		p_(lit.second), c_(lit.first) {
	}
	constexpr WeightedPermutation(T c, Permutation<N> const &p) :
		p_(p), c_(c) {
	}
	constexpr WeightedPermutation(WeightedPermutation const &) = default;
	constexpr WeightedPermutation(WeightedPermutation &&) = default;
	constexpr operator std::pair<T, int>() const {
		return std::pair(c_, Sn.find(p_));
	}
	template <int... Is>
	constexpr auto deleteAt() const {
		return WeightedPermutation<T, N - sizeof...(Is)>(c_, p_.template deleteAt<Is...>());
	}
	constexpr WeightedPermutation &operator=(WeightedPermutation const &) = default;
	constexpr WeightedPermutation &operator=(WeightedPermutation &&) = default;
	constexpr WeightedPermutation &operator*=(T const &factor) {
		*this = *this * factor;
		return *this;
	}
	constexpr WeightedPermutation &operator/=(T const &factor) {
		*this = *this / factor;
		return *this;
	}
	constexpr WeightedPermutation operator+() const {
		return *this;
	}
	constexpr WeightedPermutation operator-() const {
		return WeightedPermutation(-c_, p_);
	}
	constexpr WeightedPermutation operator*(T const &factor) const {
		return WeightedPermutation(factor * c_, p_);
	}
	constexpr WeightedPermutation operator/(T const &factor) const {
		return WeightedPermutation(c_ / factor, p_);
	}
	friend constexpr WeightedPermutation operator*(T const &factor, WeightedPermutation &term) {
		return term * factor;
	}
	template <int M>
	constexpr WeightedPermutation<T, N + M> operator*(WeightedPermutation<T, M> const &B) const {
		return WeightedPermutation<T, N + M>(c_ * B.c_, concatenate(p_, B.p_));
	}
	friend std::ostream &operator<<(std::ostream &os, WeightedPermutation const &p) {
		if (p.c_ != T(1)) {
			os << p.c_ << "*";
		}
		os << p.p_;
		return os;
	}
	constexpr T weight() const {
		return c_;
	}
	constexpr auto permutation() const {
		return p_;
	}
	template <int I, int J>
	friend constexpr auto transpose(WeightedPermutation<T, N> wP) {
		std::swap(wP.p_[I], wP.p_[J]);
		return wP;
	}
	template <typename, int>
	friend class WeightedPermutation;

private:
	static constexpr SymmetricGroup<N> Sn{};
	Permutation<N> p_;
	T c_;
};
