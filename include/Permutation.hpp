/*
 * Permutation.hpp
 *
 *  Created on: Dec 26, 2025
 *      Author: dmarce1
 */

#pragma once

#include "Definitions.hpp"
#include "Math.hpp"

#include <algorithm>
#include <array>
#include <bitset>
#include <initializer_list>
#include <iostream>
#include <numeric>
#include <vector>

struct Cycle : public std::vector<int> {
	friend std::ostream &operator<<(std::ostream &os, Cycle const &p) {
		os << "(";
		os << std::to_string(p[0]);
		for (size_t i = 1; i < p.size(); i++) {
			os << " ";
			os << std::to_string(p[i]);
		}
		os << ")";
		return os;
	}
};

template <int N>
struct Permutation : public std::array<int, N> {
	using base_type = std::array<int, N>;
	using parity_type = int8_t;
	constexpr Permutation() {
		std::iota(base_type::begin(), base_type::end(), 0);
	}
	template <std::integral Int>
	constexpr Permutation(std::initializer_list<Int> const &list) {
		int i = 0;
		for (int v : list) {
			base_type::operator[](i++) = v;
		}
	}
	constexpr Permutation operator*=(Permutation const &B) const {
		*this = *this * B;
		return *this;
	}
	constexpr Permutation operator*(Permutation const &B) const {
		auto const &A = *this;
		Permutation C;
		for (int i = 0; i < N; i++) {
			C[i] = A[B[i]];
		}
		return C;
	}
	template <Indexed<size_t> Container>
	constexpr auto apply(Container const &A) const {
		auto const &P = *this;
		Container B = A;
		for (int i = 0; i < N; i++) {
			B[i] = A[P[i]];
		}
		return B;
	}
	constexpr parity_type parity() const {
		auto const &P = *this;
		parity_type p = +1;
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				if (P[i] > P[j]) {
					p = -p;
				}
			}
		}
		return p;
	}
	auto cycles() const {
		std::bitset<N> visited{};
		std::vector<Cycle> cycles{};
		auto const &P = *this;
		for (int n = 0; n < N; n++) {
			if (visited[n]) {
				continue;
			}
			Cycle cycle{};
			cycle.push_back(n);
			visited[n] = true;
			for (int k = P[n]; k != n; k = P[k]) {
				cycle.push_back(k);
				visited[k] = true;
			}
			cycles.push_back(cycle);
		}
		return cycles;
	}
	friend constexpr Permutation inverse(Permutation const &P) {
		Permutation iP{};
		for (int i = 0; i < N; i++) {
			iP[P[i]] = i;
		}
		return iP;
	}
	friend std::ostream &operator<<(std::ostream &os, Permutation const &p) {
		os << "(";
		os << std::to_string(p[0]);
		for (int i = 1; i < N; i++) {
			os << " ";
			os << std::to_string(p[i]);
		}
		os << ")";
		return os;
	}
	constexpr bool operator==(Permutation<N> const &other) const {
		for (int n = 0; n < N; n++) {
			if (base_type::operator[](n) != other[n]) {
				return false;
			}
		}
		return true;
	}
	constexpr bool operator!=(Permutation<N> const &other) const {
		return !(*this == other);
	}
};

template <int N1, int N2>
constexpr auto concatenate(Permutation<N1> const &P1, Permutation<N2> const &P2) {
	constexpr int N = N1 + N2;
	Permutation<N> P;
	for (int n = 0; n < N1; n++) {
		P[n] = P1[n];
	}
	for (int n = 0; n < N2; n++) {
		P[n + N1] = P2[n] + N1;
	}
	return P;
}

template <typename T, int N>
struct PermutationTerm {
	constexpr PermutationTerm(Permutation<N> const &p) {
		p_ = p;
		c_ = 1.0;
	}
	constexpr PermutationTerm(T c, Permutation<N> const &p) {
		p_ = p;
		c_ = 1.0;
	}
	constexpr PermutationTerm() = default;
	constexpr PermutationTerm(PermutationTerm const &) = default;
	constexpr PermutationTerm(PermutationTerm &&) = default;
	constexpr PermutationTerm &operator=(PermutationTerm const &) = default;
	constexpr PermutationTerm &operator=(PermutationTerm &&) = default;
	constexpr PermutationTerm &operator*=(T const &factor) {
		*this = *this * factor;
		return *this;
	}
	constexpr PermutationTerm &operator/=(T const &factor) {
		*this = *this / factor;
		return *this;
	}
	constexpr PermutationTerm operator+() const {
		return *this;
	}
	constexpr PermutationTerm operator-() const {
		return PermutationTerm(-c_, p_);
	}
	constexpr PermutationTerm operator*(T const &factor) const {
		return PermutationTerm(factor * c_, p_);
	}
	constexpr PermutationTerm operator/(T const &factor) const {
		return PermutationTerm(factor / c_, p_);
	}
	friend constexpr PermutationTerm operator*(T const &factor, PermutationTerm &term) {
		return term * factor;
	}

private:
	double c_;
	Permutation<N> p_;
};
// template <typename T, int N, int M>
// struct PermutationSum {
//	constexpr PermutationSum() = default;
//	constexpr PermutationSum(PermutationSum const &) = default;
//	constexpr PermutationSum(PermutationSum &&) = default;
//	constexpr PermutationSum &operator=(PermutationSum const &) = default;
//	constexpr PermutationSum &operator=(PermutationSum &&) = default;
//	constexpr PermutationSum operator+() const {
//		return *this;
//	}
//	constexpr PermutationSum operator-() const {
//		PermutationSum result;
//		for (int m = 0; m < M; m++) {
//			result.ps_[m] = -ps_[m];
//		}
//		return result;
//	}
//	constexpr PermutationSum operator*(T const &factor) const {
//		PermutationSum result;
//		for (int m = 0; m < M; m++) {
//			result.ps_[m] = factor * ps_[m];
//		}
//		return result;
//	}
//	constexpr PermutationSum operator/(T const &factor) const {
//		return *this * (T(1) / factor);
//	}
//	template <int L>
//	constexpr PermutationSum operator*(PermutationSum<T, N, L> const &other) const {
//		PermutationSum result;
//		for (int m = 0; m < M; m++) {
//			result.ps_[m] = factor * ps_[m];
//		}
//		return result;
//	}
//
// private:
//	std::array<PermutationTerm<T, N>, M> ps_;
// };
