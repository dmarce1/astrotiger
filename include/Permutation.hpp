/*
 * Permutation.hpp
 *
 *  Created on: Dec 26, 2025
 *      Author: dmarce1
 */

#pragma once

#include "Definitions.hpp"
#include "Indices.hpp"
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
struct Permutation : public std::array<Index, N> {
	constexpr Permutation() = default;
	template <typename T>
	constexpr Permutation(std::initializer_list<T> const &list) {
		std::copy(list.begin(), list.end(), this->begin());
	}
	template <std::integral auto M>
	constexpr Permutation(Indices<N, M> const &indices) {
		std::copy(indices.begin(), indices.end(), this->begin());
	}
};

template <int N>
constexpr auto identityPermutation() {
	Permutation<N> p;
	std::iota(p.begin(), p.end(), 0);
	return p;
}

template <int N>
constexpr Permutation<N> operator*(Permutation<N> const &A, Permutation<N> const &B) {
	Permutation<N> C;
	for (int i = 0; i < N; i++) {
		C[i] = A[B[i]];
	}
	return C;
}

template <int N>
constexpr Permutation<N> &operator*=(Permutation<N> &A, Permutation<N> const &B) {
	A = A * B;
	return A;
}

template <int N, Indexed<size_t> Container>
constexpr auto apply(Permutation<N> const &P, Container const &A) {
	Container B = A;
	for (int i = 0; i < N; i++) {
		B[i] = A[P[i]];
	}
	return B;
}

template <int N, int I, int... Is>
constexpr auto deleteAt(Permutation<N> const &A) {
	Permutation<N - 1> B;
	int const delVal = A[I];
	for (int i = 0; i < I; i++) {
		B[i] = (A[i] > delVal) ? A[i] - 1 : A[i];
	}
	for (int i = I; i < N - 1; i++) {
		B[i] = (A[i + 1] > delVal) ? A[i + 1] - 1 : A[i + 1];
	}
	if constexpr (sizeof...(Is) > 0) {
		return B.template deleteAt<Is...>();
	} else {
		return B;
	}
}

template <typename T, int N>
concept PermutationType = std::is_base_of_v<std::array<int, N>, T>;

template <typename Begin, typename End>
constexpr Sign parity(Begin const &begin, End const &end) {
	Sign p = +1;
	for (auto it1 = begin; it1 != end; it1++) {
		for (auto it2 = it1 + 1; it2 != end; it2++) {
			if (*it1 == *it2) {
				return Sign(0);
			} else if (*it1 > *it2) {
				p = -p;
			}
		}
	}
	return p;
}

template <int N>
constexpr auto parity(Permutation<N> const &p) {
	return parity(p.begin(), p.end());
}

template <int N>
constexpr Permutation<N> inverse(Permutation<N> const &P) {
	Permutation<N> iP{};
	for (int i = 0; i < N; i++) {
		iP[P[i]] = i;
	}
	return iP;
}

template <int N>
std::ostream &operator<<(std::ostream &os, Permutation<N> const &p) {
	os << "(";
	os << std::to_string(p[0]);
	for (int i = 1; i < N; i++) {
		os << " ";
		os << std::to_string(p[i]);
	}
	os << ")";
	return os;
}

template <int N>
constexpr bool operator==(Permutation<N> const &A, Permutation<N> const &B) {
	for (int n = 0; n < N; n++) {
		if (A[n] != B[n]) {
			return false;
		}
	}
	return true;
}

template <int N>
constexpr bool operator!=(Permutation<N> const &A, Permutation<N> const &B) {
	return !(A == B);
}

template <std::integral... Args>
constexpr auto createPermutation(Args... args) {
	constexpr int N = sizeof...(Args);
	Permutation<N> P = {args...};
	return P;
}

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
