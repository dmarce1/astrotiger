/******************************************************************************
 Copyright (C) 2026  Dominic C. Marcello
 *******************************************************************************/
#pragma once

#include "Permutation.hpp"
#include <array>
#include <bitset>
#include <utility>

template <typename S, typename E, int N>
using LinearCombination = std::array<std::pair<S, E>, N>;

template <typename S, typename E, size_t N>
std::ostream &operator<<(std::ostream &os, std::array<std::pair<S, E>, N> const &linco) {
	constexpr int R = E().size();
	for (unsigned i = 0; i < N; i++) {
		auto co = linco[i].first;
		if (co > S(0)) {
			if (i) {
				os << " + ";
			}
		} else {
			co = -co;
			os << " - ";
		}
		if (co != S(1)) {
			os << co << "*";
		}
		os << Permutation<R>(linco[i].second);
	}
	return os;
}

template <typename S, typename E, int N, LinearCombination<S, E, N> linCo>
consteval auto countUnique() {
	int count = 0;
	for (int n1 = 0; n1 < N; n1++) {
		int unique = 1;
		for (int n2 = n1 + 1; n2 < N; n2++) {
			if (linCo[n1].second == linCo[n2].second) {
				unique = 0;
				break;
			}
		}
		count += unique;
	}
	return count;
}

template <typename S, typename E, int N, LinearCombination<S, E, N> linCo1>
consteval auto combineDuplicates() {
	constexpr auto M = countUnique<S, E, N, linCo1>();
	std::bitset<N> visited{};
	LinearCombination<S, E, M> linCo2;
	for (int n1 = 0, m = 0; n1 < N; n1++) {
		if (visited[n1]) {
			continue;
		}
		linCo2[m] = linCo1[n1];
		for (int n2 = n1 + 1; n2 < N; n2++) {
			if (linCo1[n1].second == linCo1[n2].second) {
				linCo2[m].first += linCo1[n2].first;
				visited[n2] = true;
			}
		}
		m++;
	}
	return linCo2;
}

template <typename S, typename E, int N, LinearCombination<S, E, N> linCo>
consteval auto countZeros() {
	int count = 0;
	for (int n = 0; n < N; n++) {
		if (linCo[n].first == S(0)) {
			count++;
		}
	}
	return count;
}

template <typename S, typename E, int N, LinearCombination<S, E, N> linCo1>
consteval auto removeZeros() {
	constexpr auto M = N - countZeros<S, E, N, linCo1>();
	LinearCombination<S, E, M> linCo2;
	for (int n = 0, m = 0; n < N; n++) {
		if (linCo1[n].first != S(0)) {
			linCo2[m++] = linCo1[n];
		}
	}
	return linCo2;
}

template <typename S, typename E, int N, LinearCombination<S, E, N> linCo>
consteval auto canonicalize() {
	constexpr auto uniques = combineDuplicates<S, E, N, linCo>();
	return removeZeros<S, E, uniques.size(), uniques>();
}

template <typename S, typename E, E V>
consteval auto ele2linco() {
	LinearCombination<S, E, 1> linCo;
	linCo[0].first = S(1);
	linCo[0].second = V;
	return linCo;
}

template <typename S, typename E, int N, LinearCombination<S, E, N> linCo1>
consteval auto negate() {
	LinearCombination<S, E, N> linCo2 = linCo1;
	for (int n = 0; n < N; n++) {
		linCo2[n].first = -linCo2[n].first;
	}
	return linCo2;
}

namespace detail {
template <typename S, typename E, int N1, int N2, LinearCombination<S, E, N1> linCo1, LinearCombination<S, E, N2> linCo2>
consteval auto plus() {
	constexpr auto N3 = N1 + N2;
	LinearCombination<S, E, N3> linCo3;
	for (int n1 = 0; n1 < N1; n1++) {
		linCo3[n1] = linCo1[n1];
	}
	for (int n2 = 0; n2 < N2; n2++) {
		linCo3[N1 + n2] = linCo2[n2];
	}
	return linCo3;
}
} // namespace detail
  // namespace detail

template <typename S, typename E, int N1, int N2, LinearCombination<S, E, N1> linCo1, LinearCombination<S, E, N2> linCo2>
consteval auto plus() {
	constexpr auto linCo3 = detail::plus<S, E, N1, N2, linCo1, linCo2>();
	return canonicalize<S, E, N1 + N2, linCo3>();
}

template <typename S, typename E, int N1, int N2, LinearCombination<S, E, N1> linCo1, LinearCombination<S, E, N2> linCo2>
consteval auto minus() {
	return plus<S, E, N1, N2, linCo1, negate<S, E, N2, linCo2>()>();
}

template <typename S, typename E, int N, LinearCombination<S, E, N> linCo1>
consteval auto multiplies(S const &scale) {
	LinearCombination<S, E, N> linCo2 = linCo1;
	for (int n = 0; n < N; n++) {
		linCo2[n].first *= scale;
	}
	return linCo2;
}

namespace detail {
template <typename S, typename E1, typename E2, int N1, int N2, LinearCombination<S, E1, N1> linCo1, LinearCombination<S, E2, N2> linCo2>
consteval auto rawDirectProduct() {
	using E3 = decltype(concatenate(E1{}, E2{}));
	constexpr int N3 = N1 + N2;
	LinearCombination<S, E3, N3> linCo3;
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			linCo3[n1 * N2 + n2].first = linCo1[n1].first * linCo2[n2].first;
			linCo3[n1 * N2 + n2].second = linCo1[n1].second * linCo2[n2].second;
		}
	}
	return linCo3;
}
} // namespace detail

template <typename S, typename E1, typename E2, int N1, int N2, LinearCombination<S, E1, N1> linCo1, LinearCombination<S, E2, N2> linCo2>
consteval auto directProduct() {
	using E3 = decltype(concatenate(E1{}, E2{}));
	constexpr int N3 = N1 + N2;
	constexpr auto linCo3 = detail::rawDirectProduct<S, E1, E2, N1, N2, linCo1, linCo2>();
	return canonicalize<S, E3, N3, linCo3>();
}
