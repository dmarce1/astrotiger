/*
 * YoungTableau.hpp
 *
 *  Created on: Dec 27, 2025
 *      Author: dmarce1
 */

#pragma once

#include "IntegerPartition.hpp"
#include "Math.hpp"
#include "Matrix.hpp"
#include "Multidices.hpp"
#include "Permutation.hpp"

#include <utility>

template <int... λ>
	requires(nonIncreasing<λ...>())
struct YoungTableau : public Permutation<(0 + ... + λ)> {
	using base_type = Permutation<(0 + ... + λ)>;
	constexpr operator base_type() const {
		return static_cast<base_type const &>(*this);
	}
	static constexpr auto size() {
		return (0 + ... + λ);
	}
	constexpr int operator()(int i, int j) const {
		return base_type::operator[](start_[i] + j);
	}
	constexpr int &operator()(int i, int j) {
		return base_type::operator[](start_[i] + j);
	}
	constexpr YoungTableau() {
		std::iota(base_type::begin(), base_type::end(), 0);
	}
	static constexpr int rowPermutationCount() {
		return (1 * ... * factorial<λ>());
	}
	static constexpr int columnPermutationCount() {
		return YoungTableau{}.conj().rowPermutationCount();
	}
	static constexpr int permutationCount() {
		return rowPermutationCount() * columnPermutationCount();
	}
	constexpr YoungTableau nextRowPermutation() const {
		YoungTableau next = *this;
		int *b, *e;
		int r = 0;
		do {
			b = &(next(r, 0));
			e = &(next(r, nCol[r]));
			r++;
		} while (!std::next_permutation(b, e) && (r < nRow));
		return next;
	}
	constexpr YoungTableau nextColumnPermutation() const {
		return (conj().nextRowPermutation()).conj();
	}
	friend std::ostream &operator<<(std::ostream &os, YoungTableau<λ...> const &yt) {
		auto printHLine = [&](int row) {
			os << "+";
			for (int m = 0; m < nCol[std::max(row - 1, 0)]; m++) {
				os << "---+";
			}
			os << '\n';
		};
		for (int n = 0; n < nRow; n++) {
			printHLine(n);
			os << "|";
			for (int m = 0; m < nCol[n]; m++) {
				os << " " << yt(n, m) << " |";
			}
			os << '\n';
		}
		printHLine(nRow - 1);

		return os;
	}
	constexpr bool operator==(YoungTableau const &other) const {
		return base_type::operator==(static_cast<base_type const &>(other));
	}
	constexpr bool operator!=(YoungTableau const &other) const {
		return base_type::operator!=(static_cast<base_type const &>(other));
	}
	constexpr auto conj() const;

private:
	static constexpr auto genStarts() {
		std::array<int, nRow + 1> starts{};
		for (int n = 0; n < nRow; n++) {
			starts[n + 1] = starts[n] + nCol[n];
		}
		return starts;
	};
	static constexpr int nRow = sizeof...(λ);
	static constexpr auto nCol = std::array<int, nRow>{λ...};
	static constexpr auto start_ = genStarts();
};

template <int... λ>
	requires(nonIncreasing<λ...>())
constexpr auto YoungTableau<λ...>::conj() const {
	auto const create = []<int... δ>(std::integer_sequence<int, δ...>) {
		constexpr auto conjPart = IntegerPartition<λ...>::conj();
		return YoungTableau<conjPart[δ]...>{};
	};
	constexpr auto seq = std::make_integer_sequence<int, nCol[0]>{};
	auto conjTab = create(seq);
	for (int r = 0; r < nRow; r++) {
		for (int c = 0; c < nCol[r]; c++) {
			conjTab(c, r) = (*this)(r, c);
		}
	}
	return conjTab;
}

template <typename T, int N>
constexpr Matrix<T, N, N> permutationMatrix(Permutation<N> const &p) {
	Matrix<T, N, N> P{};
	for (int n = 0; n < N; n++) {
		P(n, p[n]) = T(1);
	}
	return P;
}

namespace detail {

template <int... λ>
consteval auto YoungSymmetrizerPermutations() {
	constexpr int R = (0 + ... + λ);
	YoungTableau<λ...> υ{};
	constexpr int rowCount = υ.rowPermutationCount();
	constexpr int colCount = υ.columnPermutationCount();
	constexpr int size = rowCount * colCount;
	std::array<std::pair<int8_t, Permutation<R>>, size> sum;
	int i = 0;
	YoungTableau<λ...> σ{};
	do {
		auto στ = σ;
		do {
			sum[i++] = std::pair<int8_t, Permutation<R>>(στ.parity() / σ.parity(), στ);
			στ = στ.nextColumnPermutation();
		} while (στ != σ);
		σ = σ.nextRowPermutation();
	} while (σ != υ);
	return sum;
}

template <int D, int... λ>
consteval auto YoungSymmetrizerMatrix() {
	using Type = double;
	constexpr int R = (0 + ... + λ);
	constexpr int N = ipow(D, R);
	Matrix<Type, N, N> Ω{};
	auto const sum = YoungSymmetrizerPermutations<λ...>();
	using itype = Multidices<R, D>;
	for (auto const &q : sum) {
		for (auto r = itype::ibegin(); r != itype::iend(); r++) {
			auto const c = q.second.apply(r);
			Ω(Index(r), Index(c)) += Type(q.first);
		}
	}
	return Ω;
}
} // namespace detail

template <int D, int... λ>
consteval auto YoungSymmetrizer() {
	using Type = double;
	constexpr auto zero = Type(0);
	constexpr auto one = Type(1);
	constexpr auto A = detail::YoungSymmetrizerMatrix<D, λ...>();
	constexpr auto Alit = A.literal();
	constexpr int N = A.rowCount();
	constexpr auto reducedA = rankReduce<Type, N, N, Alit>();
	constexpr int M = reducedA.rowCount();
	auto B = pseudoinverse(reducedA);
	for (int m = 0; m < M; m++) {
		Type factor = zero;
		for (int n = 0; n < N; n++) {
			if (B(n, m) != zero) {
				factor = std::max(factor, one / std::abs(B(n, m)));
			}
		}
		for (int n = 0; n < N; n++) {
			B(n, m) *= factor;
		}
	}
	return B;
}
