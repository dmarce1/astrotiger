/*
 * YoungTableau.hpp
 *
 *  Created on: Dec 27, 2025
 *      Author: dmarce1
 */

#pragma once

#include "Definitions.hpp"
#include "Indices.hpp"
#include "IntegerPartition.hpp"
#include "Math.hpp"
#include "Matrix.hpp"
#include "Rational.hpp"
#include "WeightedPermutation.hpp"

#include <utility>

template <int... λ>
	requires(sizeof...(λ) > 0 && nonIncreasing<λ...>())
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
	constexpr YoungTableau(IntegerPartition<λ...> const &) {
		std::iota(base_type::begin(), base_type::end(), 0);
	}
	constexpr YoungTableau(base_type const &init) :
		base_type(init) {
	}
	static constexpr int rowPermutationCount() {
		return (1 * ... * factorial<λ>());
	}
	static constexpr int colPermutationCount() {
		return YoungTableau{}.conj().rowPermutationCount();
	}
	static constexpr int permutationCount() {
		return rowPermutationCount() * colPermutationCount();
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
	constexpr YoungTableau nextColPermutation() const {
		return (conj().nextRowPermutation()).conj();
	}
	constexpr bool isStandard() const {
		for (int r = 0; r + 1 < nRow; r++) {
			for (int c = 0; c < nCol[r]; c++) {
				if ((*this)(r, c) >= (*this)(r + 1, c)) {
					return false;
				}
			}
		}
		for (int r = 0; r < nRow; r++) {
			for (int c = 0; c + 1 < nCol[r]; c++) {
				if ((*this)(r, c) >= (*this)(r, c + 1)) {
					return false;
				}
			}
		}
		return true;
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
	static constexpr int standardCount() {
		int den = 1;
		for (int r1 = 0; r1 < nRow; r1++) {
			for (int c = 0; c < nCol[r1]; c++) {
				int h = nCol[r1] - c;
				for (int r2 = r1 + 1; r2 < nRow; r2++) {
					if (c >= nCol[r2]) {
						break;
					}
					h++;
				}
				den *= h;
			}
		}
		return factorial(size()) / den;
	}
	static constexpr int semistandardCount(int dimCount) {
		int den = 1;
		int num = 1;
		for (int r1 = 0; r1 < nRow; r1++) {
			for (int c = 0; c < nCol[r1]; c++) {
				int h = nCol[r1] - c;
				num *= dimCount + c - r1;
				for (int r2 = r1 + 1; r2 < nRow; r2++) {
					if (c >= nCol[r2]) {
						break;
					}
					h++;
				}
				den *= h;
			}
		}
		return num / den;
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
	}
	static constexpr int nRow = sizeof...(λ);
	static constexpr auto nCol = std::array<int, nRow>{λ...};
	static constexpr auto start_ = genStarts();
};

template <int... λ>
	requires (sizeof...(λ) > 0 && nonIncreasing<λ...>())
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
