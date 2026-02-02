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

#include <functional>
#include <utility>
//	auto const createView = [&](auto &data, auto const &idx) {
//		return (idx | std::views::transform([&](int i) -> int & {
//					return data[i];
//				}));
//	};

template <int... λ>
	requires(nonIncreasing<λ...>())
struct YoungTableau {
	template <std::integral auto... Is>
	constexpr YoungTableau() {
		int i = 0;
		((data_[i++] = Is), ...);
	}
	constexpr YoungTableau() :
		data_{} {
	}
	static constexpr int rowCount() {
		return Λ.count();
	}
	static constexpr int colCount() {
		return Λ[0];
	}
	static constexpr int rowLength(int r) {
		return Λ[r];
	}
	static constexpr int colLength(int c) {
		return conjΛ[c];
	}
	static constexpr auto size() {
		return Λ.size();
	}
	static constexpr int rowPermutationCount() {
		int count = 1;
		for (int i = 0; i < rowCount(); i++) {
			count *= factorial(Λ[i]);
		}
		return count;
	}
	static constexpr int colPermutationCount() {
		int count = 1;
		for (int i = 0; i < colCount(); i++) {
			count *= factorial(conjΛ[i]);
		}
		return count;
	}
	static constexpr int permutationCount() {
		return rowPermutationCount() * colPermutationCount();
	}
	constexpr operator Permutation<size()>() const {
		Permutation<size()> p;
		std::copy_n(data_.begin(), size(), p.begin());
		return p;
	}
	constexpr bool nextRowPermutation() {
		if (colCount() == 1) {
			return false;
		}
		int begin = 0;
		auto const it = data_.begin();
		int i, j, r, end;
		for (r = 0; r < rowCount(); r++) {
			end = begin + Λ[r];
			for (i = end - 2; i >= begin; i--) {
				if (data_[i] < data_[i + 1]) {
					break;
				}
			}
			if (i >= begin) {
				break;
			}
			std::reverse(it + begin, it + end);
			begin = end;
		}
		if (r == rowCount()) {
			return false;
		}
		for (j = end - 1;; j--) {
			if (data_[j] > data_[i]) {
				break;
			}
		}
		std::swap(data_[i], data_[j]);
		std::reverse(it + i + 1, it + end);
		return true;
	}
	constexpr Sign nextColPermutation() {
		if (rowCount() == 1) {
			return 0;
		}
		constexpr auto idx = [] {
			std::array<int, size()> idx;
			int i = 0;
			for (int c = 0; c < colCount(); c++) {
				int colIdx = c;
				for (int r = 0; r < conjΛ[c]; r++) {
					idx[i++] = colIdx;
					colIdx += Λ[r];
				}
			}
			return idx;
		}();
		int i, j, c, end;
		int begin = 0;
		Sign sgn = +1;
		for (c = 0; c < colCount(); c++) {
			int const n = conjΛ[c];
			end = begin + n;
			for (i = end - 2; i >= begin; i--) {
				if (data_[idx[i]] < data_[idx[i + 1]]) {
					break;
				}
			}
			if (i >= begin) {
				break;
			}
			for (int i1 = begin, i2 = end - 1; i1 < i2; i1++, i2--) {
				std::swap(data_[idx[i1]], data_[idx[i2]]);
				sgn = -sgn;
			}
			begin += n;
		}
		if (c == colCount()) {
			return 0;
		}
		for (j = end - 1;; j--) {
			if (data_[idx[j]] > data_[idx[i]]) {
				break;
			}
		}
		std::swap(data_[idx[i]], data_[idx[j]]);
		sgn = -sgn;
		for (int i1 = i + 1, i2 = end - 1; i1 < i2; i1++, i2--) {
			std::swap(data_[idx[i1]], data_[idx[i2]]);
			sgn = -sgn;
		}
		return sgn;
	}
	constexpr bool isSemistandard() const {
		for (int r = 0; r < rowCount(); r++) {
			for (int c = 1; c < rowLength(r); c++) {
				if ((*this)(r, c) < (*this)(r, c - 1)) {
					return false;
				}
			}
		}
		for (int c = 0; c < colCount(); c++) {
			for (int r = 1; r < colLength(c); r++) {
				if ((*this)(r, c) <= (*this)(r - 1, c)) {
					return false;
				}
			}
		}
		return true;
	}
	constexpr Index operator[](int i) const {
		return data_[i];
	}
	constexpr Index &operator[](int i) {
		return data_[i];
	}
	constexpr Index operator()(int r, int c) const {
		return data_[rowIndex[r] + c];
	}
	constexpr Index &operator()(int r, int c) {
		return data_[rowIndex[r] + c];
	}
	constexpr auto begin() const {
		return data_.cbegin();
	}
	constexpr auto end() const {
		return data_.cend();
	}
	constexpr auto begin() {
		return data_.begin();
	}
	constexpr auto end() {
		return data_.end();
	}
	template <std::integral auto D>
	constexpr explicit operator Indices<size(), D>() const {
		return Indices<size(), D>(data_);
	}
	template <std::integral auto D>
	constexpr Index index() const {
		return Index(Indices<size(), D>(data_));
	}
	friend std::ostream &operator<<(std::ostream &os, YoungTableau const &Y) {
		for (int r = 0; r < rowCount(); r++) {
			for (int c = 0; c < rowLength(r); c++) {
				os << Y(r, c) << " ";
			}
			os << std::endl;
		}
		os << std::endl;
		return os;
	}

private:
	static constexpr IntegerPartition<λ...> Λ{};
	static constexpr auto conjΛ = Λ.conj();
	static constexpr auto rowIndex = []() {
		std::array<Index, rowCount()> starts;
		starts[0] = 0;
		for (int i = 0; i < rowCount() - 1; i++) {
			starts[i + 1] = starts[i] + Λ[i];
		}
		return starts;
	}();
	std::array<Index, Λ.size()> data_{};
};

template <typename T>
struct IsYoungTableau {
	static constexpr bool value = false;
};

template <IntegerPartitionType auto Λ>
struct IsYoungTableau<YoungTableau<Λ>> {
	static constexpr bool value = true;
};

template <typename T>
concept YoungTableauType = IsYoungTableau<std::remove_cvref_t<T>>::value;

template <int... λ>
	requires(nonIncreasing<λ...>())
constexpr auto genYoungPermutations() {
	YoungTableau<λ...> Y{};
	auto const rank = Y.size();
	auto const count = Y.permutationCount();
	using SignedPermutation = std::pair<Sign, Permutation<rank>>;
	std::array<SignedPermutation, count> permutations;
	Sign parity = 1;
	std::iota(Y.begin(), Y.end(), 0);
	for (int i = 0; i < count; i++) {
		permutations[i] = SignedPermutation(parity, Y);
		parity *= Y.nextColPermutation();
		if (parity == 0) {
			parity = +1;
			Y.nextRowPermutation();
		}
	}
	return permutations;
}

template <int... λ>
	requires(nonIncreasing<λ...>())
constexpr auto genHookLengths() {
	constexpr IntegerPartition<λ...> Λ{};
	std::array<int, Λ.size()> h;
	constexpr auto conjΛ = Λ.conj();
	int i = 0;
	for (int r = 0; r < Λ.count(); r++) {
		int η = Λ[r] + conjΛ[0] - 1 - r;
		h[i++] = η;
		for (int c = 0; c + 1 < Λ[r]; c++) {
			η += conjΛ[c + 1] - conjΛ[c] - 1;
			h[i++] = η;
		}
	}
	return h;
}

template <int... λ>
	requires(nonIncreasing<λ...>())
constexpr auto genHookLengths(IntegerPartition<λ...> const &Λ) {
	return genHookLengths<λ...>();
}

template <IntegerPartitionType auto Λ, std::integral auto D>
constexpr auto countSemistandardTableau() {
	constexpr auto conjΛ = Λ.conj();
	constexpr auto h = genHookLengths(Λ);
	int n = 1;
	int d = 1;
	int i = 0;
	for (int r = 0; r < Λ.count(); r++) {
		for (int c = 0; c < Λ[r]; c++) {
			n *= D + c - r;
			d *= h[i++];
		}
	}
	return n / d;
}

template <IntegerPartitionType auto Λ, std::integral auto D>
constexpr auto genSemistandardTableau() {
	constexpr int count = countSemistandardTableau<Λ, D>();
	auto Y = createYoungTableau(Λ);
	std::array<decltype(Y), count> tabs;
	if (Y.rowCount() > D) {
		return tabs;
	}
	int i = 0;
	while (true) {
		if (Y.isSemistandard()) {
			tabs[i++] = Y;
		}
		if (i == count) {
			break;
		}
		int j = 0;
		while (++Y[j] >= D) {
			Y[j++] = 0;
		}
	}
	return tabs;
}

template <std::integral auto... λ>
constexpr auto createYoungTableau(IntegerPartition<λ...> const &Λ) {
	return YoungTableau<λ...>{};
}

template <std::integral auto... λ>
constexpr auto genYoungPermutations(IntegerPartition<λ...> const &Λ) {
	return genYoungPermutations<λ...>();
}
