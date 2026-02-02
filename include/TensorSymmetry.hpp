#pragma once
//
// #include "Indices.hpp"
// #include "Matrix.hpp"
// #include "Permutation.hpp"
// #include "Rational.hpp"
#include "SparseMatrix.hpp"
// #include "SparseVector.hpp"
#include "Young.hpp"
//
// #include <array>
// #include <bitset>
// #include <concepts>
// #include <tuple>
// #include <type_traits>
// #include <utility>
//
// template <IntegerPartitionType auto Λ>
// constexpr auto hookLengths() {
//	std::array<int, Λ.size()> h;
//	constexpr auto conjΛ = Λ.conj();
//	int i = 0;
//	for (int r = 0; r < Λ.count(); r++) {
//		int η = Λ[r] + conjΛ[0] - 1 - r;
//		h[i++] = η;
//		for (int c = 0; c + 1 < Λ[r]; c++) {
//			η += conjΛ[c + 1] - conjΛ[c] - 1;
//			h[i++] = η;
//		}
//	}
//	return h;
//}
//
// template <IntegerPartitionType auto Λ, std::integral auto D>
// constexpr auto countSemistandardTableau() {
//	constexpr auto conjΛ = Λ.conj();
//	constexpr auto h = hookLengths<Λ>();
//	int n = 1;
//	int d = 1;
//	int i = 0;
//	for (int r = 0; r < Λ.count(); r++) {
//		for (int c = 0; c < Λ[r]; c++) {
//			n *= D + c - r;
//			d *= h[i++];
//		}
//	}
//	return n / d;
//}
//
//

struct SymmetryElement {
	Index row;
	Index col;
	double value;
	friend std::ostream &operator<<(std::ostream &os, SymmetryElement const &ele) {
		os << "(" << ele.row << ", " << ele.col << "): " << ele.value;
		return os;
	}
};

template <std::integral auto D, IntegerPartitionType auto Λ>
constexpr auto symmetrizer() {
	constexpr auto Y = createYoungTableau(Λ);
	auto const lambda = [&Y]<int phase>() {
		constexpr auto basis = genSemistandardTableau<Λ, D>();
		constexpr auto ps = genYoungPermutations(Λ);
		constexpr int R = basis.size();
		constexpr int O = Y.size();
		constexpr int N = ipow(D, O);
		using IndexType = Indices<O, D>;
		auto const pCount = Y.permutationCount();
		SparseMatrix<double> V(R, N);
		int i = 0;
		std::bitset<N> isBasis{};
		for (auto const &b : basis) {
			isBasis[b.template index<D>()] = true;
		}
		for (auto ns = IndexType::ibegin(); ns != IndexType ::iend(); ns++) {
			auto const n = Index(ns);
			if (isBasis[n]) {
				SparseVector<double> vRow(N);
				for (auto p : ps) {
					auto const ms = apply<O>(p.second, ns);
					auto const m = Index(ms);
					vRow[m] += p.first;
				}
				for (auto p : ps) {
					auto const ms = apply<O>(p.second, ns);
					auto const m = Index(ms);
					if (vRow[m] > 0) {
						vRow[m] = 1;
					} else if (vRow[m] < 0) {
						vRow[m] = -1;
					}
				}
				V[i++] = std::move(vRow);
			}
		}
		if constexpr (phase == 0) {
			return V.density();
		} else if constexpr (phase == 1) {
			return transpose(V);
		}
	};
	constexpr int eleCount = lambda.template operator()<0>();
	auto const Q = lambda.template operator()<1>();
	std::array<SymmetryElement, eleCount> elements;
	SymmetryElement ele;
	int i = 0;
	for (unsigned r = 0; r < Q.rowCount(); r++) {
		ele.row = r;
		for (auto const &col : Q[r]) {
			ele.col = col.first;
			ele.value = col.second;
			elements[i++] = ele;
		}
	}
	return elements;
}
