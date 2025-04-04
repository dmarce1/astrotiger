/*
 * SymmetricGroup.hpp
 *
 *  Created on: Apr 4, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_SYMMETRICGROUP_HPP_
#define INCLUDE_SYMMETRICGROUP_HPP_

#include <bitset>

#include "Numbers.hpp"
#include "Permutation.hpp"

template<size_t N>
constexpr auto genFactorialArray() {
	std::array<size_t, N> factorial;
	factorial[0] = 1;
	for (size_t n = 1; n < N; n++) {
		factorial[n] = n * factorial[n - 1];
	}
	return factorial;
}

template<size_t N>
constexpr size_t factoradicToIndex(Permutation<N> P) {
	static constexpr auto factorial = genFactorialArray<N>();
	size_t index = 0;
	for (size_t i = 0; i < N; i++) {
		size_t const j = N - 1 - i;
		for (size_t k = 0; k < j; k++) {
			if (P[k] > P[j]) {
				P[k]--;
			}
		}
		index += factorial[j] * (P[j] - 1);
	}
	return index;
}

template<size_t N>
using SymmetryPermutation = std::pair<int, Permutation<N>>;

template<size_t N>
struct Symmetries {
	static constexpr size_t M = Math::nFactorial<size_t>(N);
	std::array<bool, M> symmetries;
	std::array<bool, M> antisymmetries;
};

template<size_t N>
constexpr size_t genSymmetryPermutationCount(Symmetries<N> const &S) {
	size_t sum = 0;
	for( size_t m = 0; m < S.symmetries.size(); m++) {
		sum += size_t(S.symmetries[m]);
		sum += size_t(S.antisymmetries[m]);
	}
	return sum;
}

template<size_t N, SymmetryPermutation<N> ...P>
constexpr auto genSymmetries() {
	Symmetries<N> S;
	((((P.first > 0) ? S.symmetries[factoradicToIndex<N>(P.second)] : S.antisymmetries[factoradicToIndex<N>(P.second)]) = true),...);
	return S;
}

template<size_t N, Symmetries<N> S>
constexpr auto genSymmetryPermutations() {
	static constexpr size_t M = genSymmetryPermutationCount(S);
	std::array<SymmetryPermutation<N>, M> permutations;
	size_t j = 0;
	for (size_t i = 0; i < S.symmetries.size(); i++) {
		Permutation<N> thisPermutation = Permutation<N>::identity;
		if (S.symmetries[i] || S.antisymmetries[i]) {
			permutations[j].first = S.antisymmetries[i] ? -1 : +1;
			permutations[j].second = thisPermutation;
			j++;
		}
		thisPermutation = thisPermutation.next();
	}
	return permutations;
}

//template<size_t N>
//struct SymmetricGroup {
//	SymmetricGroup() {
//		for (size_t i = 0; i < N; i++) {
//			P[i] = i + 1;
//		}
//		for (size_t i = 0; i < N - 1; i++) {
//			P[i + 1] = P[i].next();
//		}
//	}
//private:
//	static constexpr size_t Size = Math::nFactorial(N);
//	std::array<Permutation<N>, Size> P;
//};

#endif /* INCLUDE_SYMMETRICGROUP_HPP_ */
