///******************************************************************************
// Copyright (C) 2024  Dominic C. Marcello
// *******************************************************************************/
//
//#pragma once
//
//#include "Definitions.hpp"
//
//#include <algorithm>
//#include <array>
//#include <cstddef>
//#include <cstdint>
//#include <functional>
//#include <numeric>
//#include <ostream>
//#include <type_traits>
//
//template<size_t N>
//using Permutation = std::array<uint8_t, N>;
//
//template<size_t N>
//auto identity() {
//	Permutation<N> p;
//	std::iota(p.begin(), p.end());
//}
//
//template<size_t N1, size_t N2>
//constexpr auto operator+(Permutation<N1> const &A, Permutation<N2> const &B) {
//	Permutation<N1 + N2> C;
//	for (size_t i = 0; i < N1; i++) {
//		C[i] = A[i];
//	}
//	for (size_t i = 0; i < N2; i++) {
//		C[N1 + i] = B[i];
//	}
//	return C;
//}
//
//template<size_t N>
//constexpr auto operator*(Permutation<N> const &A, Permutation<N> const &B) {
//	Permutation<N> C;
//	for (size_t i = 0; i < N; i++) {
//		C[i] = A[B[i]];
//	}
//	return C;
//}
//
//template<size_t N>
//constexpr auto inverse(Permutation<N> const &A) {
//	Permutation<N> C;
//	for (size_t i = 0; i < N; i++) {
//		C[A[i]] = i;
//	}
//	return C;
//}
//
