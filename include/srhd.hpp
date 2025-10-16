/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_SRHD_SRHD_HPP_
#define INCLUDE_SRHD_SRHD_HPP_

#include "type_traits.hpp"

template<typename, int>
struct GasConserved;

template<typename, int>
struct GasPrimitive;

template<typename Type, int dimensionCount>
using Tensor2 = SquareMatrix<Type, dimensionCount, SymmetryType::symmetric>;

template<typename Type, int dimensionCount>
using Tensor1 = Vector<Type, dimensionCount>;

template<typename Type, int dimensionCount>
using Tensor0 = Type;

template<typename Type, int dimensionCount>
inline auto space2spaceTime(Tensor0<Type, dimensionCount> const &s, Tensor1<Type, dimensionCount> const &v, Tensor2<Type, dimensionCount> const &t) {
	Tensor2<Type, dimensionCount + 1> τ;
	for (int j = 0; j < dimensionCount; j++) {
		for (int k = 0; k < dimensionCount; k++) {
			τ(j, k) = t(k, j);
		}
		τ(j, dimensionCount) = τ(dimensionCount, j) = v[j];
	}
	τ(dimensionCount, dimensionCount) = s;
	return τ;
}

template<typename Type, int dimensionCount>
inline auto space2spaceTime(Tensor0<Type, dimensionCount> const &s, Tensor1<Type, dimensionCount> const &v) {
	Tensor1<Type, dimensionCount + 1> V;
	for (int j = 0; j < dimensionCount; j++) {
		V[j] = v[j];
	}
	V[dimensionCount] = s;
	return V;
}

template<typename Type, int dimensionCount>
inline auto spaceTime2Tensor0(Tensor2<Type, dimensionCount> const &τ) {
	Tensor0<Type, dimensionCount - 1> s;
	s = τ(dimensionCount - 1, dimensionCount - 1);
	return s;
}

template<typename Type, int dimensionCount>
inline auto spaceTime2Tensor1(Tensor2<Type, dimensionCount> const &τ) {
	Tensor1<Type, dimensionCount - 1> v;
	for (int j = 0; j + 1 < dimensionCount; j++) {
		v[j] = τ(j, dimensionCount - 1);
	}
	return v;
}

template<typename Type, int dimensionCount>
inline auto spaceTime2Tensor2(Tensor2<Type, dimensionCount> const &τ) {
	Tensor2<Type, dimensionCount - 1> t;
	for (int j = 0; j + 1 < dimensionCount; j++) {
		t[j] = τ(j, dimensionCount - 1);
	}
	return t;
}

#include "conserved.hpp"
#include "primitive.hpp"

#endif /* INCLUDE_SRHD_SRHD_HPP_ */
