/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_SRHD_SRHD_HPP_
#define INCLUDE_SRHD_SRHD_HPP_

#include "type_traits.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"



#pragma GCC diagnostic pop

template<int>
struct GasConserved;

template<int>
struct GasPrimitive;

template<typename Type, int NDIM>
using Tensor2 = SquareMatrix<Type, NDIM, SymmetryType::symmetric>;

template<typename Type, int NDIM>
using Tensor1 = Vector<Type, NDIM>;

template<typename Type, int NDIM>
using Tensor0 = Type;

template<typename Type, int NDIM>
inline auto space2spaceTime(Tensor0<Type, NDIM> const &s, Tensor1<Type, NDIM> const &v, Tensor2<Type, NDIM> const &t) {
	Tensor2<Type, NDIM + 1> τ;
	for (int j = 0; j < NDIM; j++) {
		for (int k = 0; k < NDIM; k++) {
			τ(j, k) = t(k, j);
		}
		τ(j, NDIM) = τ(NDIM, j) = v[j];
	}
	τ(NDIM, NDIM) = s;
	return τ;
}

template<typename Type, int NDIM>
inline auto space2spaceTime(Tensor0<Type, NDIM> const &s, Tensor1<Type, NDIM> const &v) {
	Tensor1<Type, NDIM + 1> V;
	for (int j = 0; j < NDIM; j++) {
		V[j] = v[j];
	}
	V[NDIM] = s;
	return V;
}

template<typename Type, int NDIM>
inline auto spaceTime2Tensor0(Tensor2<Type, NDIM> const &τ) {
	Tensor0<Type, NDIM - 1> s;
	s = τ(NDIM - 1, NDIM - 1);
	return s;
}

template<typename Type, int NDIM>
inline auto spaceTime2Tensor1(Tensor2<Type, NDIM> const &τ) {
	Tensor1<Type, NDIM - 1> v;
	for (int j = 0; j + 1 < NDIM; j++) {
		v[j] = τ(j, NDIM - 1);
	}
	return v;
}

template<typename Type, int NDIM>
inline auto spaceTime2Tensor2(Tensor2<Type, NDIM> const &τ) {
	Tensor2<Type, NDIM - 1> t;
	for (int j = 0; j + 1 < NDIM; j++) {
		t[j] = τ(j, NDIM - 1);
	}
	return t;
}

#include "conserved.hpp"
#include "primitive.hpp"

#endif /* INCLUDE_SRHD_SRHD_HPP_ */
