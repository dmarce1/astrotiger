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

template<typename Type, int ndim>
using Tensor2 = SquareMatrix<Type, ndim, true>;

template<typename Type, int ndim>
using Tensor1 = Vector<Type, ndim>;

template<typename Type, int ndim>
using Tensor0 = Type;

template<typename Type, int ndim>
inline auto space2spaceTime(Tensor0<Type, ndim> const &s, Tensor1<Type, ndim> const &v, Tensor2<Type, ndim> const &t) {
	Tensor2<Type, ndim + 1> τ;
	for (int j = 0; j < ndim; j++) {
		for (int k = 0; k < ndim; k++) {
			τ(j, k) = t(k, j);
		}
		τ(j, ndim) = τ(ndim, j) = v[j];
	}
	τ(ndim, ndim) = s;
	return τ;
}

template<typename Type, int ndim>
inline auto space2spaceTime(Tensor0<Type, ndim> const &s, Tensor1<Type, ndim> const &v) {
	Tensor1<Type, ndim + 1> V;
	for (int j = 0; j < ndim; j++) {
		V[j] = v[j];
	}
	V[ndim] = s;
	return V;
}

template<typename Type, int ndim>
inline auto spaceTime2Tensor0(Tensor2<Type, ndim> const &τ) {
	Tensor0<Type, ndim - 1> s;
	s = τ(ndim - 1, ndim - 1);
	return s;
}

template<typename Type, int ndim>
inline auto spaceTime2Tensor1(Tensor2<Type, ndim> const &τ) {
	Tensor1<Type, ndim - 1> v;
	for (int j = 0; j + 1 < ndim; j++) {
		v[j] = τ(j, ndim - 1);
	}
	return v;
}

template<typename Type, int ndim>
inline auto spaceTime2Tensor2(Tensor2<Type, ndim> const &τ) {
	Tensor2<Type, ndim - 1> t;
	for (int j = 0; j + 1 < ndim; j++) {
		t[j] = τ(j, ndim - 1);
	}
	return t;
}

#include "conserved.hpp"
#include "primitive.hpp"

#endif /* INCLUDE_SRHD_SRHD_HPP_ */
