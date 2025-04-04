/*
 * Vector.hpp
 *
 *  Created on: Dec 11, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_VECTOR_HPP_
#define INCLUDE_VECTOR_HPP_

#include <array>
#include <cstddef>

#include "Matrix.hpp"
#include "Utilities.hpp"
#include "Zobrist.hpp"

namespace Math {

template<typename Type, int Ndim>
using Vector = Matrix<Type, Ndim, 1>;

template<typename Type, int Ndim>
Vector<Type, Ndim> vectorMaximum(Vector<Type, Ndim> const &A, Vector<Type, Ndim> const &B) {
	Vector<Type, Ndim> result;
	for (int n = 0; n < Ndim; n++) {
		result[n] = A[n] < B[n] ? B[n] : A[n];
	}
	return result;
}

template<typename Type, int Ndim>
Type vectorMaximum(Vector<Type, Ndim> const &A) {
	Type maxA = A[0];
	for (int n = 1; n < Ndim; n++) {
		if (maxA < A[n]) {
			maxA = A[n];
		}
	}
	return maxA;
}

template<typename Type, int Ndim>
Vector<Type, Ndim> vectorMinimum(Vector<Type, Ndim> const &A, Vector<Type, Ndim> const &B) {
	Vector<Type, Ndim> result;
	for (int n = 0; n < Ndim; n++) {
		result[n] = A[n] < B[n] ? A[n] : B[n];
	}
	return result;
}

template<typename Type, int Ndim>
Type vectorMinimum(Vector<Type, Ndim> const &A) {
	Type minA = A[0];
	for (int n = 1; n < Ndim; n++) {
		if (minA > A[n]) {
			minA = A[n];
		}
	}
	return minA;
}

template<typename Type, int Ndim>
Vector<Type, Ndim> operator*(Type const &scalar, Vector<Type, Ndim> const &vector) {
	return vector * scalar;
}

template<typename Type, int Ndim>
Vector<Type, Ndim> zeroVector() {
	Vector<Type, Ndim> result;
	std::fill(result.begin(), result.end(), Type(0));
	return result;
}

template<typename Type, int Ndim>
Vector<Type, Ndim> unitVector(int dim) {
	Vector<Type, Ndim> result;
	std::fill(result.begin(), result.end(), Type(0));
	result[dim] = Type(1);
	return result;
}

template<typename Type, int Ndim>
Vector<Type, Ndim> constantVector(Type constant) {
	Vector<Type, Ndim> result;
	std::fill(result.begin(), result.end(), constant);
	return result;
}

template<typename Type, int Ndim>
Type vectorDotProduct(Vector<Type, Ndim> const &A, Vector<Type, Ndim> const &B) {
	Type product = Type(0);
	for (int d = 0; d < Ndim; d++) {
		product += A[d] * B[d];
	}
	return product;
}

template<typename Type, int Ndim>
Type vectorNorm(Vector<Type, Ndim> const &A) {
	return vectorDotProduct(A, A);
}

template<typename Type, int Ndim>
Type vectorMagnitude(Vector<Type, Ndim> const &A) {
	return sqrt(vectorNorm(A));
}

template<typename Type, int Ndim>
Type vectorSum(Vector<Type, Ndim> const &A) {
	Type sum = Type(0);
	for (int d = 0; d < Ndim; d++) {
		sum += A[d];
	}
	return sum;
}

template<typename Type, int Ndim>
Type vectorProduct(Vector<Type, Ndim> const &A) {
	Type product = Type(1);
	for (int d = 0; d < Ndim; d++) {
		product *= A[d];
	}
	return product;
}

template<typename RType, typename IType, int Ndim>
constexpr std::enable_if<std::is_integral_v<RType>, IType>::type Pow(Vector<RType, Ndim> const &x, Vector<IType, Ndim> n) {
	RType xn = RType(1);
	for (int d = 0; d < Ndim; d++) {
		RType xm = x[d];
		while (n[d]) {
			if (n[d] & 1) {
				xn *= xm;
			}
			xm *= xm;
			n[d] >>= 1;
		}
	}
	return xn;
}

template<typename Type, int Ndim>
SquareMatrix<Type, Ndim> vectorTensorProduct(Vector<Type, Ndim> const &A, Vector<Type, Ndim> const &B) {
	SquareMatrix<Type, Ndim> C;
	for (int n = 0; n < Ndim; n++) {
		for (int k = 0; k < Ndim; k++) {
			C[n, k] = A[n] * B[k];
		}
	}
	return C;
}

template<typename Type, int Ndim>
Vector<Type, Ndim> vectorIota(Type init = Type(0)) {
	static constexpr Type one = Type(1);
	Vector<Type, Ndim> I;
	I[0] = init;
	for (int k = 0; k < Ndim - 1; k++) {
		I[k + 1] = I[k] + one;
	}
	return I;
}

template<typename T, int D, T MaxT>
struct vectorHashKey {
	size_t operator()(Vector<T, D> const &vector) const {
		static ZobristGenerator<T> const zobristGenerator(MaxT);
		return zobristGenerator(vector.begin(), vector.end());
	}
};

}
#endif /* INCLUDE_VECTOR_HPP_ */
