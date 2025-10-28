/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
*******************************************************************************/

#ifndef INCLUDE_MATRIX_FWD_HPP_
#define INCLUDE_MATRIX_FWD_HPP_



enum class SymmetryType : int {
	antisymmetric = -1, asymmetric = 0, symmetric = +1
};

template<typename, int>
struct Vector;

template<typename >
struct IsVector {
	static constexpr bool value = false;
};

template<typename Type, int size>
struct IsVector<Vector<Type, size>> {
	static constexpr bool value = true;
};

template<typename Type, int rowCount, int columnCount, SymmetryType symmetry = SymmetryType::asymmetric>
struct Matrix;

template<typename Type, int count, SymmetryType symmetry = SymmetryType::asymmetric>
using SquareMatrix = Matrix<Type, count, count, symmetry>;

template<typename >
struct IsMatrix {
	static constexpr bool value = false;
};

template<typename Type, int rowCount, int columnCount, SymmetryType symmetry>
struct IsMatrix<Matrix<Type, rowCount, columnCount, symmetry>> {
	static constexpr bool value = true;
};





#endif /* INCLUDE_MATRIX_FWD_HPP_ */
