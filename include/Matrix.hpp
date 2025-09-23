#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <stacktrace>
#include <type_traits>
#include "ContainerArithmetic.hpp"

template<typename Type, int R1, int C1>
struct Matrix;

template<typename, int>
struct ColumnVector;

template<typename, int>
struct RowVector;

template<typename Type, int Ndim>
using SquareMatrix = Matrix<Type, Ndim, Ndim>;

template<typename Type, int Ndim>
struct DiagonalMatrix;

template<typename T>
struct IsMatrix {
	static constexpr bool value = false;
};

template<typename T, int Nr, int Nc>
struct IsMatrix<Matrix<T, Nr, Nc>> {
	static constexpr bool value = true;
};

template<typename T, int N>
struct IsMatrix<ColumnVector<T, N>> {
	static constexpr bool value = true;
};

template<typename T, int N>
struct IsMatrix<RowVector<T, N>> {
	static constexpr bool value = true;
};

template<typename T, int N>
struct IsMatrix<SquareMatrix<T, N>> {
	static constexpr bool value = true;
};

template<typename T, int N>
struct IsMatrix<DiagonalMatrix<T, N>> {
	static constexpr bool value = true;
};

template<typename T>
struct IsSquareMatrix {
	static constexpr bool value = false;
};

template<typename T, int N>
struct IsSquareMatrix<SquareMatrix<T, N>> {
	static constexpr bool value = true;
};

template<typename T1, int R1, int C1>
struct Matrix {
	static constexpr std::size_t size() {
		return C1 * R1;
	}
	static constexpr std::size_t rowCount() {
		return R1;
	}
	static constexpr std::size_t columnCount() {
		return C1;
	}

	constexpr Matrix() {
	}

	constexpr Matrix(std::array<std::array<T1, C1>, R1> const &initList) :
			values { } {
		for (int n = 0; n < R1; n++) {
			for (int m = 0; m < C1; m++) {
				(*this)(n, m) = initList[n][m];
			}
		}
	}

	constexpr Matrix(std::initializer_list<T1> initList) {
		int i = 0;
		for (auto const &v : initList) {
			values[i++] = v;
		}
	}

	constexpr Matrix(std::initializer_list<std::initializer_list<T1>> init) {
		int i = 0;
		for (auto const &row : init) {
			for (auto const &val : row) {
				values[i++] = val;
			}
		}
	}

	constexpr Matrix(T1 const &init) {
		for (std::size_t i = 0; i < size(); i++) {
			values[i] = init;
		}
	}

	constexpr Matrix(Matrix const &other) :
			values(other.values) {
	}
	constexpr Matrix(Matrix &&other) :
			values(std::move(other.values)) {
	}

	constexpr Matrix& operator=(Matrix const &other) {
		values = other.values;
		return *this;
	}
	constexpr Matrix& operator=(Matrix &&other) {
		values = std::move(other.values);
		return *this;
	}

	constexpr T1& operator()(int n, int m) {
		return values[n * C1 + m];
	}
	constexpr T1 const& operator()(int n, int m) const {
		return values[n * C1 + m];
	}

	constexpr Matrix& operator+=(Matrix const &A) {
		*this = *this + A;
		return *this;
	}
	constexpr Matrix& operator-=(Matrix const &A) {
		*this = *this - A;
		return *this;
	}
	template<typename T2, std::enable_if<std::is_same<T1, decltype(T1() * T2())>::value, int> = 0>
	constexpr Matrix& operator*=(T2 const &a) {
		*this = *this * a;
		return *this;
	}
	template<typename T2, std::enable_if<std::is_same<T1, decltype(T1() / T2())>::value, int> = 0>
	constexpr Matrix& operator/=(T2 const &a) {
		*this = *this / a;
		return *this;
	}
	template<typename T2, typename T3 = decltype(T1() * T2()), std::enable_if_t<!IsMatrix<T2>::value, int> = 0>
	constexpr Matrix<T3, R1, C1> operator*(T2 const &a) const {
		Matrix<T3, R1, C1> B;
		for (std::size_t k = 0; k < size(); k++) {
			B.values[k] = a * values[k];
		}
		return B;
	}
	template<typename T2, typename T3 = decltype(T1() / T2())>
	constexpr Matrix<T3, R1, C1> operator/(T2 const &a) const {
		Matrix<T3, R1, C1> B;
		for (std::size_t k = 0; k < size(); k++) {
			B.values[k] = values[k] / a;
		}
		return B;
	}

	constexpr Matrix operator+() const {
		return *this;
	}
	constexpr Matrix operator-() const {
		Matrix B;
		for (std::size_t k = 0; k < size(); k++) {
			B.values[k] = -values[k];
		}
		return B;
	}

	template<typename T2, typename T3 = decltype(T1() + T2())>
	constexpr auto operator+(Matrix<T2, R1, C1> const &A) const {
		Matrix<T3, R1, C1> B;
		for (std::size_t k = 0; k < size(); k++) {
			B.values[k] = values[k] + A.values[k];
		}
		return B;
	}

	template<typename T2, typename T3 = decltype(T1() + T2())>
	constexpr auto operator-(Matrix<T2, R1, C1> const &A) const {
		Matrix<T3, R1, C1> B;
		for (std::size_t k = 0; k < size(); k++) {
			B.values[k] = values[k] - A.values[k];
		}
		return B;
	}

	constexpr bool operator==(Matrix const &A) const {
		return values == A.values;
	}
	constexpr bool operator!=(Matrix const &A) const {
		return !(*this == A);
	}

	static constexpr Matrix identity() {
		static_assert(R1 == C1, "Identity called for non-square matrix");
		Matrix I;
		for (int n = 0; n < R1; n++) {
			for (int m = 0; m < C1; m++) {
				I(n, m) = T1(n == m);
			}
		}
		return I;
	}

	static constexpr Matrix zero() {
		Matrix Z;
		std::fill(Z.begin(), Z.end(), T1(0));
		return Z;
	}

	template<typename T2, typename T3 = decltype(T1() * T2()), std::enable_if_t<!IsMatrix<T2>::value, int> = 0>
	friend Matrix<T3, R1, C1> operator*(T2 const &a, Matrix const &B) {
		return B * a;
	}
	constexpr auto begin() {
		return values.begin();
	}
	constexpr auto end() {
		return values.end();
	}
	constexpr auto begin() const {
		return values.begin();
	}
	constexpr auto end() const {
		return values.end();
	}
	template<typename T2>
	operator Matrix<T2, R1, C1>() const {
		Matrix<T1, R1, C1> const &A = *this;
		Matrix<T2, R1, C1> B;
		for (int n = 0; n < R1; n++) {
			for (int m = 0; m < C1; m++) {
				B(n, m) = T2(A(n, m));
			}
		}
		return B;
	}
	template<typename, int, int>
	friend class Matrix;
private:
	std::array<T1, size()> values;
};

template<typename T1, int N>
struct ColumnVector: public Matrix<T1, N, 1> {
	using base_type = Matrix<T1, N, 1>;
	ColumnVector() = default;
	ColumnVector(ColumnVector const &other) = default;
	ColumnVector(ColumnVector &&other) = default;
	ColumnVector(T1 const &other) :
			base_type(other) {
	}
	ColumnVector(std::initializer_list<T1> const &list) :
			base_type(list) {
	}
	template<typename T2>
	ColumnVector(Matrix<T2, N, 1> const &other) :
			base_type(other) {
	}
	template<typename T2>
	ColumnVector(Matrix<T2, N, 1> &&other) :
			base_type(base_type(other)) {
	}
	ColumnVector& operator=(ColumnVector const &other) = default;
	ColumnVector& operator=(ColumnVector &&other) = default;
	ColumnVector& operator=(base_type const &other) {
		base_type::operator=(other);
		return *this;
	}
	ColumnVector& operator=(base_type &&other) {
		base_type::operator=(std::move(other));
		return *this;
	}
	ColumnVector& operator=(std::initializer_list<T1> const &list) {
		*this = ColumnVector(list);
		return *this;
	}
	T1 operator[](int i) const {
		return base_type::operator()(i, 0);
	}
	T1& operator[](int i) {
		return base_type::operator()(i, 0);
	}
	template<typename T2>
	operator ColumnVector<T2, N>() const {
		ColumnVector<T1, N> const &A = *this;
		ColumnVector<T2, N> B;
		for (int n = 0; n < N; n++) {
			B[n] = A[n];
		}
		return B;
	}
};

template<typename T1, int N>
struct RowVector: public Matrix<T1, N, 1> {
	using base_type = Matrix<T1, 1, N>;
	RowVector() = default;
	RowVector(RowVector const &other) = default;
	RowVector(base_type const &other) :
			base_type(other) {
	}
	RowVector& operator=(RowVector const &other) = default;
	RowVector& operator=(base_type const &other) {
		base_type::operator=(other);
		return *this;
	}
	T1 operator[](int i) const {
		return base_type::operator()(0, i);
	}
	T1& operator[](int i) {
		return base_type::operator()(0, i);
	}
};

template<typename T1, typename T2, int N, int M, int L, typename T3 = decltype(T1() * T2())>
Matrix<T3, N, L> operator*(Matrix<T1, N, M> const &A, Matrix<T2, M, L> const &B) {
	Matrix<T3, N, L> C;
	for (int n = 0; n < N; n++) {
		for (int l = 0; l < L; l++) {
			C(n, l) = A(n, 0) * B(0, l);
			for (int m = 1; m < M; m++) {
				C(n, l) += A(n, m) * B(m, l);
			}
		}
	}
	return C;
}

template<typename Type, int N>
Type matrixDiagonalDominance(SquareMatrix<Type, N> const &A) {
	Type pSum = Type(0);
	Type nSum = Type(0);
	for (int n = 0; n < N; n++) {
		pSum += std::abs(A(n, n));
		for (int m = 0; m < N; m++) {
			nSum += std::abs(A(n, m));
		}
	}
	nSum -= pSum;
	return (pSum - nSum) / (pSum + nSum);
}

template<typename Type, int R, int C>
auto matrixRow(Matrix<Type, R, C> const &A, int r) {
	Matrix<Type, 1, C> row;
	for (int c = 0; c < C; c++) {
		row(0, c) = A(r, c);
	}
	return row;
}

template<typename Type, int R, int C>
Matrix<Type, R, 1> matrixColumn(Matrix<Type, R, C> const &A, int c) {
	Matrix<Type, R, 1> col;
	for (int r = 0; r < R; r++) {
		col(r, 0) = A(r, c);
	}
	return col;
}

template<typename T, int N>
T matrixInverseAndDeterminant(SquareMatrix<T, N> &A) {
	T constexpr zero(0), one(1);
	T det = one;
	auto D = SquareMatrix<T, N>::identity();
	for (int i = 0; i < N; ++i) {
		T pivot = A(i, i);
		if (pivot == zero) {
			int k = i;
			do {
				k++;
				if (k >= N)
					std::abort();
				pivot = A(k, i);
			} while (pivot == zero);
			bool singular = true;
			for (int j = 0; j < N; ++j) {
				if (A(i, j) != zero)
					singular = false;
				std::swap(A(i, j), A(k, j));
				std::swap(D(i, j), D(k, j));
			}
			det = -det;
			if (singular)
				return zero;
		}
		T invP = one / pivot;
		for (int j = 0; j < N; ++j) {
			A(i, j) *= invP;
			D(i, j) *= invP;
		}
		det *= pivot;
		for (int r = 0; r < N; ++r) {
			if (r == i)
				continue;
			T f = A(r, i);
			for (int j = 0; j < N; ++j) {
				A(r, j) -= f * A(i, j);
				D(r, j) -= f * D(i, j);
			}
		}
	}
	A = D;
	return det;
}

template<typename Type, int Ndim>
Type matrixTrace(SquareMatrix<Type, Ndim> const &A) {
	Type sum = A(0, 0);
	for (int i = 1; i < Ndim; i++) {
		sum += A(i, i);
	}
	return sum;
}

template<typename Type, int R, int C>
constexpr Matrix<Type, C, R> matrixTranspose(Matrix<Type, R, C> const &B) {
	Matrix<Type, C, R> A;
	for (int r = 0; r < R; r++) {
		for (int c = 0; c < C; c++) {
			A(c, r) = B(r, c);
		}
	}
	return A;
}

template<typename Type, int Ndim>
constexpr SquareMatrix<Type, Ndim - 1> subMatrix(SquareMatrix<Type, Ndim> const &A, int row, int col) {
	SquareMatrix<Type, Ndim - 1> M;
	for (int r = 0; r < row; r++) {
		for (int c = 0; c < col; c++) {
			M(r, c) = A(r, c);
		}
		for (int c = col; c < Ndim - 1; c++) {
			M(r, c) = A(r, c + 1);
		}
	}
	for (int r = row; r < Ndim - 1; r++) {
		for (int c = 0; c < col; c++) {
			M(r, c) = A(r + 1, c);
		}
		for (int c = col; c < Ndim - 1; c++) {
			M(r, c) = A(r + 1, c + 1);
		}
	}
	return M;
}

template<typename T, int N>
constexpr T matrixDeterminant(SquareMatrix<T, N> const &A);

template<typename T, int N>
constexpr T matrixCofactor(SquareMatrix<T, N> const &A, int r, int c) {
	if constexpr (N > 1) {
		if ((r + c) % 2 == 0) {
			return matrixDeterminant(subMatrix(A, r, c));
		} else {
			return -matrixDeterminant(subMatrix(A, r, c));
		}
	} else {
		return A(0, 0);
	}
}

template<typename T, int N>
constexpr SquareMatrix<T, N> matrixCofactor(SquareMatrix<T, N> const &A) {
	SquareMatrix<T, N> C;
	for (int r = 0; r < N; r++) {
		for (int c = 0; c < N; c++) {
			C(r, c) = matrixCofactor(A, r, c);
		}
	}
	return C;
}

template<typename T, int N>
constexpr T matrixDeterminant(SquareMatrix<T, N> const &A) {
	if constexpr (N > 1) {
		T sum = A(0, 0) * matrixCofactor(A, 0, 0);
		for (int c = 1; c < N; c++) {
			auto const ac = A(0, c);
			if (ac != T(0)) {
				T cof = matrixCofactor(A, 0, c);
				T const tmp = A(0, c) * cof;
				sum += tmp;
			}
		}
		return sum;
	} else {
		return A(0, 0);
	}
}

template<typename T, int N>
constexpr SquareMatrix<T, N> matrixAdjoint(SquareMatrix<T, N> const &A) {
	return matrixTranspose(matrixCofactor(A));
}

template<typename T, int N>
constexpr auto matrixInverse(SquareMatrix<T, N> const &A) {
	T const det = matrixDeterminant(A);
	auto const adj = matrixAdjoint(A);
	return adj / det;
}

template<typename Type, int P, int Q, int M, int N>
Matrix<Type, P * M, Q * N> matrixKroneckerProduct(Matrix<Type, P, Q> const &A, Matrix<Type, M, N> const &B) {
	Matrix<Type, P * M, Q * N> C;
	for (int p = 0; p < P; p++) {
		int pQ = Q * p;
		for (int q = 0; q < Q; q++) {
			for (int m = 0; m < M; m++) {
				int mN = N * m;
				for (int n = 0; n < N; n++) {
					C(pQ + q, mN + n) = A(p, q) * B(n, m);
				}
			}
		}
	}
	return C;
}

template<typename T, int N>
void matrixQRDecomposition(SquareMatrix<T, N> const &A, SquareMatrix<T, N> &Q, SquareMatrix<T, N> &R) {
	T one = T(1);
	R = A;
	Q = SquareMatrix<T, N>::identity();
	for (int j = 0; j < N; j++) {
		for (int i = j + 1; i < N; i++) {
			SquareMatrix<T, N> G = SquareMatrix<T, N>::identity();
			T x = R(j, j), y = R(i, j);
			T inv = one / std::sqrt(x * x + y * y);
			T c = x * inv, s = -y * inv;
			G(j, j) = G(i, i) = c;
			G(i, j) = s;
			G(j, i) = -s;
			R = G * R;
			Q *= matrixTranspose(G);
		}
	}
	Q = matrixTranspose(Q);
}

template<typename T, int N>
std::array<T, N> matrixEigenvalues(SquareMatrix<T, N> const &A, int maxIters = 1000, T tol = 1e-12) {
	SquareMatrix<T, N> Ak = A;
	SquareMatrix<T, N> Q, R;
	for (int iter = 0; iter < maxIters; ++iter) {
		matrixQRDecomposition(Ak, Q, R);
		Ak = R * Q;
		std::cout << toString(Ak) << "\n";
		bool converged = true;
		for (int i = 1; i < N && converged; ++i) {
			for (int j = 0; j < i; ++j) {
				if (std::abs(Ak(i, j)) > tol) {
					printf("%i %i %e\n", i, j, std::abs(Ak(i, j)));
					converged = false;
					break;
				}
			}
		}
		if (converged) {
			break;
		}
	}
	std::array<T, N> eigenvalues;
	for (int i = 0; i < N; ++i) {
		eigenvalues[i] = Ak(i, i);
	}
	return eigenvalues;
}

template<typename T, int N>
constexpr auto matrixLUDecompose(SquareMatrix<T, N> A) {
	for (int n = 0; n < N - 1; n++) {
		T inv = T(1) / A(n, n);
		for (int k = n + 1; k < N; k++) {
			A(k, n) *= inv;
		}
		for (int j = n + 1; j < N; j++) {
			for (int k = n + 1; k < N; k++) {
				A(k, j) -= A(k, n) * A(n, j);
			}
		}
	}
	return A;
}

template<typename T, int N, int M>
void matrixLUMultiply(SquareMatrix<T, N> const &LU, Matrix<T, N, M> &X) {
	for (int k = 0; k < N; k++) {
		for (int m = 0; m < M; m++) {
			X(k, m) *= LU(k, k);
		}
		for (int j = k + 1; j < N; j++) {
			for (int m = 0; m < M; m++) {
				X(k, m) += LU(k, j) * X(j, m);
			}
		}
	}
	for (int k = N - 1; k > 0; k--) {
		for (int j = 0; j < k; j++) {
			for (int m = 0; m < M; m++) {
				X(k, m) += LU(k, j) * X(j, m);
			}
		}
	}
}

template<typename T, int N>
void matrixLURecompose(SquareMatrix<T, N> &A) {
	for (int n = N - 2; n >= 0; n--) {
		for (int j = n + 1; j < N; j++) {
			for (int k = n + 1; k < N; k++) {
				A(k, j) += A(k, n) * A(n, j);
			}
		}
		T a = A(n, n);
		for (int k = n + 1; k < N; k++) {
			A(k, n) *= a;
		}
	}
}

template<typename T, int N>
struct DiagonalMatrix {
	constexpr T operator()(int c, int r) const {
		if (c == r) {
			return D[r];
		} else {
			return T(0);
		}
	}
	constexpr T& operator()(int c, int r) {
		assert(c == r);
		if (c != r) {
			throw std::runtime_error("Attempt to assign to diagnal matrix off diagonal\n");
		}
		return D[r];
	}
private:
	std::array<T, N> D;
};

template<typename T, int N>
constexpr T matrixDeterminant(DiagonalMatrix<T, N> const &A) {
	T det = T(1);
	for (int r = 0; r < N; r++) {
		det *= A(r, r);
	}
	return det;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
template<typename T, int N>
constexpr DiagonalMatrix<T, N> matrixInverse(DiagonalMatrix<T, N> const &A) {
	DiagonalMatrix<T, N> Ainv;
	for (int r = 0; r < N; r++) {
		Ainv(r, r) = T(1) / A(r, r);
	}
	return Ainv;
}
#pragma GCC diagnostic pop

template<typename Type, int R, int C>
std::string toString(Matrix<Type, R, C> const &M) {
	auto formatEntry = [](Type const &value) {
		std::ostringstream oss;
		oss << std::scientific << std::setprecision(3) << value;
		return oss.str();
	};
	int maxLen = 0;
	for (int i = 0; i < R; ++i) {
		for (int j = 0; j < C; ++j) {
			std::string const entry = formatEntry(M(i, j));
			maxLen = std::max(maxLen, static_cast<int>(entry.size()));
		}
	}
	std::string line((C * (maxLen + 3) + 1), '-');
	line += "\n";
	std::string out = line;
	for (int i = 0; i < R; ++i) {
		for (int j = 0; j < C; ++j) {
			std::string cellStr = formatEntry(M(i, j));
			// pad on the left to align
			while (static_cast<int>(cellStr.size()) < maxLen) {
				cellStr = " " + cellStr;
			}
			out += "| " + cellStr + " ";
		}
		out += "|\n" + line;
	}

	return out;
}

template<typename Type, int R, int C>
std::string toMathematica(Matrix<Type, R, C> const &M) {
	using std::to_string;
	std::string out = "A =: {\n";
	for (int i = 0; i < R; i++) {
		out += "\t{";
		for (int j = 0; j < C; j++) {
			out += to_string(M(i, j));
			if (j + 1 < C)
				out += ",";
		}
		out += "}";
		if (i + 1 < R)
			out += ",";
		out += "\n";
	}
	out += "}";
	return out;
}

template<typename T1, typename T2, int N>
auto vectorDotProduct(Matrix<T1, N, 1> const &v1, Matrix<T2, N, 1> const &v2) {
	return (matrixTranspose(v1) * v2)(0, 0);
}

template<typename T1, typename T2, int N, typename T3 = decltype(T1() * T2())>
auto vectorDyadicProduct(Matrix<T1, N, 1> const &v1, Matrix<T2, N, 1> const &v2) {
	return v1 * matrixTranspose(v2);
}

template<typename T1, int N>
auto vectorMagnitude(Matrix<T1, N, 1> const &v) {
	return sqrt(vectorDotProduct(v, v));
}

//template<typename T, auto N, int M, typename Container>
//inline constexpr auto operator*(Matrix<T, N, M> const &A, Container const &B) {
//	Container C;
//	for (int n = 0; n < N; n++) {
//		C[n] = A(n, 0) * B[0];
//		for (int m = 1; m < M; m++) {
//			C[n] += A(n, m) * B[m];
//		}
//	}
//	return C;
//}
//
//template<typename T, auto N, typename Container>
//inline constexpr auto operator*(DiagonalMatrix<T, N> const &A, Container const &B) {
//	Container C;
//	for (int row = 0; row < N; row++) {
//		int const &col = row;
//		C[row] = A(row, col) * B[col];
//	}
//	return C;
//}

