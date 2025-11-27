#pragma once

#include <bitset>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stack>

#include "matrix_fwd.hpp"
#include "math.hpp"
#include "vector.hpp"

#define ASSIGN_BINARY_OP(op)                                            \
   template<typename OtherType>                                         \
   constexpr SignedReference& operator op##= (OtherType const &value) { \
      ref_ op##= signedValue(value);                                    \
      return *this;                                                     \
   }

struct MatrixException: public std::runtime_error {
	explicit MatrixException(std::string const &msg) :
			std::runtime_error(msg) {
	}
};

template<typename TypeA, typename TypeB, int rowCount, int commonCount, int columnCount, SymmetryType symmetryA, SymmetryType symmetryB>
constexpr auto operator*(Matrix<TypeA, rowCount, commonCount, symmetryA> const &a, Matrix<TypeB, commonCount, columnCount, symmetryB> const &b);

template<typename TypeA, int count, SymmetryType symmetryA>
constexpr auto symmetrize(Matrix<TypeA, count, count, symmetryA> const &a);

template<typename Type, int count, SymmetryType symmetry>
inline Type trace(SquareMatrix<Type, count, symmetry> const &a);

template<typename Type, int rowCount, int columnCount, SymmetryType symmetry>
inline auto constexpr transpose(Matrix<Type, rowCount, columnCount, symmetry> const &a);

template<typename Type, int count, SymmetryType symmetry>
constexpr Type determinant(SquareMatrix<Type, count, symmetry> const &a);

template<typename Type, int count, SymmetryType symmetry>
constexpr auto comatrix(SquareMatrix<Type, count, symmetry> const &a);

template<typename Type, int count, SymmetryType symmetry>
constexpr auto adjoint(SquareMatrix<Type, count, symmetry> const &a);

template<typename Type, int count, SymmetryType symmetry>
constexpr auto inverse(SquareMatrix<Type, count, symmetry> const &a);

template<typename, int, SymmetryType symmetry = SymmetryType::symmetric>
constexpr auto identity();

template<typename TypeA, typename TypeB, int rowCount, int columnCount>
constexpr auto operator*(Vector<TypeA, rowCount> const&, Vector<TypeB, columnCount> const&);

template<typename Type>
struct SignedReference {
	constexpr SignedReference() = default;
	constexpr SignedReference(SignedReference const&) = default;
	constexpr SignedReference(SignedReference&&) = default;
	constexpr SignedReference(int sign, Type &ref) :
			sign_(sign), ref_(ref) {
	}
	template<typename OtherType>
	constexpr SignedReference& operator=(OtherType const &value) {
		ref_ = signedValue(value);
		return *this;
	}
	template<typename OtherType>
	constexpr SignedReference& operator=(OtherType &&value) {
		ref_ = signedValue(std::move(value));
		return *this;
	}
	constexpr operator Type() const {
		return signedValue(ref_);
	}
	ASSIGN_BINARY_OP(+)
	ASSIGN_BINARY_OP(-)
	ASSIGN_BINARY_OP(*)
	ASSIGN_BINARY_OP(/)
private:
	constexpr Type signedValue(Type const &value) const {
		constexpr Type zero(0);
		return Type((sign_ == 0) ? zero : ((sign_ > 0) ? value : -value));
	}
	int const sign_;
	Type &ref_;
};

template<typename Type, int rowCount, int columnCount, SymmetryType symmetry>
struct Matrix {
	using reference = SignedReference<Type>;
	constexpr Matrix() = default;
	constexpr Matrix(Matrix&&) = default;
	constexpr Matrix(Matrix const&) = default;
	template<typename OtherType>
	constexpr Matrix(Matrix<OtherType, rowCount, columnCount, symmetry> const &other) :
			μ_(other.μ_) {
	}
	constexpr Matrix(std::array<Vector<Type, rowCount>, columnCount> const &init) {
		for (int n = 0; n < rowCount; n++) {
			for (int m = 0; m < columnCount; m++) {
				this->operator()(n, m) = init[m][n];
			}
		}
	}
	constexpr Matrix(std::initializer_list<Type> init) :
			μ_(init) {
	}
	constexpr Matrix(Type const &init) :
			μ_(init) {
	}
	constexpr Matrix& operator=(Matrix const&) = default;
	constexpr Matrix& operator=(Matrix&&) = default;
	template<typename OtherType>
	constexpr Matrix& operator=(Matrix<OtherType, rowCount, columnCount> const &other) {
		μ_ = other.μ_;
		return *this;
	}
	constexpr Type operator()(int n, int m) const {
		auto const sgn = sign(n, m);
		auto const value = μ_[flatIndex(n, m)];
		return (sgn > 0) ? value : ((sgn < 0) ? -value : Type(0));
	}
	constexpr reference operator()(int n, int m) {
		return reference(sign(n, m), μ_[flatIndex(n, m)]);
	}
	constexpr auto getRow(int n) const {
		Vector<Type, columnCount> row;
		for (int m = 0; m < columnCount; m++) {
			row[m] = this->operator()(n, m);
		}
		return row;
	}
	constexpr Matrix& setRow(int n, Vector<Type, columnCount> const &row) {
		for (int m = 0; m < columnCount; m++) {
			this->operator()(n, m) = row[m];
		}
		return *this;
	}
	constexpr auto getColumn(int m) const {
		Vector<Type, rowCount> col;
		for (int n = 0; n < rowCount; n++) {
			col[n] = this->operator()(n, m);
		}
		return col;
	}
	constexpr Matrix& setColumn(int m, Vector<Type, rowCount> const &col) {
		for (int n = 0; n < rowCount; n++) {
			this->operator()(n, m) = col[n];
		}
		return *this;
	}
	constexpr Matrix& swapColumns(int nA, int nB) {
		if (nA == nB) {
			return *this;
		}
		auto &a = *this;
		for (int n = 0; n < rowCount; n++) {
			Type const tmp = a(n, nA);
			a(n, nA) = a(n, nB);
			a(n, nB) = tmp;
		}
		return *this;
	}
	constexpr Matrix& swapRows(int nA, int nB) {
		if (nA == nB) {
			return *this;
		}
		auto &a = *this;
		for (int m = 0; m < columnCount; m++) {
			Type const tmp = a(nA, m);
			a(nA, m) = a(nB, m);
			a(nB, m) = tmp;
		}
		return *this;
	}
	template<SymmetryType otherSymmetry, std::enable_if_t<(symmetry == SymmetryType::asymmetric) || (symmetry == otherSymmetry), int> = 0>
	friend constexpr Matrix& operator+=(Matrix<Type, rowCount, columnCount, symmetry> &A, Matrix<Type, rowCount, columnCount, otherSymmetry> const &B) {
		A = A + B;
		return A;
	}
	template<SymmetryType otherSymmetry, std::enable_if_t<(symmetry == SymmetryType::asymmetric) || (symmetry == otherSymmetry), int> = 0>
	friend constexpr Matrix& operator-=(Matrix<Type, rowCount, columnCount, symmetry> &A, Matrix<Type, rowCount, columnCount, otherSymmetry> const &B) {
		A = A - B;
		return A;
	}
	constexpr Matrix& operator*=(Type const &a) {
		μ_ *= a;
		return *this;
	}
	constexpr Matrix& operator/=(Type const &a) {
		μ_ /= a;
		return *this;
	}
	constexpr bool operator==(Matrix const &A) const {
		return μ_ == A.μ_;
	}
	constexpr bool operator!=(Matrix const &A) const {
		return !(*this == A);
	}
	template<typename OtherType, typename = std::enable_if_t<!IsMatrix<OtherType>::value && !IsVector<OtherType>::value>>
	friend constexpr auto operator*(Matrix const &b, OtherType const &c) {
		using ReturnType = decltype(Type() * OtherType());
		Matrix<ReturnType, rowCount, columnCount, symmetry> a;
		a.μ_ = b.μ_ * c;
		return a;
	}
	template<typename OtherType, typename = std::enable_if_t<!IsMatrix<OtherType>::value && !IsVector<OtherType>::value>>
	friend constexpr auto operator*(OtherType const &c, Matrix const &b) {
		return b * c;
	}
	template<typename OtherType, typename = std::enable_if_t<!IsMatrix<OtherType>::value>>
	friend constexpr auto operator/(Matrix const &b, OtherType const &c) {
		using ReturnType = decltype(Type() / OtherType());
		Matrix<ReturnType, rowCount, columnCount> a;
		a.μ_ = b.μ_ / c;
		return a;
	}
	friend constexpr Matrix operator+(Matrix const &a) {
		return a;
	}
	friend constexpr Matrix operator-(Matrix a) {
		a.μ_ = -a.μ_;
		return a;
	}
	template<SymmetryType otherSymmetry>
	friend constexpr auto operator+(Matrix const &a, Matrix<Type, rowCount, columnCount, otherSymmetry> const &b) {
		if constexpr ((symmetry == SymmetryType::symmetric) && (otherSymmetry == SymmetryType::symmetric)) {
			Matrix<Type, rowCount, columnCount, SymmetryType::symmetric> c;
			for (int l = 0; l < rowCount; l++) {
				for (int m = 0; m <= l; m++) {
					c(l, m) = a(l, m) + b(l, m);
				}
			}
			return c;
		} else if constexpr ((symmetry == SymmetryType::antisymmetric) && (otherSymmetry == SymmetryType::antisymmetric)) {
			Matrix<Type, rowCount, columnCount, SymmetryType::antisymmetric> c;
			for (int l = 0; l < rowCount; l++) {
				for (int m = 0; m < l; m++) {
					c(l, m) = a(l, m) + b(l, m);
				}
			}
			return c;
		} else {
			Matrix<Type, rowCount, columnCount> c;
			for (int l = 0; l < rowCount; l++) {
				for (int m = 0; m < columnCount; m++) {
					c(l, m) = a(l, m) + b(l, m);
				}
			}
			return c;
		}
	}
	template<SymmetryType otherSymmetry>
	friend constexpr auto operator-(Matrix const &a, Matrix<Type, rowCount, columnCount, otherSymmetry> const &b) {
		if constexpr ((symmetry == SymmetryType::symmetric) && (otherSymmetry == SymmetryType::symmetric)) {
			Matrix<Type, rowCount, columnCount, SymmetryType::symmetric> c;
			for (int l = 0; l < rowCount; l++) {
				for (int m = 0; m <= l; m++) {
					c(l, m) = a(l, m) - b(l, m);
				}
			}
			return c;
		} else if constexpr ((symmetry == SymmetryType::antisymmetric) && (otherSymmetry == SymmetryType::antisymmetric)) {
			Matrix<Type, rowCount, columnCount, SymmetryType::antisymmetric> c;
			for (int l = 0; l < rowCount; l++) {
				for (int m = 0; m < l; m++) {
					c(l, m) = a(l, m) - b(l, m);
				}
			}
			return c;
		} else {
			Matrix<Type, rowCount, columnCount> c;
			for (int l = 0; l < rowCount; l++) {
				for (int m = 0; m < columnCount; m++) {
					c(l, m) = a(l, m) - b(l, m);
				}
			}
			return c;
		}
	}
	constexpr auto getSubMatrix(int p, int q) const {
		auto const &a = *this;
		Matrix<Type, rowCount - 1, columnCount - 1> sub;
		for (int n = 0; n + 1 < rowCount; n++) {
			int const j = int(n >= p);
			for (int m = 0; m + 1 < columnCount; m++) {
				int const k = int(m >= q);
				sub(n, m) = a(n + j, m + k);
			}
		}
		return sub;
	}
	constexpr auto setSubmatrix(int p, int q, Matrix<Type, rowCount - 1, columnCount - 1> const &sub) {
		auto &a = *this;
		for (int n = 0, j = 0; j + 1 < rowCount; j++, n++) {
			if (n == p) {
				n++;
			}
			for (int m = 0, k = 0; k + 1 < columnCount; k++, m++) {
				if (m == q) {
					m++;
				}
				a(n, m) = sub(j, k);
			}
		}
		return *this;
	}
	constexpr Type minor(int, int) const;
	constexpr Type cofactor(int, int) const;
	static constexpr size_t size() {
		if constexpr (symmetry == SymmetryType::symmetric) {
			static_assert(rowCount == columnCount);
			constexpr auto count = rowCount;
			return ((count + 1) * count) >> 1;
		} else if constexpr (symmetry == SymmetryType::antisymmetric) {
			static_assert(rowCount == columnCount);
			constexpr auto count = rowCount;
			return ((count - 1) * count) >> 1;
		} else {
			return columnCount * rowCount;
		}
	}
	constexpr bool isNearZero(int n, int m) const {
		using std::abs;
		using std::sqrt;
		auto const z = sqrt(std::numeric_limits<Type>::epsilon()) * frobeniusNorm(*this);
		return bool(abs((*this)(n, m)) < z);
	}
	template<typename, int, int, SymmetryType>
	friend class Matrix;
private:
	constexpr int flatIndex(int n, int m) const {
		using std::min;
		using std::max;
		if constexpr (symmetry == SymmetryType::symmetric) {
			if (n < m) {
				return flatIndex(m, n);
			} else {
				return ((n * (n + 1)) >> 1) + m;
			}
		} else if constexpr (symmetry == SymmetryType::antisymmetric) {
			if (n < m) {
				return flatIndex(m, n);
			} else {
				return ((n * (n - 1)) >> 1) + m;
			}
		} else {
			return n * columnCount + m;
		}
	}
	constexpr int sign(int n, int m) const {
		using std::min;
		using std::max;
		if constexpr (symmetry == SymmetryType::antisymmetric) {
			if (n > m) {
				return +1;
			} else if (n < m) {
				return -1;
			} else {
				return 0;
			}
		} else {
			return +1;
		}
	}
	Vector<Type, size()> μ_;
}
;

constexpr SymmetryType multSymmetry(SymmetryType a, SymmetryType b) {
	if (a == SymmetryType::symmetric) {
		if (b == SymmetryType::symmetric) {
			return SymmetryType::symmetric;
		} else if (b == SymmetryType::antisymmetric) {
			return SymmetryType::antisymmetric;
		}
	} else if (a == SymmetryType::antisymmetric) {
		if (b == SymmetryType::symmetric) {
			return SymmetryType::antisymmetric;
		} else if (b == SymmetryType::antisymmetric) {
			return SymmetryType::symmetric;
		}
	}
	return SymmetryType::asymmetric;
}

template<typename TypeA, typename TypeB, int rowCount, int commonCount, int columnCount, SymmetryType symmetryA, SymmetryType symmetryB>
constexpr auto operator*(Matrix<TypeA, rowCount, commonCount, symmetryA> const &a, Matrix<TypeB, commonCount, columnCount, symmetryB> const &b) {
	using ReturnType = decltype(TypeA() * TypeB());
	Matrix<ReturnType, rowCount, columnCount> c;
	for (int l = 0; l < rowCount; l++) {
		for (int m = 0; m < columnCount; m++) {
			c(l, m) = ReturnType(0);
			for (int n = 0; n < commonCount; n++) {
				c(l, m) += a(l, n) * b(n, m);
			}
		}
	}
	return c;
}

template<typename Type, int rowCount, int columnCount, SymmetryType symmetry>
constexpr Type Matrix<Type, rowCount, columnCount, symmetry>::minor(int row, int col) const {
	return determinant(getSubMatrix(row, col));
}

template<typename Type, int rowCount, int columnCount, SymmetryType symmetry>
constexpr Type Matrix<Type, rowCount, columnCount, symmetry>::cofactor(int row, int col) const {
	return Type(nonepow(row + col)) * minor(row, col);
}

template<typename Type, int count, SymmetryType symmetry>
inline Type trace(SquareMatrix<Type, count, symmetry> const &a) {
	Type traceA = Type(0);
	for (int n = 0; n < count; n++) {
		traceA += a(n, n);
	}
	return traceA;
}

template<typename Type, int rowCount, int columnCount, SymmetryType symmetry>
inline constexpr auto transpose(Matrix<Type, rowCount, columnCount, symmetry> const &a) {
	if constexpr (symmetry == SymmetryType::symmetric) {
		return a;
	} else if constexpr (symmetry == SymmetryType::antisymmetric) {
		return -a;
	} else {
		Matrix<Type, columnCount, rowCount> transposeA;
		for (int n = 0; n < columnCount; n++) {
			for (int m = 0; m < rowCount; m++) {
				transposeA(n, m) = a(m, n);
			}
		}
		return transposeA;
	}
}

template<typename Type, int count, SymmetryType symmetry>
constexpr Type determinant(SquareMatrix<Type, count, symmetry> const &a) {
	if constexpr (count > 1) {
		Type det(0);
		for (int n = 0; n < count; n++) {
			Type const az = a(0, n);
			if (az != Type(0)) {
				det += a(0, n) * a.cofactor(0, n);
			}
		}
		return det;
	} else {
		return a(0, 0);
	}
}

template<typename Type, int count, SymmetryType symmetry>
constexpr auto comatrix(SquareMatrix<Type, count, symmetry> const &a) {
	SquareMatrix<Type, count> coA;
	for (int n = 0; n < count; n++) {
		for (int m = 0; m < count; m++) {
			coA(n, m) = a.cofactor(n, m);
		}
	}
	return coA;
}

template<typename Type, int count, SymmetryType symmetry>
constexpr auto adjoint(SquareMatrix<Type, count, symmetry> const &a) {
	return transpose(comatrix(a));
}

template<typename Type, int count, SymmetryType symmetry>
constexpr auto inverse(SquareMatrix<Type, count, symmetry> const &a) {
	auto const det = determinant(a);
	if (det == Type(0)) {
		throw MatrixException("Matrix inversion failed: determinant = 0.");
	}
	return adjoint(a) / det;
}

template<typename Type, int count, SymmetryType symmetry>
constexpr auto identity() {
	static_assert(symmetry != SymmetryType::antisymmetric);
	SquareMatrix<Type, count, symmetry> I;
	for (int n = 0; n < count; n++) {
		I(n, n) = Type(1);
		for (int m = 0; m < n; m++) {
			I(n, m) = I(m, n) = Type(0);
		}
	}
	return I;
}

template<typename TypeA, typename TypeB, int rowCount, int columnCount, SymmetryType symmetry>
constexpr auto operator*(Matrix<TypeA, rowCount, columnCount, symmetry> const &A, Vector<TypeB, columnCount> const &B) {
	using ReturnType = decltype(TypeA() * TypeB());
	Vector<ReturnType, rowCount> C;
	for (int l = 0; l < rowCount; l++) {
		C[l] = A(l, 0) * B[0];
		for (int m = 1; m < columnCount; m++) {
			C[l] += A(l, m) * B[m];
		}
	}
	return C;
}

template<typename TypeA, typename TypeB, int rowCount, int columnCount, SymmetryType symmetry>
constexpr auto operator*(Vector<TypeA, rowCount> const &A, Matrix<TypeB, rowCount, columnCount, symmetry> const &B) {
	using ReturnType = decltype(TypeA() * TypeB());
	Vector<ReturnType, columnCount> C;
	for (int l = 0; l < columnCount; l++) {
		C[l] = A[0] * B(0, l);
		for (int m = 1; m < rowCount; m++) {
			C[l] += A[m] * B(m, l);
		}
	}
	return C;
}

template<typename TypeA, typename TypeB, int rowCount, int columnCount>
constexpr auto operator*(Vector<TypeA, rowCount> const &a, Vector<TypeB, columnCount> const &b) {
	using ReturnType = decltype(TypeA() * TypeB());
	Matrix<ReturnType, rowCount, columnCount> c;
	for (int n = 0; n < rowCount; n++) {
		for (int m = 0; m < columnCount; m++) {
			c(n, m) = a[n] * b[m];
		}
	}
	return c;
}

template<typename Type, int count>
constexpr auto sqr(Vector<Type, count> const &a) {
	using ReturnType = decltype(Type() * Type());
	SquareMatrix<ReturnType, count, SymmetryType::symmetric> b;
	for (int n = 0; n < count; n++) {
		for (int m = 0; m <= n; m++) {
			b(n, m) = a[n] * a[m];
		}
	}
	return b;
}

template<typename Type, int count, SymmetryType symmetry>
constexpr auto symmetrize(SquareMatrix<Type, count, symmetry> const &a) {
	constexpr auto zero = Type(0);
	constexpr auto one = Type(1);
	constexpr auto two = Type(2);
	constexpr auto half = one / two;
	SquareMatrix<Type, count, SymmetryType::symmetric> c;
	if constexpr (symmetry == SymmetryType::symmetric) {
		c = a;
	} else if constexpr (symmetry == SymmetryType::symmetric) {
		c = SquareMatrix<Type, count, SymmetryType::symmetric>(zero);
	} else {
		for (int n = 0; n < count; n++) {
			c(n, n) = a(n, n);
			for (int m = 0; m < n; m++) {
				c(n, m) = half * (a(n, m) + a(m, n));
			}
		}
	}
	return c;
}

template<typename Type, int count, SymmetryType symmetry>
constexpr auto antisymmetrize(SquareMatrix<Type, count, symmetry> const &a) {
	constexpr auto zero = Type(0);
	constexpr auto one = Type(1);
	constexpr auto two = Type(2);
	constexpr auto half = one / two;
	SquareMatrix<Type, count, SymmetryType::antisymmetric> c;
	if constexpr (symmetry == SymmetryType::symmetric) {
		c = SquareMatrix<Type, count, SymmetryType::symmetric>(zero);
	} else if constexpr (symmetry == SymmetryType::symmetric) {
		c = a;
	} else {
		for (int n = 0; n < count; n++) {
			for (int m = 0; m < n; m++) {
				c(n, m) = half * (a(n, m) - a(m, n));
			}
		}
	}
	return c;
}

template<typename Type, int N, int ncol, SymmetryType symmetry>
constexpr auto frobeniusNorm(Matrix<Type, N, ncol, symmetry> const &a) {
	using std::sqrt;
	Type n2 = Type(0);
	for (int j = 0; j < N; j++) {
		for (int k = 0; k < ncol; k++) {
			n2 += a(j, k) * a(j, k);
		}
	}
	return sqrt(n2);
}

template<typename Type, int count>
constexpr auto diagonal(Vector<Type, count> const &λ) {
	auto A = identity<Type, count>();
	for (int i = 0; i < count; i++) {
		A(i, i) = λ[i];
	}
	return A;
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

template<typename Type, int rowCount, int columnCount, SymmetryType symmetry>
std::ostream& operator<<(std::ostream &os, Matrix<Type, rowCount, columnCount, symmetry> const &M) {
	using std::max;
	auto const isZero = [M](Type n) {
		using std::abs;
		using std::sqrt;
		auto const z = sqrt(std::numeric_limits < Type > ::epsilon()) * frobeniusNorm(M);
		return bool(abs(n) < z);
	};
	auto formatEntry = [isZero](Type const &value) {
		std::ostringstream oss;
		if (isZero(value)) {
			oss << "";
		} else {
			oss << std::defaultfloat << std::setprecision(6) << value;
		}
		return oss.str();
	};
	size_t maxLen = 0;
	for (int i = 0; i < rowCount; ++i) {
		for (int j = 0; j < columnCount; ++j) {
			std::string const entry = formatEntry(M(i, j));
			maxLen = max(maxLen, entry.size());
		}
	}
	std::string line((columnCount * (maxLen + 3) + 1), '-');
	line += "\n";
	std::string out = line;
	for (int i = 0; i < rowCount; ++i) {
		for (int j = 0; j < columnCount; ++j) {
			std::string cellStr = formatEntry(M(i, j));
			while (cellStr.size() < maxLen) {
				cellStr = " " + cellStr;
			}
			out += "| " + cellStr + " ";
		}
		out += "|\n" + line;
	}
	os << out;
	return os;
}

template<typename T, auto N, auto ncol, auto P, auto Q, SymmetryType sym1, SymmetryType sym2>
constexpr auto kroneckerProduct(Matrix<T, N, ncol, sym1> const &A, Matrix<T, P, Q, sym2> const &B) {
	Matrix<T, N * P, ncol * Q> C;
	for (int n = 0; n < N; n++) {
		for (int p = 0; p < P; p++) {
			for (int m = 0; m < ncol; m++) {
				for (int q = 0; q < Q; q++) {
					C(P * n + p, Q * m + q) = A(n, m) * B(p, q);
				}
			}
		}
	}
	return C;
}

template<typename T, auto N>
struct LUDecomposition {
	SquareMatrix<T, N> L = SquareMatrix<T, N>(T(0));
	SquareMatrix<T, N> U = SquareMatrix<T, N>(T(0));
	SquareMatrix<T, N> P = identity<T, N, SymmetryType::asymmetric>();
};

/*  0 0 1 0
 *  0 1 0 0
 *  1 0 0 0
 *  0 0 0 1
 */
template<typename T, auto N, int M = 0>
auto luFactorize(SquareMatrix<T, N> const A, SquareMatrix<T, N> const P = identity<T, N, SymmetryType::asymmetric>()) {
	using std::abs;
	using std::numeric_limits;
	using std::sqrt;
	constexpr auto eps = sqrt(numeric_limits < T > ::epsilon());
	if constexpr (N > M) {
		auto const tiny = eps * frobeniusNorm(A);
		for (int p = M; p < N; p++) {
			if (abs(A(M, p)) >= tiny) {
				auto B = A;
				auto Q = P;
				B.swapRows(M, p);
				Q.swapRows(M, p);
				for (int n = M + 1; n < N; n++) {
					B(n, M) /= B(M, M);
					for (int m = M + 1; m < N; m++) {
						B(n, m) = B(n, m) - B(n, M) * B(M, m);
					}
				}
				auto const rc = luFactorize<T, N, M + 1>(B, Q);
				if (std::get<0>(rc)) {
					return rc;
				}
			}
		}
		assert(M > 0);
		return std::tuple(false, P, A);
	}
	return std::tuple(true, P, A);
}

template<typename T, auto N>
auto luDecomposition(SquareMatrix<T, N> A) {
	using std::abs;
	using std::numeric_limits;
	using std::sqrt;
	using std::swap;
	std::stack<SquareMatrix<T, N>> stack;
	std::array<int, N> pivot;
	LUDecomposition<T, N> lu;
	auto [success, P, LU] = luFactorize(A);
	lu.P = P;
	for (int j = 0; j < N; j++) {
		for (int k = 0; k <= j; k++) {
			lu.L(j, k) = LU(j, k);
		}
		for (int k = j; k < N; k++) {
			lu.U(j, k) = LU(j, k);
		}
		lu.L(j, j) = T(1);
	}
//	std::cout << transpose(lu.P);
//	std::cout << lu.L;
//	std::cout << lu.U;
//	std::cout << A;
//	std::cout << lu.L * lu.U *  transpose(lu.P);
	return lu;
}

template<typename T, auto N>
constexpr auto permutationMatrix(Vector<int, N> p) {
	SquareMatrix<T, N> P(T(0));
	for (int n = 0; n < N; n++) {
		P(n, p[n]) = T(1);
	}
	return P;
}

template<typename T, auto N, auto ncol>
auto inPlaceDecomposition(Matrix<T, N, ncol> A) {
	using std::min;
	using std::max;
	constexpr auto nmin = min(N, ncol);
	constexpr auto nmax = max(N, ncol);
	constexpr auto I = identity<T, nmax, SymmetryType::asymmetric>();
	SquareMatrix<T, nmin> LU;
//	printf( "%i %i\n", N, ncol);

	for (int n = 0; n < nmin; n++) {
		for (int k = 0; k < nmin; k++) {
			LU(n, k) = A(n, k);
		}
	}
	auto const lu = luDecomposition(LU);
	auto L = I;
	auto U = I;
	Matrix<T, N, ncol> R;
	for (int n = 0; n < nmin; n++) {
		for (int k = 0; k < nmin; k++) {
			L(n, k) = lu.L(n, k);
			U(n, k) = lu.U(n, k);
			R(n, k) = lu.P(k, n);
		}
		for (int k = N; k < ncol; k++) {
			R(n, k) = A(n, k);
		}
	}
	return std::tuple(R, L, U);
}

template<typename T, int N>
using SymmetricMatrix = SquareMatrix<T, N, SymmetryType::symmetric>;

#undef ASSIGN_BINARY_OP
