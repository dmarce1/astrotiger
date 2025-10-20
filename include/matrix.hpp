#pragma once

#include "math.hpp"
#include "vector.hpp"

#include <iomanip>
#include <stdexcept>
#include <string>

#include "matrix_fwd.hpp"

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

template<typename TypeA, typename TypeB, int count, SymmetryType symmetryA, SymmetryType symmetryB>
constexpr auto symmetrize(Matrix<TypeA, count, count, symmetryA> const &a, Matrix<TypeB, count, count, symmetryB> const &b);

template<typename Type, int count, SymmetryType symmetry>
inline Type trace(SquareMatrix<Type, count, symmetry> const &a);

template<typename Type, int rowCount, int columnCount, SymmetryType symmetry>
inline auto transpose(Matrix<Type, rowCount, columnCount, symmetry> const &a);

template<typename Type, int count, SymmetryType symmetry>
inline Type determinant(SquareMatrix<Type, count, symmetry> const &a);

template<typename Type, int count, SymmetryType symmetry>
constexpr auto comatrix(SquareMatrix<Type, count, symmetry> const &a);

template<typename Type, int count, SymmetryType symmetry>
constexpr auto adjoint(SquareMatrix<Type, count, symmetry> const &a);

template<typename Type, int count, SymmetryType symmetry>
constexpr auto inverse(SquareMatrix<Type, count, symmetry> const &a);

template<typename, int>
constexpr auto identity();

template<typename Type, int rowCount, int columnCount>
std::string toString(Matrix<Type, rowCount, columnCount> const &M);

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
	constexpr auto getColumn(int m) const {
		Vector<Type, rowCount> c;
		for (int n = 0; n < rowCount; n++) {
			c[n] = this->operator()(n, m);
		}
		return c;
	}
	constexpr auto getRow(int n) const {
		Vector<Type, columnCount> c;
		for (int m = 0; m < columnCount; m++) {
			c[m] = this->operator()(n, m);
		}
		return c;
	}
	constexpr Matrix& setColumn(int m, Vector<Type, rowCount> const &c) {
		for (int n = 0; n < rowCount; n++) {
			this->operator()(n, m) = c[n];
		}
		return *this;
	}
	constexpr Matrix& setRow(int n, Vector<Type, columnCount> const &r) {
		for (int m = 0; m < columnCount; m++) {
			this->operator()(n, m) = r[n];
		}
		return *this;
	}
	constexpr Matrix& swapColumns(int nA, int nB) {
		using std::swap;
		auto &a = *this;
		for (int n = 0; n < columnCount; n++) {
			swap(a(n, nB), a(n, nA));
		}
		return *this;
	}
	constexpr Matrix& swapRows(int nA, int nB) {
		using std::swap;
		auto &a = *this;
		for (int m = 0; m < rowCount; m++) {
			swap(a(nA, m), a(nB, m));
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
	constexpr auto submatrix(int p, int q) const {
		auto const &a = *this;
		Matrix<Type, rowCount - 1, columnCount - 1> sub;
		for (int n = 0, j = 0; j + 1 < rowCount; j++, n++) {
			if (n == p) {
				n++;
			}
			for (int m = 0, k = 0; k + 1 < columnCount; k++, m++) {
				if (m == q) {
					m++;
				}
				sub(j, k) = a(n, m);
			}
		}
		return sub;
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
	return determinant(submatrix(row, col));
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
inline auto transpose(Matrix<Type, rowCount, columnCount, symmetry> const &a) {
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
inline Type determinant(SquareMatrix<Type, count, symmetry> const &a) {
	if constexpr (count > 1) {
		Type det(0);
		for (int n = 0; n < count; n++) {
			det += a(0, n) * a.cofactor(0, n);
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

template<typename Type, int count>
constexpr auto identity() {
	SquareMatrix<Type, count, SymmetryType::symmetric> I;
	for (int n = 0; n < count; n++) {
		I(n, n) = Type(1);
		for (int m = 0; m < n; m++) {
			I(n, m) = I(m, n) = Type(0);
		}
	}
	return I;
}

template<typename Type, int rowCount, int columnCount>
std::string toString(Matrix<Type, rowCount, columnCount> const &M) {
	using std::max;
	auto formatEntry = [](Type const &value) {
		std::ostringstream oss;
		oss << std::scientific << std::setprecision(3) << value;
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
	return out;
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
		C[l] = A(l, 0) * B[0];
		for (int m = 1; m < rowCount; m++) {
			C[l] += A(l, m) * B[m];
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
constexpr auto symmetric(SquareMatrix<Type, count, symmetry> const &a) {
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
constexpr auto antisymmetric(SquareMatrix<Type, count, symmetry> const &a) {
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

template<typename Type, int count>
SquareMatrix<Type, count, SymmetryType::symmetric> diagonal(Vector<Type, count> &λ) {
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

#undef ASSIGN_BINARY_OP
