#pragma once

#include "math.hpp"
#include "vector.hpp"

#include <iomanip>
#include <stdexcept>
#include <string>

struct SingularMatrixError: public std::runtime_error {
	explicit SingularMatrixError(std::string const &msg) :
			std::runtime_error(msg) {
	}
};

template<typename Type, int rowCount, int columnCount>
struct Matrix;

template<typename Type, int count>
using SquareMatrix = Matrix<Type, count, count>;

template<typename >
struct IsMatrix {
	static constexpr bool value = false;
};

template<typename Type, int rowCount, int columnCount>
struct IsMatrix<Matrix<Type, rowCount, columnCount>> {
	static constexpr bool value = true;
};

template<typename Type, int rowCount, int columnCount>
struct Matrix {
	constexpr Matrix() = default;
	constexpr Matrix(Matrix&&) = default;
	constexpr Matrix(Matrix const&) = default;
	template<typename OtherType>
	constexpr Matrix(Matrix<OtherType, rowCount, columnCount> const &other) :
			μ_(other.μ_) {
	}
	constexpr Matrix(std::array<std::array<Type, columnCount>, rowCount> const &init) {
		auto &α = *this;
		for (int n = 0; n < rowCount; n++) {
			for (int m = 0; m < columnCount; m++) {
				α(n, m) = init[n][m];
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
	constexpr Type const& operator()(int n, int m) const {
		return μ_[n * columnCount + m];
	}
	constexpr Type& operator()(int n, int m) {
		return μ_[n * columnCount + m];
	}
	constexpr Matrix& operator+=(Matrix const &A) {
		μ_ += A.μ_;
		return *this;
	}
	constexpr Matrix& operator-=(Matrix const &A) {
		μ_ -= A.μ_;
		return *this;
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
	template<typename OtherType, typename = std::enable_if_t<!IsMatrix<OtherType>::value>>
	friend constexpr auto operator*(Matrix const &b, OtherType const &c) {
		using ReturnType = decltype(Type() * OtherType());
		Matrix<ReturnType, rowCount, columnCount> a;
		a.μ_ = b.μ_ * c;
		return a;
	}
	template<typename OtherType, typename = std::enable_if_t<!IsMatrix<OtherType>::value>>
	friend constexpr auto operator*(OtherType const &c, Matrix const &b) {
		return b * c;
	}
	template<typename OtherType, typename ReturnType = decltype(Type() / OtherType())>
	friend constexpr auto operator/(Matrix const &b, OtherType const &c) {
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
	friend constexpr Matrix operator+(Matrix a, Matrix const &b) {
		a.μ_ += b.μ_;
		return a;
	}
	friend constexpr Matrix operator-(Matrix a, Matrix const &b) {
		a.μ_ -= b.μ_;
		return a;
	}
	constexpr Matrix<Type, rowCount - 1, columnCount - 1> submatrix(int row, int column) const {
		static_assert(rowCount == columnCount);
		constexpr int count = rowCount - 1;
		auto const &a = *this;
		Matrix<Type, count, count> sub;
		for (int rowFrom, rowTo = 0; rowTo < count; rowTo++, rowFrom++) {
			if (rowFrom == row) {
				rowFrom++;
			}
			for (int columnFrom, columnTo = 0; columnTo < count; columnTo++, columnFrom++) {
				if (columnFrom == column) {
					columnFrom++;
				}
				sub(rowTo, columnTo) = a(rowFrom, columnFrom);
			}
		}
		return sub;
	}
	constexpr Type minor(int, int) const;
	constexpr Type cofactor(int, int) const;
	static constexpr size_t size() {
		return columnCount * rowCount;
	}
	template<typename, int, int>
	friend class Matrix;
private:
	Vector<Type, size()> μ_;
};

template<typename Type, int count>
inline Type determinant(SquareMatrix<Type, count> const&);

template<typename TypeA, typename TypeB, int rowCount, int commonCount, int columnCount>
constexpr auto operator*(Matrix<TypeA, rowCount, commonCount> const &a, Matrix<TypeB, commonCount, columnCount> const &b) {
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

template<typename Type, int rowCount, int columnCount>
constexpr Type Matrix<Type, rowCount, columnCount>::minor(int row, int col) const {
	return determinant(submatrix(row, col));
}

template<typename Type, int rowCount, int columnCount>
constexpr Type Matrix<Type, rowCount, columnCount>::cofactor(int row, int col) const {
	return nonepow(row + col) * minor(row, col);
}

template<typename Type, int count>
inline Type trace(SquareMatrix<Type, count> const &a) {
	Type traceA = Type(0);
	for (int n = 0; n < count; n++) {
		traceA += a(n, n);
	}
	return traceA;
}

template<typename Type, int rowCount, int columnCount>
inline auto transpose(Matrix<Type, rowCount, columnCount> const &a) {
	Matrix<Type, columnCount, rowCount> transposeA;
	for (int n = 0; n < rowCount; n++) {
		for (int m = 0; m < columnCount; m++) {
			transposeA(m, n) = a(n, m);
		}
	}
	return transposeA;
}

template<typename Type, int count>
inline Type determinant(SquareMatrix<Type, count> const &a) {
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

template<typename Type, int count>
constexpr SquareMatrix<Type, count> comatrix(SquareMatrix<Type, count> const &a) {
	SquareMatrix<Type, count> coA;
	for (int n = 0; n < count; n++) {
		for (int m = 0; m < count; m++) {
			coA(n, m) = a.cofactor(n, m);
		}
	}
	return coA;
}

template<typename Type, int count>
constexpr SquareMatrix<Type, count> adjoint(SquareMatrix<Type, count> const &a) {
	return transpose(comatrix(a));
}

template<typename Type, int count>
constexpr SquareMatrix<Type, count> inverse(SquareMatrix<Type, count> const &a) {
	auto const det = determinant(a);
	if (det == Type(0)) {
		throw SingularMatrixError("Matrix inversion failed: determinant = 0.");
	}
	return adjoint(a) / det;
}

template<typename Type, int count>
constexpr SquareMatrix<Type, count> identity() {
	SquareMatrix<Type, count> I;
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

template<typename Type1, typename Type2, int rowCount, int columnCount, typename ReturnType = decltype(Type1() * Type2())>
constexpr Matrix<ReturnType, rowCount, columnCount> operator*(Vector<Type1, rowCount> const &a, Vector<Type2, columnCount> const &b) {
	Matrix<ReturnType, rowCount, columnCount> outer;
	for (int n = 0; n < rowCount; n++) {
		for (int m = 0; m < columnCount; m++) {
			outer(n, m) = a[n] * b[m];
		}
	}
	return outer;
}

