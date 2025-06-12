#pragma once

#include <iomanip>
#include <ios>
#include <string>

template<typename Type, int Size>
struct SquareMatrix;


template<typename Type, int Size>
struct SquareMatrix {
	static constexpr int size() {
		return Size * Size;
	}
	static constexpr int rowCount() {
		return Size;
	}
	static constexpr int columnCount() {
		return Size;
	}
	constexpr SquareMatrix() {
	}
	constexpr SquareMatrix(std::array<std::array<Type, Size>, Size> const &initList) :
			values { } {
		for (int n = 0; n < Size; n++) {
			for (int m = 0; m < Size; m++) {
				(*this)(n, m) = initList[n][m];
			}
		}
	}
	constexpr SquareMatrix(std::initializer_list<Type> initList) {
		int i = 0;
		for (auto const &v : initList) {
			values[i++] = v;
		}
	}
	constexpr SquareMatrix(std::initializer_list<std::initializer_list<Type>> init) {
		int i = 0;
		for (auto const &row : init) {
			for (auto const &val : row) {
				values[i++] = val;
			}
		}
	}
	constexpr SquareMatrix(Type const &init) :
			values { } {
		for (int i = 0; i < size(); i++) {
			values[i] = init;
		}
	}
	constexpr SquareMatrix(SquareMatrix const &other) :
			values(other.values) {
	}
	constexpr SquareMatrix(SquareMatrix &&other) :
			values(std::move(other.values)) {
	}

	constexpr SquareMatrix& operator=(SquareMatrix const &other) {
		values = other.values;
		return *this;
	}
	constexpr SquareMatrix& operator=(SquareMatrix &&other) {
		values = std::move(other.values);
		return *this;
	}
	constexpr Type& operator()(int n, int m) {
		return values[n * Size + m];
	}
	constexpr Type const& operator()(int n, int m) const {
		return values[n * Size + m];
	}
	constexpr SquareMatrix& operator+=(SquareMatrix const &A) {
		*this = *this + A;
		return *this;
	}
	constexpr SquareMatrix& operator-=(SquareMatrix const &A) {
		*this = *this - A;
		return *this;
	}
	constexpr SquareMatrix& operator*=(Type const &a) {
		*this = *this * a;
		return *this;
	}
	constexpr SquareMatrix& operator/=(Type const &a) {
		*this = *this / a;
		return *this;
	}
	constexpr SquareMatrix operator*(Type const &a) const {
		SquareMatrix B;
		for (int k = 0; k < size(); k++) {
			B.values[k] = a * values[k];
		}
		return B;
	}
	constexpr SquareMatrix operator/(Type const &a) const {
		constexpr Type one = Type(1);
		SquareMatrix B;
		Type const aInv = one / a;
		for (int k = 0; k < size(); k++) {
			B.values[k] = aInv * values[k];
		}
		return B;
	}
	constexpr SquareMatrix operator+() const {
		return *this;
	}
	constexpr SquareMatrix operator-() const {
		SquareMatrix B;
		for (int k = 0; k < size(); k++) {
			B.values[k] = -values[k];
		}
		return B;
	}
	constexpr SquareMatrix operator+(SquareMatrix const &A) const {
		SquareMatrix B;
		for (int k = 0; k < size(); k++) {
			B.values[k] = values[k] + A.values[k];
		}
		return B;
	}
	constexpr SquareMatrix operator-(SquareMatrix const &A) const {
		SquareMatrix B;
		for (int k = 0; k < size(); k++) {
			B.values[k] = values[k] - A.values[k];
		}
		return B;
	}
	constexpr bool operator==(SquareMatrix const &A) const {
		return values == A.values;
	}
	constexpr bool operator!=(SquareMatrix const &A) const {
		return !(*this == A);
	}
	static constexpr SquareMatrix identity() {
		SquareMatrix I;
		for (int n = 0; n < Size; n++) {
			for (int m = 0; m < Size; m++) {
				I(n, m) = Type(n == m);
			}
		}
		return I;
	}
	static constexpr SquareMatrix zero() {
		SquareMatrix Z;
		std::fill(Z.begin(), Z.end(), Type(0));
		return Z;
	}
	friend SquareMatrix operator*(Type const &a, SquareMatrix const &B) {
		return B * a;
	}
private:
	std::array<Type, size()> values;
};


template<typename T, int N, int L>
SquareMatrix<T, N> operator*(SquareMatrix<T, N> const &A, SquareMatrix<T, N> const &B) {
	SquareMatrix<T, N> C;
	for (int n = 0; n < N; n++) {
		for (int l = 0; l < N; l++) {
			C(n, l) = T(0);
			for (int m = 0; m < N; m++) {
				C(n, l) += A(n, m) * B(m, l);
			}
		}
	}
	return C;
}

template<typename T, int N>
SquareMatrix<T, N> operator*=(SquareMatrix<T, N> &A, SquareMatrix<T, N> const &C) {
	SquareMatrix<T, N> B = A;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			A(i, j) = T(0);
			for (int k = 0; k < N; k++) {
				A(i, j) += B(i, k) * C(k, j);
			}
		}
	}
	return A;
}

template<typename Type, int N>
constexpr SquareMatrix<Type, N> matrixTranspose(SquareMatrix<Type, N> const &B) {
	SquareMatrix<Type, N> A;
	for (int r = 0; r < N; r++) {
		for (int c = 0; c < N; c++) {
			A(c, r) = B(r, c);
		}
	}
	return A;
}

template<typename T, int N>
constexpr SquareMatrix<T, N - 1> subMatrix(SquareMatrix<T, N> const &A, int row, int col) {
	SquareMatrix<T, N - 1> M;
	for (int r = 0; r < row; r++) {
		for (int c = 0; c < col; c++)
			M(r, c) = A(r, c);
		for (int c = col; c < N - 1; c++)
			M(r, c) = A(r, c + 1);
	}
	for (int r = row; r < N - 1; r++) {
		for (int c = 0; c < col; c++)
			M(r, c) = A(r + 1, c);
		for (int c = col; c < N - 1; c++)
			M(r, c) = A(r + 1, c + 1);
	}
	return M;
}

template<typename T, int N>
constexpr T matrixDeterminant(SquareMatrix<T, N> const &A);

template<typename T, int N>
constexpr T matrixCofactor(SquareMatrix<T, N> const &A, int r, int c) {
	if constexpr (N > 1) {
		T sgn = ((r + c) & 1) ? -T(1) : T(1);
		return sgn * matrixDeterminant(subMatrix(A, r, c));
	} else {
		return A(0, 0);
	}
}

template<typename T, int N>
constexpr SquareMatrix<T, N> matrixCofactorMatrix(SquareMatrix<T, N> const &A) {
	SquareMatrix<T, N> C;
	for (int r = 0; r < N; r++) {
		for (int c = 0; c < N; c++) {
			C(r, c) = matrixCofactorMatrix(A, r, c);
		}
	}
	return C;
}

template<typename T, int N>
constexpr T matrixDeterminant(SquareMatrix<T, N> const &A) {
	if constexpr (N > 1) {
		T sum = T(0);
		for (int c = 0; c < N; c++) {
			sum += A(0, c) * matrixCofactor(A, 0, c);
		}
		return sum;
	} else {
		return A(0, 0);
	}
}

template<typename T, int N>
constexpr SquareMatrix<T, N> matrixAdjoint(SquareMatrix<T, N> const &A) {
	return matrixTranspose(matrixCofactorMatrix(A));
}

template<typename T, int N>
constexpr auto matrixInverse(SquareMatrix<T, N> const &A) {
	return matrixAdjoint(A) / matrixDeterminant(A);
}


template<typename T, int N>
std::string toString(SquareMatrix<T, N> const &M) {
	auto formatEntry = [](T const &value) {
		std::ostringstream oss;
		oss << std::scientific << std::setprecision(3) << value;
		return oss.str();
	};
	int maxLen = 0;
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			std::string const entry = formatEntry(M(i, j));
			maxLen = std::max(maxLen, static_cast<int>(entry.size()));
		}
	}
	std::string line((N * (maxLen + 3) + 1), '-');
	line += "\n";
	std::string out = line;
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			std::string cellStr = formatEntry(M(i, j));
			while ((int) cellStr.size() < maxLen) {
				cellStr = " " + cellStr;
			}
			out += "| " + cellStr + " ";
		}
		out += "|\n" + line;
	}

	return out;
}

