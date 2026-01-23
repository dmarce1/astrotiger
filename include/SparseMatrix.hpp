#pragma once

#include <concepts>
#include <map>
#include <numeric>
#include <array>
#include <algorithm>
#include <ostream>
#include <sstream>
#include <string>

#include "SparseVector.hpp"

template<typename T>
struct IsSparseMatrix: public std::false_type {
};

template<typename T>
struct IsSparseMatrix<SparseMatrix<T>> : public std::true_type {
};

template<typename T>
concept SparseMatrixType = IsSparseMatrix<T>::value;

template<typename T>
struct SparseMatrix {
	constexpr SparseMatrix() = default;
	constexpr SparseMatrix(SparseMatrix const&) = default;
	constexpr SparseMatrix(SparseMatrix&&) = default;
	constexpr ~SparseMatrix() = default;
	constexpr SparseMatrix& operator=(SparseMatrix const&) = default;
	constexpr SparseMatrix& operator=(SparseMatrix&&) = default;
	struct reference {
		constexpr reference(SparseMatrix<T> &ref, int row, int col) :
				ref_(ref), row_(row), col_(col) {
		}
		constexpr reference& operator=(T value) {
			ref_.rows_[row_][col_] = value;
			return *this;
		}
		constexpr reference& operator*=(T value) {
			ref_.rows_[row_][col_] *= value;
			return *this;
		}
		constexpr reference& operator+=(T value) {
			ref_.rows_[row_][col_] += value;
			return *this;
		}
		constexpr reference& operator-=(T value) {
			ref_.rows_[row_][col_] -= value;
			return *this;
		}
		constexpr operator T() const {
			return ref_.rows_[row_][col_];
		}
	private:
		SparseMatrix<T> &ref_;
		int row_;
		int col_;
	};
	friend struct reference;
	constexpr SparseMatrix(int n, int m) :
			rows_(n), N_(n), M_(m) {
	}
	constexpr SparseMatrix(int n) :
			rows_(n), N_(n), M_(n) {
	}
	template<size_t L>
	constexpr SparseMatrix(int n, int m, std::array<std::pair<std::pair<int, int>, T>, L> lits) :
			rows_(n), N_(n), M_(m) {
		for (auto const &lit : lits) {
			(*this)(lit.first.first, lit.first.second) = lit.second;
		}
	}
	constexpr SparseMatrix& operator+=(SparseMatrix const &other) {
		for (int r = 0; r < N_; r++) {
			rows_[r] += other.rows_[r];
		}
		return *this;
	}
	constexpr SparseMatrix& operator-=(SparseMatrix const &other) {
		for (int r = 0; r < N_; r++) {
			rows_[r] -= other.rows_[r];
		}
		return *this;
	}
	constexpr SparseMatrix& operator*=(SparseMatrix const &other) {
		*this = *this * other;
		return *this;
	}
	constexpr SparseMatrix& operator*=(T scale) {
		for (auto &row : rows_) {
			row *= scale;
		}
		return *this;
	}
	constexpr SparseMatrix& operator/=(T scale) {
		*this *= (one / scale);
		return *this;
	}
	constexpr reference operator()(int row, int col) {
		return reference(*this, row, col);
	}
	constexpr T operator()(int row, int col) const {
		return rows_[row][col];
	}
	constexpr SparseMatrix operator+() const {
		return *this;
	}
	constexpr SparseMatrix operator-() const {
		return *this * (-one);
	}
	constexpr SparseMatrix operator+(SparseMatrix const &other) const {
		auto result = *this;
		result += other;
		return result;
	}
	constexpr SparseMatrix operator-(SparseMatrix const &other) const {
		auto result = *this;
		result -= other;
		return result;
	}
	constexpr SparseMatrix operator*(T scale) const {
		auto result = *this;
		result *= scale;
		return result;
	}
	friend constexpr SparseMatrix operator*(SparseMatrix const &A, SparseMatrix const &B) {
		SparseMatrix < T > C(A.N_, B.M_);
		auto const trB = transpose(B);
		for (int n = 0; n < A.N_; n++) {
			for (int m = 0; m < B.M_; m++) {
				C.rows_[n].set(m, A[n] * trB[m]);
			}
		}
		return C;
	}
	friend constexpr SparseMatrix operator*(T scale, SparseMatrix const &a) {
		return a * scale;
	}
	constexpr SparseMatrix operator/(T scale) const {
		return *this * (one / scale);
	}
	constexpr int density() const {
		int count = 0;
		for (int i = 0; i < N_; i++) {
			count += rows_[i].density();
		}
		return count;
	}
	constexpr int sparsity() const {
		return N_ * M_ - density();
	}
	constexpr auto rowCount() const {
		return N_;
	}
	constexpr auto colCount() const {
		return M_;
	}
	constexpr bool isSquare() const {
		return N_ == M_;
	}
	friend constexpr SparseMatrix transpose(SparseMatrix const &A) {
		SparseMatrix<T> trA(A.M_, A.N_);
		auto const &src = A.rows_;
		auto &dst = trA.rows_;
		for (int ri = 0; ri < A.N_; ri++) {
			for (auto ci = src[ri].begin(); ci != src[ri].end(); ci++) {
				dst[ci->first][ri] = ci->second;
			}
		}
		return trA;
	}
	friend constexpr auto rankReduce(SparseMatrix<T>& A, SparseMatrixType auto&... B) {
		auto const N = A.rowCount();
		int p = 0;
		for (int n = 0; n < N; n++) {
			for (int j = p; j < N; j++) {
				if (A(j, n) != T(0)) {
					(std::swap(B[j], B[p]),...);
					std::swap(A[j], A[p]);
					break;
				}
			}
			if (A(p, n) != T(0)) {
				((B[p] /= A(p, n)),...);
				A[p] /= A(p, n);
				for (int j = p + 1; j < N; j++) {
					if (A(j, n) != T(0)) {
						((B[j] -= A(j, n) * B[p]),...);
						A[j] -= A(j, n) * A[p];
					}
				}
				p++;
			}
		}
		A.rows_.resize(p);
		A.N_ = p;
		return p;
	}
	friend constexpr void gaussJordanElimination(SparseMatrix& A, SparseMatrixType auto&... B) {
		rankReduce(A, B...);
		for (int n = A.rowCount() - 1; n >= 0; n--) {
			for (int j = 0; j < n; j++) {
				if (A(j, n) != zero) {
					((B[j] -= A(j, n) * B[n]),...);
					A[j] -= A(j, n) * A[n];
				}
			}
		}
	}
	friend constexpr std::ostream& operator<<(std::ostream &os, SparseMatrix<T> const &A) {
		NUMERICAL_CONSTANTS(T);
		auto const N = A.rowCount();
		auto const M = A.colCount();
		SparseMatrix<int> widths(N, M);
		int maxWidth = 1;
		for (int n = 0; n < N; n++) {
			for (int m = 0; m < M; m++) {
				if (A(n, m) != zero) {
					std::ostringstream ss;
					ss << A(n, m);
					widths(n, m) = ss.str().size();
					maxWidth = std::max(maxWidth, (int) ss.str().size());
				}
			}
		}
		auto printHorizontal = [&os, maxWidth, M]() {
			for (int m = 0; m < M; m++) {
				os << '+' << std::string(maxWidth, '-');
			}
			os << "+\n";
		};
		for (int n = 0; n < N; n++) {
			printHorizontal();
			for (int m = 0; m < M; m++) {
				os << '|';
				if (A(n, m) != zero) {
					int const padding = maxWidth - widths(n, m);
					os << std::string(padding / 2, ' ') << A(n, m) << std::string((padding + 1) / 2, ' ');
				} else {
					os << std::string(maxWidth, ' ');
				}
			}
			os << "|\n";
		}
		printHorizontal();
		return os;
	}
	constexpr auto literal() const {
		std::vector<std::pair<std::pair<int, int>, T>> lits;
		lits.reserve(density());
		for (int ri = 0; ri < N_; ri++) {
			auto const &row = rows_[ri];
			for (auto it = row.begin(); it != row.end(); it++) {
				lits.push_back( { { ri, it->first }, it->second });
			}
		}
		return lits;
	}
	static constexpr SparseMatrix identity(int n) {
		SparseMatrix I(n, n);
		for (int m = 0; m < n; m++) {
			I(m, m) = one;
		}
		return I;
	}
	friend struct reference;
private:
	constexpr SparseVector<T>& operator[](int i) {
		return rows_[i];
	}
	constexpr SparseVector<T> operator[](int i) const {
		return rows_[i];
	}
	NUMERICAL_CONSTANTS (T);
	std::vector<SparseVector<T>> rows_ { };
	int N_ { };
	int M_ { };
};

template<typename T>
constexpr SparseMatrix<T> inverse(SparseMatrix<T> A) {
	ASSERT(A.isSquare());
	auto iA = SparseMatrix<T>::identity(A.rowCount());
	gaussJordanElimination(A, iA);
	return iA;
}

template<typename T>
constexpr SparseMatrix<T> leftPseudoinverse(SparseMatrix<T> A) {
	auto const trA = transpose(A);
	return inverse(trA * A) * trA;
}

template<typename T>
constexpr SparseMatrix<T> rightPseudoinverse(SparseMatrix<T> A) {
	auto const trA = transpose(A);
	return trA * inverse(A * trA);
}

