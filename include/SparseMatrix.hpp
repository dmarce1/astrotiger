#pragma once

#include <algorithm>
#include <array>
#include <concepts>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <ostream>
#include <sstream>
#include <string>
#include <type_traits>

#include "Definitions.hpp"
#include "SparseVector.hpp"

template <typename T>
struct IsSparseMatrix : public std::false_type {};

template <typename T>
struct IsSparseMatrix<SparseMatrix<T>> : public std::true_type {};

template <typename T>
concept SparseMatrixType = IsSparseMatrix<T>::value;

template <typename T>
struct SparseMatrix {
	constexpr SparseMatrix() = default;
	constexpr SparseMatrix(SparseMatrix const &) = default;
	constexpr SparseMatrix(SparseMatrix &&) = default;
	constexpr ~SparseMatrix() = default;
	constexpr SparseMatrix &operator=(SparseMatrix const &) = default;
	constexpr SparseMatrix &operator=(SparseMatrix &&) = default;
	struct reference {
		constexpr reference(SparseMatrix<T> &ref, int row, int col) :
			ref_(ref), row_(row), col_(col) {
		}
		constexpr reference &operator=(T value) {
			ref_.rows_[row_][col_] = value;
			return *this;
		}
		constexpr reference &operator=(reference const &value) {
			ref_.rows_[row_][col_] = value;
			return *this;
		}
		constexpr reference &operator*=(T value) {
			ref_.rows_[row_][col_] *= value;
			return *this;
		}
		constexpr reference &operator+=(T value) {
			ref_.rows_[row_][col_] += value;
			return *this;
		}
		constexpr reference &operator-=(T value) {
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
		rows_(n, SparseVector<T>(m)), N_(n), M_(m) {
	}
	constexpr SparseMatrix(int n) :
		rows_(n, SparseVector<T>(n)), N_(n), M_(n) {
	}
	template <std::integral auto N, std::integral auto M>
	constexpr SparseMatrix(std::array<std::array<T, N>, M> const &a) {
		N_ = N;
		M_ = M;
		rows_.resize(N_);
		for (int n = 0; n < N; n++) {
			for (int m = 0; m < M; m++) {
				auto const anm = a[n][m];
				if (!isZero(anm)) {
					rows_[n][m] = a[n][m];
				}
			}
		}
	}
	template <size_t L>
	constexpr SparseMatrix(int n, int m, std::array<std::pair<std::pair<int, int>, T>, L> lits) :
		rows_(n, SparseVector<T>(m)), N_(n), M_(m) {
		for (auto const &lit : lits) {
			(*this)(lit.first.first, lit.first.second) = lit.second;
		}
	}
	constexpr SparseMatrix &operator+=(SparseMatrix const &other) {
		ASSERT((N_ == other.N_) && (M_ == other.M_));
		for (int r = 0; r < N_; r++) {
			rows_[r] += other.rows_[r];
		}
		return *this;
	}
	constexpr SparseMatrix &operator-=(SparseMatrix const &other) {
		ASSERT((N_ == other.N_) && (M_ == other.M_));
		for (int r = 0; r < N_; r++) {
			rows_[r] -= other.rows_[r];
		}
		return *this;
	}
	constexpr SparseMatrix &operator*=(SparseMatrix const &other) {
		ASSERT(isSquare() && (N_ == other.N_) && (M_ == other.M_));
		*this = *this * other;
		return *this;
	}

	constexpr SparseMatrix &operator*=(T scale) {
		for (auto &row : rows_) {
			row *= scale;
		}
		return *this;
	}
	constexpr SparseMatrix &operator/=(T scale) {
		*this *= (one / scale);
		return *this;
	}
	constexpr const SparseVector<T>& operator[](int i) const {
		return rows_[i];
	}
	constexpr SparseVector<T> &operator[](int i) {
		return rows_[i];
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
		ASSERT((N_ == other.N_) && (M_ == other.M_));
		auto result = *this;
		result += other;
		return result;
	}
	constexpr SparseMatrix operator-(SparseMatrix const &other) const {
		ASSERT((N_ == other.N_) && (M_ == other.M_));
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
		ASSERT(A.M_ == B.N_);
		SparseMatrix<T> C(A.N_, B.M_);
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
	constexpr void setRowCount(int n) {
		rows_.resize(n);
		N_ = n;
	}
	constexpr SparseMatrix &addRow(SparseVector<T> &&v) {
		rows_.push_back(std::move(v));
		N_++;
		return *this;
	}
	constexpr SparseMatrix &addRow() {
		addRow(SparseVector<T>(M_));
		return *this;
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
	template <typename U>
	constexpr operator SparseMatrix<U>() const {
		SparseMatrix<U> A(N_, M_);
		for (int n = 0; n < N_; n++) {
			A.rows_[n] = static_cast<SparseVector<U>>(rows_[n]);
		}
		return A;
	}
	friend constexpr auto normalize(SparseMatrix const &B) {
		SparseMatrix A(B.N_, B.M_);
		for (int n = 0; n < B.N_; n++) {
			A.rows_[n] = normalize(B.rows_[n]);
		}
		return A;
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
	constexpr T determinant() const {
		auto const N = rowCount();
		auto A = *this;
		T det = one;
		int k = 0;
		for (int n = 0; n < N; n++) {
			int p = -1;
			for (int j = k; j < N; j++) {
				if (!isZero(A(j, n))) {
					p = j;
					break;
				}
			}
			if (p == -1) {
				return zero;
			}
			if (p != k) {
				std::swap(A[p], A[k]);
				det = -det;
			}
			T const Akn = A(k, n);
			det *= Akn;
			for (int j = k + 1; j < N; j++) {
				T const Ajn = A(j, n);
				if (!isZero(Ajn)) {
					A[j] -= A[k] * (Ajn / Akn);
				}
			}
			k++;
		}
		return det;
	}
	constexpr auto pivots() const {
		auto const N = rowCount();
		auto A = *this;
		int k = 0;
		for (int n = 0; n < N; n++) {
			for (int j = k; j < N; j++) {
				T const Aⱼₙ = A(j, n);
				if (!isZero(Aⱼₙ)) {
					std::swap(A[j], A[k]);
					break;
				}
			}
			T const Aₖₙ = A(k, n);
			if (!isZero(Aₖₙ)) {
				A[k] /= Aₖₙ;
				for (int j = k + 1; j < N; j++) {
					T const Aⱼₙ = A(j, n);
					if (!isZero(Aⱼₙ)) {
						A[j] -= Aⱼₙ * A[k];
					}
				}
				k++;
			}
		}
		std::vector<int> pivots;
		pivots.reserve(N);
		for (int n = 0; n < N; n++) {
			auto const it = A[n].begin();
			if (it != A[n].end()) {
				pivots.push_back(it->first);
			}
		}
		return pivots;
	}
	constexpr auto rank() const {
		return pivots().size();
	}
	constexpr auto basis() const {
		auto const p = this->pivots();
		int const rank = p.size();
		auto const trA = transpose(*this);
		SparseMatrix φ(rank, rowCount());
		for (int n = 0; n < rank; n++) {
			φ[n] = trA[p[n]];
		}
		φ = transpose(φ);
		return φ;
	}
	friend constexpr void gaussJordanElimination(SparseMatrix &A, SparseMatrixType auto &B) {
		auto const N = A.rowCount();
		int k = 0;
		for (int n = 0; n < N; n++) {
			for (int j = k; j < N; j++) {
				T const Aⱼₙ = A(j, n);
				if (!isZero(Aⱼₙ)) {
					std::swap(B[j], B[k]);
					std::swap(A[j], A[k]);
					break;
				}
			}
			T const Aₖₙ = A(k, n);
			if (!isZero(Aₖₙ)) {
				B[k] /= Aₖₙ;
				A[k] /= Aₖₙ;
				for (int j = k + 1; j < N; j++) {
					T const Aⱼₙ = A(j, n);
					if (!isZero(Aⱼₙ)) {
						B[j] -= Aⱼₙ * B[k];
						A[j] -= Aⱼₙ * A[k];
					}
				}
				k++;
			}
		}
		for (int n = k - 1; n >= 0; n--) {
			for (int j = 0; j < n; j++) {
				T const Aⱼₙ = A(j, n);
				if (!isZero(Aⱼₙ)) {
					B[j] -= Aⱼₙ * B[n];
					A[j] -= Aⱼₙ * A[n];
				}
			}
		}
	}
	friend constexpr SparseVector<T> operator*(SparseMatrix<T> const &A, SparseVector<T> const &x) {
		ASSERT(A.colCount() == x.size());
		auto const N = A.rowCount();
		SparseVector<T> v(A.rowCount());
		for (int n = 0; n < N; n++) {
			v[n] = A.rows_[n] * x;
		}
		return v;
	}
	friend constexpr SparseVector<T> operator*(SparseVector<T> const &x, SparseMatrix<T> const &A) {
		ASSERT(A.rowCount() == x.size());
		return transpose(A) * x;
	}
	friend constexpr std::ostream &operator<<(std::ostream &os, SparseMatrix<T> const &A) {
		NUMERICAL_CONSTANTS(T);
		auto const N = A.rowCount();
		auto const M = A.colCount();
		SparseMatrix<int> widths(N, M);
		int maxWidth = 1;
		os << N << "x" << M << std::endl;
		for (int n = 0; n < N; n++) {
			for (int m = 0; m < M; m++) {
				if (A(n, m) != zero) {
					std::ostringstream ss;
					ss << A(n, m);
					widths(n, m) = ss.str().size();
					maxWidth = std::max(maxWidth, (int)ss.str().size());
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
				lits.push_back({{ri, it->first}, it->second});
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
	template <typename>
	friend struct SparseMatrix;
	friend struct reference;

private:
	static constexpr bool isZero(T const &value) {
		if constexpr (std::is_floating_point<T>::value) {
			constexpr T ε = std::sqrt(std::numeric_limits<T>::epsilon());
			return ((-ε < value) && (value < +ε));
		} else {
			return value == zero;
		}
	}
	NUMERICAL_CONSTANTS(T);
	std::vector<SparseVector<T>> rows_{};
	int N_{};
	int M_{};
};

template <typename T>
constexpr SparseMatrix<T> inverse(SparseMatrix<T> A) {
	ASSERT(A.isSquare());
	auto iA = SparseMatrix<T>::identity(A.rowCount());
	gaussJordanElimination(A, iA);
	return iA;
}

template <typename T>
constexpr SparseMatrix<T> leftPseudoinverse(SparseMatrix<T> A) {
	auto const trA = transpose(A);
	return inverse(trA * A) * trA;
}

template <typename T>
constexpr SparseMatrix<T> rightPseudoinverse(SparseMatrix<T> A) {
	auto const trA = transpose(A);
	return trA * inverse(A * trA);
}

template <typename T>
constexpr auto kroneckerProduct(SparseMatrix<T> const &A, SparseMatrix<T> const &B) {
	int const N1 = A.rowCount();
	int const M1 = A.colCount();
	int const N2 = B.rowCount();
	int const M2 = B.colCount();
	SparseMatrix<T> C(N1 * N2, M1 * M2);
	for (int n1 = 0; n1 < N1; n1++) {
		for (int m1 = 0; m1 < M1; m1++) {
			for (int n2 = 0; n2 < N2; n2++) {
				for (int m2 = 0; m2 < M2; m2++) {
					C(N2 * n1 + n2, M2 * m1 + m2) = A(n1, m1) * B(n2, m2);
				}
			}
		}
	}
	return C;
}
