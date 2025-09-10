/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
*******************************************************************************/

#ifndef INCLUDE_MULTIDICES_HPP_
#define INCLUDE_MULTIDICES_HPP_


#include <cassert>
#include "Util.hpp"

template<int D>
struct Multidices {
	constexpr Multidices() {
		indices_.fill(0);
	}
	friend constexpr int abs(Multidices<D> const &i) {
		return std::accumulate(i.indices_.begin(), i.indices_.end(), 0);
	}
	constexpr int& operator[](int i) {
		return indices_[i];
	}
	constexpr 	int operator[](int i) const {
		return indices_[i];
	}
	constexpr Multidices(Multidices const&) = default;
	constexpr Multidices(Multidices&&) = default;
	constexpr Multidices(int i) {
		*this = i;
	}
	constexpr Multidices& operator=(Multidices const&) = default;
	constexpr Multidices& operator=(Multidices&&) = default;
	constexpr Multidices& operator=(int i) {
		auto &A = *this;
		for (int d = 0; d < D - 1; d++) {
			A[d] = 0;
		}
		A[D - 1] = i;
		return *this;
	}
	constexpr Multidices operator+(Multidices const &other) const {
		Multidices A;
		for (int d = 0; d < D; d++) {
			A[d] = (*this)[d] + other[d];
		}
		return A;
	}
	constexpr Multidices operator*(Multidices other) const {
		Multidices A = *this;
		for (int d = 0; d < D; d++) {
			A[d] *= other[d];
		}
		return A;
	}
	constexpr Multidices operator*(int other) const {
		Multidices A = *this;
		for (int d = 0; d < D; d++) {
			A[d] *= other;
		}
		return A;
	}
	constexpr Multidices operator-(Multidices const &other) const {
		Multidices A;
		for (int d = 0; d < D; d++) {
			A[d] = (*this)[d] - other[d];
		}
		return A;
	}
	constexpr Multidices& operator+=(Multidices const &other) {
		*this = *this + other;
		return *this;

	}
	constexpr Multidices& operator-=(Multidices const &other) {
		*this = *this - other;
		return *this;

	}
	constexpr Multidices& operator*=(Multidices const &other) {
		*this = *this * other;
		return *this;

	}
	constexpr bool operator==(Multidices const &other) const {
		for (int d = 0; d < D; d++) {
			if (indices_[d] != other[d]) {
				return false;
			}
		}
		return true;
	}
	constexpr bool operator!=(Multidices const &other) const {
		return !(*this == other);
	}
	constexpr bool operator<(Multidices const &other) const {
		for (int d = 0; d < D; d++) {
			if (indices_[d] >= other[d]) {
				return false;
			}
		}
		return true;
	}
	constexpr bool operator<=(Multidices const &other) const {
		for (int d = 0; d < D; d++) {
			if (indices_[d] > other[d]) {
				return false;
			}
		}
		return true;
	}
	constexpr bool operator<(int const &other) const {
		for (int d = 0; d < D; d++) {
			if (indices_[d] >= other) {
				return false;
			}
		}
		return true;
	}
	constexpr bool operator<=(int const &other) const {
		for (int d = 0; d < D; d++) {
			if (indices_[d] > other) {
				return false;
			}
		}
		return true;
	}
	constexpr Multidices& operator++() {
		auto &A = *this;
		for (int d = 0; d < D - 1; d++) {
			A[d + 1] += A[d];
		}
		if constexpr (D == 1) {
			A[0]++;
		} else {
			int d = 0;
			while (++A[d] > A[d + 1]) {
				A[d] = 0;
				d++;
				if (d == D - 1) {
					for (int d = 0; d < D - 1; d++) {
						A[d] = 0;
					}
					A[D - 1]++;
					break;
				}
			}
		}
		for (int d = D - 2; d >= 0; d--) {
			A[d + 1] -= A[d];
		}
		return *this;
	}
	constexpr Multidices& operator--() {
		auto &A = *this;
		for (int d = 0; d < D - 1; d++) {
			A[d + 1] += A[d];
		}
		int d = 0;
		while (A[d]-- <= 0) {
			A[d] = 0;
			d++;
			if (d == D - 1) {
				int const next = std::max(0, A[D - 1] - 1);
				for (int d = 0; d < D; d++) {
					A[d] = next;
				}
				break;
			}
		}
		for (int d = D - 2; d >= 0; d--) {
			A[d + 1] -= A[d];
		}
		return *this;
	}
	constexpr Multidices operator++(int) {
		auto const tmp = *this;
		operator++();
		return tmp;
	}
	constexpr Multidices operator--(int) {
		auto const tmp = *this;
		operator--();
		return tmp;
	}
	constexpr operator int() const {
		auto A = *this;
		int flat = 0;
		for (int d = 0; d < D - 1; d++) {
			A[d + 1] += A[d];
		}
		for (int d = 1; d < D; d++) {
			assert(A[d - 1] <= A[d]);
		}
		for (int d = D - 1; d >= 0; d--) {
			flat += binco(A[d] + d, 1 + d);
		}
		return flat;
	}
	template<typename T>
	friend constexpr std::array<T, D> pow(std::array<T, D> const &x, Multidices const &n) {
		std::array<T, D> result;
		for (int d = 0; d < D; d++) {
			result[d] = ipow(x[d], n[d]);
		}
		return result;
	}
	friend constexpr unsigned long long factorial(Multidices const &n) {
		unsigned long long x = 1;
		for (int d = 0; d < D; d++) {
			for (long long z = 1; z <= n[d]; z++) {
				x *= z;
			}
		}
		return x;
	}
	friend constexpr unsigned long long binco(Multidices const &n, Multidices const &k) {
		unsigned long long result = 1;
		for (int d = 0; d < D; d++) {
			if (k[d] > n[d]) {
				return 0;
			}
			result *= binco(n[d], k[d]);
		}
		return result;
	}
	friend std::ostream& operator<<(std::ostream &os, Multidices const &A) {
		os << "(";
		for (int d = 0; d < D - 1; d++) {
			os << A[d] << ", ";
		}
		os << A[D - 1] << ")";
		return os;
	}
	static constexpr Multidices unit(int dim = 0) {
		Multidices I;
		I[dim] = 1;
		return I;
	}
private:
	std::array<int, D> indices_;
};

template<int N>
struct std::hash<Multidices<N>> {
	size_t operator()(Multidices<N> const &indices) const {
		constexpr size_t n = ((size_t(1) << size_t(31)) - size_t(1));
		constexpr size_t a = 48271;
		constexpr size_t b = 16807;
		std::hash<int> keyGen;
		size_t key = 42;
		for (int i = 0; i < N; i++) {
			key = ((a * key) % n) ^ ((b * keyGen(i)) % n);
		}
		return key;
	}
};


#endif /* INCLUDE_MULTIDICES_HPP_ */
