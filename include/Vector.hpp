#pragma once

#include <array>
#include <cmath>
#include <initializer_list>
#include <ostream>

template <typename T, int N>
struct Vector {
	constexpr T operator[](int i) const {
		return data_[i];
	}
	constexpr T &operator[](int i) {
		return data_[i];
	}
	constexpr Vector() = default;
	constexpr Vector(std::initializer_list<T> const &list) :
		data_(list) {
	}
	constexpr Vector(std::array<T, N> const &array) :
		data_(array) {
	}
	constexpr Vector(Vector const &) = default;
	constexpr Vector(Vector &&) = default;
	constexpr Vector &operator=(Vector const &) = default;
	constexpr Vector &operator=(Vector &&) = default;
	constexpr Vector &operator+=(Vector const &B) {
		*this = *this + B;
		return *this;
	}
	constexpr Vector &operator-=(Vector const &B) {
		*this = *this - B;
		return *this;
	}
	constexpr Vector &operator*=(T const &B) {
		*this = *this * B;
		return *this;
	}
	constexpr Vector &operator/=(T const &B) {
		return (*this *= one / B);
	}
	constexpr T dot(Vector const &A) const {
		auto const &B = *this;
		T sum = zero;
		for (int k = 0; k < N; k++) {
			sum += A[k] * B[k];
		}
		return sum;
	}
	constexpr T max() const {
		using std::max;
		auto const &A = *this;
		T m = A[0];
		for (int k = 1; k < N; k++) {
			m = max(m, A[k]);
		}
		return m;
	}
	constexpr T min() const {
		using std::min;
		auto const &A = *this;
		T m = A[0];
		for (int k = 1; k < N; k++) {
			m = min(m, A[k]);
		}
		return m;
	}
	static constexpr Vector unit(int j) {
		Vector u;
		for (int i = 0; i < N; i++) {
			u[i] = T(j == i);
		}
		return u;
	}
	friend constexpr Vector operator+(Vector const &A) {
		return A;
	}
	friend constexpr Vector operator-(Vector A) {
		for (int i = 0; i < N; i++) {
			A[i] = -A[i];
		}
		return A;
	}
	friend constexpr Vector operator+(Vector A, Vector const &B) {
		for (int i = 0; i < N; i++) {
			A[i] += B[i];
		}
		return A;
	}
	friend constexpr Vector operator-(Vector A, Vector const &B) {
		for (int i = 0; i < N; i++) {
			A[i] -= B[i];
		}
		return A;
	}
	friend constexpr Vector operator*(Vector A, T const &B) {
		for (int i = 0; i < N; i++) {
			A[i] *= B;
		}
		return A;
	}
	friend constexpr Vector operator*(T const &A, Vector B) {
		return B * A;
	}
	friend constexpr Vector operator/(Vector A, T const &B) {
		return A * (one / B);
	}
	friend constexpr T sqr(Vector const &A) {
		return A.dot(A);
	}
	friend constexpr T abs(Vector const &A) {
		using std::sqrt;
		return sqrt(sqr(A));
	}
	friend constexpr Vector norm(Vector A) {
		A /= abs(A);
		return A;
	}
	constexpr operator std::array<T, N>() const {
		return data_;
	};

private:
	static constexpr T zero = T(0);
	static constexpr T one = T(1);
	std::array<T, N> data_;
};
