/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_FWDAUTODIFF_HPP_
#define INCLUDE_FWDAUTODIFF_HPP_

#include <type_traits>
#include "Matrix.hpp"

template<typename T1>
struct FwdAutoDiff;

// ---------- Traits ----------
template<typename T>
struct IsFwdAutoDiff {
	static constexpr bool value = false;
	using type = void;
};

template<typename T>
struct IsFwdAutoDiff<FwdAutoDiff<T>> {
	static constexpr bool value = true;
	using type = T;
};

template<typename T>
inline constexpr bool IsFwdAutoDiff_v = IsFwdAutoDiff<T>::value;

template<typename T>
using FwdAutoDiffValue_t = typename IsFwdAutoDiff<T>::type;

template<typename T1>
struct FwdAutoDiff {
private:
	T1 v_, d1_, d2_;
	constexpr FwdAutoDiff(T1 value, T1 derivative, T1 derivative2) :
			v_(value), d1_(derivative), d2_(derivative2) {
	}
public:
	using value_type = T1;

	explicit constexpr operator T1() const {
		return v_;
	}

	constexpr T1 D(int d = 1) const {
		if (d == 2) {
			return d2_;
		}
		if (d == 0) {
			return v_;
		}
		return d1_;
	}

	friend constexpr FwdAutoDiff operator+(FwdAutoDiff const &A) {
		return A;
	}

	friend constexpr FwdAutoDiff operator-(FwdAutoDiff const &A) {
		FwdAutoDiff C;
		C.v_ = -A.v_;
		C.d1_ = -A.d1_;
		C.d2_ = -A.d2_;
		return C;
	}

	friend constexpr FwdAutoDiff operator+(FwdAutoDiff const &A, FwdAutoDiff const &B) {
		FwdAutoDiff C;
		C.v_ = A.v_ + B.v_;
		C.d1_ = A.d1_ + B.d1_;
		C.d2_ = A.d2_ + B.d2_;
		return C;
	}

	friend constexpr FwdAutoDiff operator-(FwdAutoDiff const &A, FwdAutoDiff const &B) {
		FwdAutoDiff C;
		C.v_ = A.v_ - B.v_;
		C.d1_ = A.d1_ - B.d1_;
		C.d2_ = A.d2_ - B.d2_;
		return C;
	}

	// --- Friend declarations (no definitions here) ---
	template<typename A, typename B>
	friend constexpr auto operator*(FwdAutoDiff<A> const&, FwdAutoDiff<B> const&);

	template<typename A, typename U, typename >
	friend constexpr FwdAutoDiff<decltype(std::declval<A>() * std::declval<U>())> operator*(FwdAutoDiff<A> const&, U const&);

	template<typename U, typename A, typename >
	friend constexpr FwdAutoDiff<decltype(std::declval<A>() * std::declval<U>())> operator*(U const&, FwdAutoDiff<A> const&);

	template<typename A, typename B>
	friend constexpr auto operator/(FwdAutoDiff<A> const&, FwdAutoDiff<B> const&);

	template<typename A, typename U>
	friend constexpr auto operator/(FwdAutoDiff<A> const&, U const&);

	template<typename U, typename A>
	friend constexpr auto operator/(U const&, FwdAutoDiff<A> const&);

	// compound ops can stay members (they use the free operators)
	template<typename Rhs>
	constexpr FwdAutoDiff& operator*=(FwdAutoDiff<Rhs> const &rhs) {
		*this = *this * rhs;
		return *this;
	}
	template<typename U, typename = std::enable_if_t<!IsFwdAutoDiff<U>::value>>
	constexpr FwdAutoDiff& operator*=(U const &rhs) {
		*this = *this * rhs;
		return *this;
	}
	template<typename Rhs>
	constexpr FwdAutoDiff& operator/=(FwdAutoDiff<Rhs> const &rhs) {
		*this = *this / rhs;
		return *this;
	}
	template<typename U, typename = std::enable_if_t<!IsFwdAutoDiff<U>::value>>
	constexpr FwdAutoDiff& operator/=(U const &rhs) {
		*this = *this / rhs;
		return *this;
	}

	constexpr FwdAutoDiff& operator+=(FwdAutoDiff const &A) {
		*this = *this + A;
		return *this;
	}

	constexpr FwdAutoDiff& operator-=(FwdAutoDiff const &A) {
		*this = *this - A;
		return *this;
	}

	friend constexpr auto sqr(FwdAutoDiff const &x) {
		return x * x;
	}

	template<typename T2>
	friend constexpr auto pow(FwdAutoDiff const &x, T2 const &n) {
		using std::pow;
		FwdAutoDiff A;
		A.v_ = pow(x.v_, n);
		A.d1_ = n * x.d1_ * A.v_ / x.v_;
		A.d2_ = n * (x.d2_ * A.v_ / x.v_ + (n - T1(1)) * sqr(x.d1_ / x.v_) * A.v_);
		return x * x;
	}

	friend constexpr FwdAutoDiff sqrt(FwdAutoDiff const &B) {
		using std::sqrt;
		FwdAutoDiff A;
		A.v_ = sqrt(B.v_);
		A.d1_ = 0.5 * B.d1_ / A.v_;
		A.d2_ = (2 * B.v_ * B.d2_ - sqr(B.d1_)) / (4 * B.v_ * A.v_);
		return A;
	}

	constexpr FwdAutoDiff() {
		v_ = T1(0.0);
		d1_ = T1(0.0);
		d2_ = T1(0.0);
	}

	constexpr FwdAutoDiff(T1 value) {
		v_ = value;
		d1_ = T1(0.0);
		d2_ = T1(0.0);
	}

	template<typename >
	friend class FwdAutoDiff;

	static constexpr FwdAutoDiff independentVariable(T1 const &var) {
		FwdAutoDiff a;
		a.v_ = var;
		a.d1_ = T1(1.0);
		a.d2_ = T1(0.0);
		return a;
	}

	static constexpr FwdAutoDiff constant(T1 const &var) {
		return FwdAutoDiff(var, 0.0, 0.0);
	}
};

template<typename A, typename B>
constexpr auto operator*(FwdAutoDiff<A> const &X, FwdAutoDiff<B> const &Y) {
	using R = decltype(std::declval<A>() * std::declval<B>());
	FwdAutoDiff<R> Z;
	Z.v_ = X.v_ * Y.v_;
	Z.d1_ = X.d1_ * Y.v_ + X.v_ * Y.d1_;
	Z.d2_ = X.d2_ * Y.v_ + X.v_ * Y.d2_ + (X.d1_ * Y.d1_) * 2; // 2 is dimensionless
	return Z;
}

// AD * scalar-like (U is not FwdAutoDiff)
template<typename A, typename U, typename = std::enable_if_t<!IsFwdAutoDiff<U>::value>>
constexpr FwdAutoDiff<decltype(std::declval<A>() * std::declval<U>())> operator*(FwdAutoDiff<A> const &X, U const &b) {
	using R = decltype(std::declval<A>() * std::declval<U>());
	FwdAutoDiff<R> Z;
	Z.v_ = X.v_ * b;
	Z.d1_ = X.d1_ * b;
	Z.d2_ = X.d2_ * b;
	return Z;
}

template<typename U, typename A, typename = std::enable_if_t<!IsFwdAutoDiff_v<U>>>
constexpr FwdAutoDiff<decltype(std::declval<U>() * std::declval<A>())> operator*(U const &a, FwdAutoDiff<A> const &Y) {
	return Y * a;
}

// ==========================
// Division
// ==========================

// AD / AD
template<typename A, typename B>
constexpr auto operator/(FwdAutoDiff<A> const &X, FwdAutoDiff<B> const &Y) {
	using R = decltype(std::declval<A>() / std::declval<B>());
	FwdAutoDiff<R> Z;
	const auto &A0 = X.v_;
	const auto &A1 = X.d1_;
	const auto &A2 = X.d2_;
	const auto &B0 = Y.v_;
	const auto &B1 = Y.d1_;
	const auto &B2 = Y.d2_;

	Z.v_ = A0 / B0;
	Z.d1_ = (A1 * B0 - A0 * B1) / (B0 * B0);

	// y = A/B.  y'' = (A'' B^2 - 2 A' B B' - A B B'' + 2 A (B')^2) / B^3
	Z.d2_ = (A2 * B0 * B0 - (A1 * B0 * B1) * 2 - A0 * B0 * B2 + A0 * (B1 * B1) * 2) / (B0 * B0 * B0);
	return Z;
}

// AD / scalar-like (U not FwdAutoDiff)
template<typename A, typename U, typename = std::enable_if_t<!IsFwdAutoDiff_v<U>>>
constexpr auto operator/(FwdAutoDiff<A> const &X, U const &b) -> FwdAutoDiff<decltype(std::declval<A>() / std::declval<U>())>
{
	using R = decltype(std::declval<A>() / std::declval<U>());
	FwdAutoDiff<R> Z;
	Z.v_ = X.v_ / b;
	Z.d1_ = X.d1_ / b;
	Z.d2_ = X.d2_ / b;
	return Z;
}

// scalar-like / AD
template<typename U, typename A, typename = std::enable_if_t<!IsFwdAutoDiff_v<U>>>
constexpr auto operator/(U const &a, FwdAutoDiff<A> const &Y) -> FwdAutoDiff<decltype(std::declval<U>() / std::declval<A>())>
{
	using R = decltype(std::declval<U>() / std::declval<A>());
	FwdAutoDiff<R> Z;
	const auto &B0 = Y.v_;
	const auto &B1 = Y.d1_;
	const auto &B2 = Y.d2_;

	Z.v_ = a / B0;
	Z.d1_ = (-a * B1) / (B0 * B0);

	// (a / B)'' = ( 2 a (B')^2 - a B B'' ) / B^3
	Z.d2_ = ((B1 * B1) * (a * 2) - a * B0 * B2) / (B0 * B0 * B0);
	return Z;
}

#endif /* INCLUDE_FWDAUTODIFF_HPP_ */
