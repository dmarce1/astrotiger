/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
*******************************************************************************/

#ifndef INCLUDE_FWDAUTODIFF_HPP_
#define INCLUDE_FWDAUTODIFF_HPP_

#include "Matrix.hpp"

template<typename T1>
struct FwdAutoDiff;


template<typename T1>
struct IsFwdAutoDiff {
	static constexpr bool value = false;
	using type = void;
};

template<typename T1>
struct IsFwdAutoDiff<FwdAutoDiff<T1>> {
	static constexpr bool value = true;
	using type = T1;
};



template<typename T1>
struct FwdAutoDiff {
private:
	T1 v_;
	T1 d1_;
	T1 d2_;
	constexpr FwdAutoDiff(T1 value, T1 derivative, T1 derivative2) :
			v_(value), d1_(derivative), d2_(derivative2) {
	}
public:
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
	template<typename T2, typename T3 = decltype(T1() * T2())>
	constexpr FwdAutoDiff<T3> operator*(FwdAutoDiff<T2> const &B) const {
		FwdAutoDiff<T3> A;
		A.v_ = v_ * B.v_;
		A.d1_ = d1_ * B.v_ + v_ * B.d1_;
		A.d2_ = d2_ * B.v_ + 2 * d1_ * B.d1_ + v_ * B.d2_;
		return A;
	}
	constexpr FwdAutoDiff& operator*=(FwdAutoDiff const &A) {
		*this = *this * A;
		return *this;
	}
	constexpr FwdAutoDiff& operator*=(T1 const &a) {
		*this = *this * a;
		return *this;
	}
	template<typename T2, typename T3 = decltype(T1() / T2())>
	constexpr FwdAutoDiff<T3> operator/(FwdAutoDiff<T2> const &B) const {
		FwdAutoDiff<T3> A;
		A.v_ = v_ / B.v_;
		A.d1_ = (d1_ * B.v_ - v_ * B.d1_) / sqr(B.v_);
		A.d2_ = (-2 * B.v_ * d1_ * B.d1_ + 2 * v_ * sqr(B.d1_) + sqr(B.v_) * d2_ - v_ * B.v_ * B.d2_) / (B.v_ * sqr(B.v_));
		return A;
	}
	constexpr FwdAutoDiff& operator/=(FwdAutoDiff const &A) {
		*this = *this / A;
		return *this;
	}
	constexpr FwdAutoDiff& operator/=(T1 const &a) {
		*this = *this / a;
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





#endif /* INCLUDE_FWDAUTODIFF_HPP_ */
