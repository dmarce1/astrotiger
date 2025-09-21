/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_UNITS_HPP_
#define INCLUDE_UNITS_HPP_

#include <array>
#include <cmath>
#include <cstdint>
#include <type_traits>

enum class UnitType : int {
	length, mass, time, temperature
};

template<int L, int M, int T, int K>
struct Unit {
	static constexpr int L_ = L;
	static constexpr int M_ = M;
	static constexpr int T_ = T;
	static constexpr int K_ = K;
};

template<typename Type>
struct isUnit {
	static constexpr bool value = false;
};

template<int L, int M, int T, int K>
struct isUnit<Unit<L, M, T, K>> {
	static constexpr bool value = true;
};

template<typename UnitA, int power, std::enable_if_t<isUnit<UnitA>::value, int> = 0>
struct UnitPower {
	using type = Unit<power * UnitA::L_, power * UnitA::M_, power * UnitA::T_, power * UnitA::K_>;
};

template<typename UnitA, std::enable_if_t<isUnit<UnitA>::value, int> = 0>
struct UnitInverse {
	using type = typename UnitPower<UnitA, -1>::type;
};

template<typename UnitA, typename UnitB, std::enable_if_t<isUnit<UnitA>::value && isUnit<UnitB>::value, int> = 0>
struct UnitProduct {
	using type = Unit<UnitA::L_ + UnitB::L_, UnitA::M_ + UnitB::M_, UnitA::T_ + UnitB::T_, UnitA::K_ + UnitB::K_>;
};

template<typename UnitA, typename UnitB, std::enable_if_t<isUnit<UnitA>::value && isUnit<UnitB>::value, int> = 0>
struct UnitQuotient {
	using type = typename UnitProduct<UnitA, typename UnitInverse<UnitB>::type>::type;
};


struct ZeroQuantity {
	constexpr ZeroQuantity() = default;
	constexpr ZeroQuantity(double z) {};
};
static constexpr ZeroQuantity zeroQ{};

struct OneQuantity {
};
static constexpr OneQuantity oneQ{};


template<typename Units = Unit<0, 0, 0, 0>, typename Type = double, Type scale = 1.0>
struct Quantity {
	static constexpr Type a_ = scale;
	static constexpr Type ia_ = Type(1) / a_;

	template<Type otherScale>
	using ScaledOther = Quantity<Units, Type, otherScale>;


	constexpr Quantity(OneQuantity = oneQ) :
			value_(Type(1)) {
	}
	constexpr Quantity(ZeroQuantity) :
			value_(Type(0)) {
	}
	template<Type otherScale>
	constexpr Quantity(ScaledOther<otherScale> const &other) {
		*this = other * other.a_ * ia_;
	}
	template<Type otherScale>
	constexpr Quantity(ScaledOther<otherScale> &&other) {
		*this = std::move(other * other.a_ * ia_);
	}
	template<Type otherScale>
	constexpr Quantity& operator=(ScaledOther<otherScale> const &other) {
		value_ = other.value_ * other.a_ * ia_;
		return *this;
	}
	template<Type otherScale>
	constexpr Quantity& operator=(ScaledOther<otherScale> &&other) {
		value_ = std::move(other.value_ * other.a_ * ia_);
		return *this;
	}
	constexpr Quantity operator+() const {
		return *this;
	}
	constexpr Quantity operator-() const {
		Quantity negation;
		negation.value_ = -value_;
		return negation;
	}
	template<Type otherScale>
	constexpr Quantity operator+(ScaledOther<otherScale> const &other) const {
		Quantity sum;
		sum.value_ = value_ + other.value_ * other.a_ * ia_;
		return sum;
	}
	constexpr Quantity operator+(Type const &other) const {
		static_assert(std::is_same<Unit<0,0,0,0>, Units>::value, "Units mismatch");
		Quantity sum;
		sum.value_ = value_ + other * ia_;
		return sum;
	}
	constexpr Quantity operator-(Type const &other) const {
		static_assert(std::is_same<Unit<0,0,0,0>, Units>::value, "Units mismatch");
		Quantity sum;
		sum.value_ = value_ - other * ia_;
		return sum;
	}
	template<Type otherScale>
	constexpr Quantity operator-(ScaledOther<otherScale> const &other) const {
		Quantity difference;
		difference.value_ = value_ - other.value_ * other.a_ * ia_;
		return difference;
	}
	template<typename OtherUnits, Type otherScale>
	constexpr auto operator*(Quantity<OtherUnits, Type, otherScale> const &other) const {
		Quantity<typename UnitProduct<Units, OtherUnits>::type, Type, a_> product;
		product.value_ = value_ * other.value_ * other.a_ * ia_;
		return product;
	}
	constexpr auto operator*(Type const &other) const {
		Quantity product;
		product.value_ = value_ * other;
		return product;
	}
	template<typename OtherUnits, Type otherScale>
	constexpr auto operator/(Quantity<OtherUnits, Type, otherScale> const &other) const {
		Quantity<typename UnitQuotient<Units, OtherUnits>::type, Type, ia_> quotient;
		quotient.value_ = value_ * a_ *  other.ia_ / other.value_;
		return quotient;
	}
	constexpr auto operator/(Type const &other) const {
		Quantity quotient;
		quotient.value_ = value_  * a_ / other;
		return quotient;
	}
	template<Type otherScale>
	constexpr Quantity operator+=(ScaledOther<otherScale> const &other) {
		value_ += other.value_ * other.a_ * ia_;
		return *this;
	}
	template<Type otherScale>
	constexpr Quantity operator-=(ScaledOther<otherScale> const &other) {
		value_ -= other.value_ * other.a_ * ia_;
		return *this;
	}
	constexpr Quantity operator*=(Type const &other) {
		value_ *= other * ia_;
		return *this;
	}
	constexpr Quantity operator/=(Type const &other) {
		value_ /= other * ia_;
		return *this;
	}
	template<Type otherScale>
	constexpr Quantity operator*=(Quantity<Unit<0, 0, 0, 0>, Type, otherScale> const &other) {
		value_ *= other.value_ * other.a_ * ia_;
		return *this;
	}
	template<Type otherScale>
	constexpr Quantity operator/=(Quantity<Unit<0, 0, 0, 0>, Type, otherScale> const &other) {
		value_ /= other.value_ * other.a_ * ia_;
		return *this;
	}
	constexpr Quantity& operator=(Type value) {
		static_assert(std::is_same<Unit<0,0,0,0>, Units>::value, "Units mismatch");
		value_ = value;
		return *this;
	}
	constexpr explicit operator Type() const {
		static_assert(std::is_same<Unit<0, 0, 0, 0>, Units>::value, "Units mismatch");
		return value_;
	}
	template<Type otherScale>
	constexpr bool operator==(ScaledOther<otherScale> const &other) const {
		return value_ == other.value_ * other.a_ * ia_;
	}
	template<Type otherScale>
	constexpr bool operator!=(ScaledOther<otherScale> const &other) const {
		return value_ != other.value_ * other.a_ * ia_;
	}
	template<Type otherScale>
	constexpr bool operator>=(ScaledOther<otherScale> const &other) const {
		return value_ >= other.value_ * other.a_ * ia_;
	}
	template<Type otherScale>
	constexpr bool operator<=(ScaledOther<otherScale> const &other) const {
		return value_ <= other.value_ * other.a_ * ia_;
	}
	template<Type otherScale>
	constexpr bool operator>(ScaledOther<otherScale> const &other) const {
		return value_ > other.value_ * other.a_ * ia_;
	}
	template<Type otherScale>
	constexpr bool operator<(ScaledOther<otherScale> const &other) const {
		return value_ < other.value_ * other.a_ * ia_;
	}
	friend constexpr auto operator*(Type const &value, Quantity const &other) {
		Quantity product;
		product.value_ = value * other.value_ * other.a_ ;
		return product;
	}
	constexpr Type value() const noexcept {
		return value_;
	}
	friend constexpr auto operator/(Type const &value, Quantity const &other) {
		Quantity<typename UnitInverse<Units>::type, Type, ia_> quotient;
		quotient.value_ = value * other.ia_ / other.value_;
		return quotient;
	}
	template<Type otherScale>
	friend constexpr Quantity operator+(Type const &a, ScaledOther<otherScale> const &b) {
		return b + a;

	}
	template<Type otherScale>
	friend constexpr Quantity operator-(Type const &a, ScaledOther<otherScale> const &b) {
		return -(b - a);
	}
	template<typename, typename OtherType, OtherType>
	friend class Quantity;
private:
	Type value_;
};

template<typename U>
struct UnitSqrt;

template<int L, int M, int T, int K>
struct UnitSqrt<Unit<L, M, T, K>> {
    using type = Unit<L / 2, M / 2, T / 2, K / 2>;
};

template<typename Units, typename Type, Type Scale>
constexpr auto sqrt(Quantity<Units, Type, Scale> const &q) {
	Quantity<typename UnitSqrt<Units>::type, Type, sqrt(Scale)> units;
	return Quantity<typename UnitSqrt<Units>::type, Type, sqrt(Scale)>(std::sqrt(q.value()) * units);
}


namespace cgs {
using One = Quantity<Unit<0, 0, 0, 0>>;
using Centimeters = Quantity<Unit<1, 0, 0, 0>>;
using Grams = Quantity<Unit<0, 1, 0, 0>>;
using Seconds = Quantity<Unit<0, 0, 1, 0>>;
using Kelvin = Quantity<Unit<0, 0, 0, 1>>;
using Ergs = Quantity<Unit<2, 1, -2, 0>>;
using ErgsPerCm3 = Quantity<Unit<-1, 1, -2, 0>>;
static constexpr One one;
static constexpr Centimeters cm;
static constexpr Grams g;
static constexpr Seconds s;
static constexpr Kelvin K;
static constexpr auto s2 = s * s;
static constexpr auto s3 = s2 * s;
static constexpr auto cm2 = cm * cm;
static constexpr auto cm3 = cm2 * cm;
static constexpr auto K4 = K * K * K * K;
static constexpr auto dyne = (g * cm) / (s * s);
static constexpr auto erg = dyne * cm;
static constexpr auto hz = one / s;
}

#endif /* INCLUDE_UNITS_HPP_ */
