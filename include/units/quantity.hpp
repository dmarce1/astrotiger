/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_UNITS_QUANTITY_HPP_
#define INCLUDE_UNITS_QUANTITY_HPP_

#include "units/units.hpp"

template<typename, typename >
struct Quantity;

template<typename Units, typename Type = double>
struct Quantity {
	explicit constexpr Quantity(Type v = Type(0)) :
			value_(v) {
	}
	constexpr Quantity(Quantity const &other) {
		*this = other;
	}
	constexpr Quantity(Quantity &&other) {
		*this = std::move(other);
	}
	constexpr Quantity& operator=(Quantity const &other) {
		value_ = other.value_;
		return *this;
	}
	constexpr Quantity& operator=(Quantity &&other) {
		value_ = std::move(other.value_);
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
	constexpr Quantity operator+(Quantity const &other) const {
		Quantity sum;
		sum.value_ = value_ + other.value_;
		return sum;
	}
	constexpr Quantity operator+(Type const &other) const {
		static_assert(std::is_same<Unit<0, 0, 0, 0, 0>, Units>::value, "Units mismatch");
		Quantity sum;
		sum.value_ = value_ + other;
		return sum;
	}
	constexpr Quantity operator-(Type const &other) const {
		static_assert(std::is_same<Unit<0, 0, 0, 0, 0>, Units>::value, "Units mismatch");
		Quantity sum;
		sum.value_ = value_ - other;
		return sum;
	}
	constexpr Quantity operator-(Quantity const &other) const {
		Quantity difference;
		difference.value_ = value_ - other.value_;
		return difference;
	}
	template<typename OtherUnits>
	constexpr auto operator*(Quantity<OtherUnits, Type> const &other) const {
		Quantity<typename UnitProduct<Units, OtherUnits>::type, Type> product;
		product.value_ = value_ * other.value_;
		return product;
	}
	constexpr auto operator*(Type const &other) const {
		Quantity product;
		product.value_ = value_ * other;
		return product;
	}
	template<typename OtherUnits>
	constexpr auto operator/(Quantity<OtherUnits, Type> const &other) const {
		Quantity<typename UnitQuotient<Units, OtherUnits>::type, Type> quotient;
		quotient.value_ = value_ / other.value_;
		return quotient;
	}
	constexpr auto operator/(Type const &other) const {
		Quantity quotient;
		quotient.value_ = value_ / other;
		return quotient;
	}
	constexpr Quantity operator+=(Quantity const &other) {
		value_ += other.value_;
		return *this;
	}
	constexpr Quantity operator-=(Quantity const &other) {
		value_ -= other.value_;
		return *this;
	}
	constexpr Quantity operator*=(Type const &other) {
		value_ *= other;
		return *this;
	}
	constexpr Quantity operator/=(Type const &other) {
		value_ /= other;
		return *this;
	}
	constexpr Quantity operator*=(Quantity<Unit<0, 0, 0, 0, 0>, Type> const &other) {
		value_ *= other.value_;
		return *this;
	}
	constexpr Quantity operator/=(Quantity<Unit<0, 0, 0, 0, 0>, Type> const &other) {
		value_ /= other.value_;
		return *this;
	}
	constexpr Quantity& operator=(Type value) {
		static_assert(std::is_same<Unit<0, 0, 0, 0, 0>, Units>::value, "Units mismatch");
		value_ = value;
		return *this;
	}
	constexpr explicit operator Type() const {
		static_assert(std::is_same<Unit<0, 0, 0, 0, 0>, Units>::value, "Units mismatch");
		return value_;
	}
	constexpr bool operator==(Quantity const &other) const {
		return value_ == other.value_;
	}
	constexpr bool operator!=(Quantity const &other) const {
		return value_ != other.value_;
	}
	constexpr bool operator>=(Quantity const &other) const {
		return value_ >= other.value_;
	}
	constexpr bool operator<=(Quantity const &other) const {
		return value_ <= other.value_;
	}
	constexpr bool operator>(Quantity const &other) const {
		return value_ > other.value_;
	}
	constexpr bool operator<(Quantity const &other) const {
		return value_ < other.value_;
	}
	friend constexpr auto operator*(Type const &value, Quantity const &other) {
		Quantity product;
		product.value_ = value * other.value_;
		return product;
	}
	constexpr Type value() const noexcept {
		return value_;
	}
	friend constexpr auto operator/(Type const &value, Quantity const &other) {
		Quantity<typename UnitInverse<Units>::type, Type> quotient;
		quotient.value_ = value / other.value_;
		return quotient;
	}
	friend constexpr Quantity operator+(Type const &a, Quantity const &b) {
		return b + a;

	}
	friend constexpr Quantity operator-(Type const &a, Quantity const &b) {
		return -(b - a);
	}
	template<typename, typename >
	friend class Quantity;
private:
	Type value_;
};


template<Rational P, typename U, typename T>
constexpr auto pow(Quantity<U, T> const &q) {
	using UnitType = typename UnitPower<U, P>::type;
	return Quantity<UnitType, T>(pow(q.value(), P));
}

template<typename U, typename T>
constexpr auto sqrt(Quantity<U, T> const &q) {
	return pow<Rational(1, 2)>(q);
}

template<typename T>
constexpr auto log(Quantity<NullUnitType, T> const &q) {
	return log(q.value());
}

template<typename T>
constexpr auto exp(Quantity<NullUnitType, T> const &q) {
	return exp(q.value());
}

#endif /* INCLUDE_UNITS_QUANTITY_HPP_ */
