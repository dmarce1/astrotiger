#pragma once

#include "units.hpp"

template<typename Units, typename Type = double>
struct Quantity {
	using units_type = Units;
	explicit constexpr Quantity() = default;
	constexpr Quantity(Type v) :
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

template<typename T>
struct IsQuantity {
	static constexpr bool value = false;
};

template<typename T1, typename T2>
struct IsQuantity<Quantity<T1, T2>> {
	static constexpr bool value = true;
};
