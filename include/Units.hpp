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

template<int L, int M, int T, int K, int N>
struct Unit {
	static constexpr int L_ = L;
	static constexpr int M_ = M;
	static constexpr int T_ = T;
	static constexpr int K_ = K;
	static constexpr int N_ = N;
};

template<typename Type>
struct isUnit {
	static constexpr bool value = false;
};

template<int L, int M, int T, int K, int N>
struct isUnit<Unit<L, M, T, K, N>> {
	static constexpr bool value = true;
};

template<typename UnitA, int power, std::enable_if_t<isUnit<UnitA>::value, int> = 0>
struct UnitPower {
	using type = Unit<power * UnitA::L_, power * UnitA::M_, power * UnitA::T_, power * UnitA::K_, power * UnitA::N_>;
};

template<typename UnitA, std::enable_if_t<isUnit<UnitA>::value, int> = 0>
struct UnitInverse {
	using type = typename UnitPower<UnitA, -1>::type;
};

template<typename UnitA, typename UnitB, std::enable_if_t<isUnit<UnitA>::value && isUnit<UnitB>::value, int> = 0>
struct UnitProduct {
	using type = Unit<UnitA::L_ + UnitB::L_, UnitA::M_ + UnitB::M_, UnitA::T_ + UnitB::T_, UnitA::K_ + UnitB::K_, UnitA::N_ + UnitB::N_>;
};

template<typename UnitA, typename UnitB, std::enable_if_t<isUnit<UnitA>::value && isUnit<UnitB>::value, int> = 0>
struct UnitQuotient {
	using type = typename UnitProduct<UnitA, typename UnitInverse<UnitB>::type>::type;
};

template<typename Units, typename Type>
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
		quotient.value_ = value  / other.value_;
		return quotient;
	}
	friend constexpr Quantity operator+(Type const &a, Quantity const &b) {
		return b + a;

	}
	friend constexpr Quantity operator-(Type const &a, Quantity const &b) {
		return -(b - a);
	}
	template<typename, typename>
	friend class Quantity;
private:
	Type value_;
};

template<typename U>
struct UnitSqrt;

template<int L, int M, int T, int K, int N>
struct UnitSqrt<Unit<L, M, T, K, N>> {
	using type = Unit<L / 2, M / 2, T / 2, K / 2, N / 2>;
};

template<typename U>
struct UnitCbrt;

template<int L, int M, int T, int K, int N>
struct UnitCbrt<Unit<L, M, T, K, N>> {
	using type = Unit<L / 3, M / 3, T / 3, K / 3, N / 3>;
};

template<typename Units, typename Type>
constexpr auto sqrt(Quantity<Units, Type> const &q) {
	using std::sqrt;
	using QuantityType = Quantity<typename UnitSqrt<Units>::type, Type>;
	QuantityType units(1.0);
	auto const sroot = sqrt(q.value());
	QuantityType const rc = sroot * units;
	return rc;
}

template<typename Units, typename Type>
constexpr auto cbrt(Quantity<Units, Type> const &q) {
	using std::cbrt;
	using QuantityType = Quantity<typename UnitCbrt<Units>::type, Type>;
	QuantityType units(1.0);
	auto const sroot = cbrt(q.value());
	QuantityType const rc = sroot * units;
	return rc;
}

template<typename Type>
constexpr auto log(Quantity<Unit<0, 0, 0, 0, 0>, Type> const &q) {
	using QuantityType = Quantity<Unit<0, 0, 0, 0, 0>, Type>;
	return QuantityType(log(q.value()));
}

template<typename Type>
constexpr auto pow(Quantity<Unit<0, 0, 0, 0, 0>, Type> const &q, double r) {
	using QuantityType = Quantity<Unit<0, 0, 0, 0, 0>, Type>;
	return QuantityType(pow(q.value(), r));
}

template<typename Type>
constexpr auto pow(Quantity<Unit<0, 0, 0, 0, 0>, Type> const &q, Quantity<Unit<0, 0, 0, 0, 0>, Type> r) {
	using QuantityType = Quantity<Unit<0, 0, 0, 0, 0>, Type>;
	return QuantityType(pow(q.value(), r.value()));
}

namespace cgs {
using One = Quantity<Unit<0, 0, 0, 0, 0>, double>;
using Centimeters = Quantity<Unit<1, 0, 0, 0, 0>, double>;
using Grams = Quantity<Unit<0, 1, 0, 0, 0>, double>;
using Seconds = Quantity<Unit<0, 0, 1, 0, 0>, double>;
using Kelvin = Quantity<Unit<0, 0, 0, 1, 0>, double>;
using Ergs = Quantity<Unit<2, 1, -2, 0, 0>, double>;
using ErgsPerCm3 = Quantity<Unit<-1, 1, -2, 0, 0>, double>;
using Mole = Quantity<Unit<0, 0, 0, 0, 1>, double>;
static constexpr One one(1.0);
static constexpr Centimeters cm(1.0);
static constexpr Grams g(1.0);
static constexpr Seconds s(1.0);
static constexpr Kelvin K(1.0);
static constexpr Mole mol(1.0);
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
