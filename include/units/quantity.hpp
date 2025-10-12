#pragma once

#include "units/units.hpp"

template<typename Units, typename Type = double>
struct Quantity {
	using units_type = Units;
	explicit constexpr Quantity() = default;
	explicit constexpr Quantity(Type v) :
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

#include "numbers/autodiff.hpp"

template<typename, int, typename >
struct FwdAutoDiff;

template<typename U1, int O, typename I>
class FwdAutoDiff {
	static constexpr size_t D = std::tuple_size<I>::value;
	using index_type = typename FwdAutoDiffDouble<O, D>::index_type;
	using value_type = Quantity<U1, double>;
	static constexpr auto nextLower(index_type α) {
		size_t dim = -1;
		for (size_t d = 0; d < D; d++) {
			if (α[d]) {
				dim = d;
				α[d]--;
				break;
			}
		}
		return std::pair(dim, α);
	}
public:
	constexpr FwdAutoDiff() :
			autodiff_() {
	}
	constexpr FwdAutoDiff(FwdAutoDiff const &other) :
			autodiff_(other.autodiff_) {
	}
	constexpr FwdAutoDiff(FwdAutoDiff &&other) :
			autodiff_(std::move(other.autodiff_)) {
	}
	constexpr FwdAutoDiff(value_type const &qty) :
			autodiff_(qty.value()) {
	}
	constexpr FwdAutoDiff& operator=(FwdAutoDiff const &other) {
		autodiff_ = other.autodiff_;
		return *this;
	}
	constexpr FwdAutoDiff& operator=(FwdAutoDiff &&other) {
		autodiff_ = std::move(other.autodiff_);
		return *this;
	}
	constexpr FwdAutoDiff& operator=(value_type const &qty) {
		autodiff_ = qty.value();
		return *this;
	}
	template<size_t dim = 0>
	static constexpr auto independent(value_type const &qty) {
		FwdAutoDiff<U1, O, I> var;
		var.autodiff_ = FwdAutoDiffDouble<O, D>::independent(qty.value(), dim);
		return var;
	}
	auto operator+() const {
		return *this;
	}
	auto operator-() const {
		FwdAutoDiff A;
		A.autodiff_ = -autodiff_;
		return A;
	}
	FwdAutoDiff& operator+=(FwdAutoDiff const &A) {
		*this = *this + A;
		return *this;
	}
	auto operator+(FwdAutoDiff A) const {
		FwdAutoDiff B;
		B.autodiff_ = autodiff_ + A.autodiff_;
		return B;
	}
	FwdAutoDiff& operator-=(FwdAutoDiff const &A) {
		*this = *this - A;
		return *this;
	}
	auto operator-(FwdAutoDiff A) const {
		FwdAutoDiff B;
		B.autodiff_ = autodiff_ - A.autodiff_;
		return B;
	}
	FwdAutoDiff& operator*=(FwdAutoDiff const &A) {
		*this = *this * A;
		return *this;
	}
	FwdAutoDiff& operator*=(value_type const &a) {
		static_assert(std::is_same<NullUnitType, typename value_type::units_type>::value, "Units mismatch");
		*this = *this * a;
		return *this;
	}
	template<typename U2>
	auto operator*(FwdAutoDiff<U2, O, I> const &A) const {
		using U3 = typename UnitProduct<U1, U2>::type;
		FwdAutoDiff<U3, O, I> B;
		B.autodiff_ = autodiff_ * A.autodiff_;
		return B;
	}
	template<typename U2>
	auto operator*(Quantity<U2, double> const &a) const {
		using U3 = typename UnitProduct<U1, U2>::type;
		FwdAutoDiff<U3, O, I> B;
		B.autodiff_ = autodiff_ * a.value();
		return B;
	}
	auto operator*(double a) const {
		FwdAutoDiff<U1, O, I> B;
		B.autodiff_ = autodiff_ * a;
		return B;
	}
	template<typename U2>
	friend auto operator*(Quantity<U2, double> b, FwdAutoDiff A) {
		return A * b;
	}
	FwdAutoDiff& operator/=(FwdAutoDiff const &A) {
		*this = *this / A;
		return *this;
	}
	FwdAutoDiff& operator/=(value_type const &a) {
		static_assert(std::is_same<NullUnitType, typename value_type::units_types>::value, "Units mismatch");
		*this = *this / a;
		return *this;
	}
	template<typename U2>
	auto operator/(FwdAutoDiff<U2, O, I> const &A) const {
		using U3 = typename UnitQuotient<U1, U2>::type;
		FwdAutoDiff<U3, O, I> B;
		B.autodiff_ = autodiff_ / A.autodiff_;
		return B;
	}
	template<typename U2>
	auto operator/(Quantity<U2, double> const &a) const {
		using U3 = typename UnitQuotient<U1, U2>::type;
		FwdAutoDiff<U3, O, I> B;
		B.autodiff_ = autodiff_ / a.value();
		return B;
	}
	auto operator/(double a) const {
		FwdAutoDiff<U1, O, I> B;
		B.autodiff_ = autodiff_ / a;
		return B;
	}
	template<typename U2>
	friend auto operator/(Quantity<U2, double> b, FwdAutoDiff A) {
		return b.value() / A;
	}
	friend auto operator/(double b, FwdAutoDiff A) {
		return A / b;
	}
	friend auto exp(FwdAutoDiff x) {
		static_assert(std::is_same<NullUnitType, typename value_type::units_type>::value, "Units mismatch");
		x.autodiff_ = exp(x.autodiff_);
		return x;
	}
	friend auto log(FwdAutoDiff x) {
		static_assert(std::is_same<NullUnitType, typename value_type::units_type>::value, "Units mismatch");
		x.autodiff_ = log(x.autodiff_);
		return x;
	}
	template<Rational power>
	friend auto pow(FwdAutoDiff const &x) {
		using U2 = typename UnitPower<U1, power>::type;
		FwdAutoDiff<U2, O, I> B;
		B.autodiff_ = pow<power>(x.autodiff_);
		return B;
	}
	friend auto sqrt(FwdAutoDiff const &x) {
		using U2 = typename UnitPower<U1, Rational(1, 2)>::type;
		FwdAutoDiff<U2, O, I> B;
		B.autodiff_ = sqrt(x.autodiff_);
		return B;
	}
	friend auto cbrt(FwdAutoDiff const &x) {
		using U2 = typename UnitPower<U1, Rational(1, 3)>::type;
		FwdAutoDiff<U2, O, I> B;
		B.autodiff_ = cbrt(x.autodiff_);
		return B;
	}
	template<index_type α>
	auto get() const {
		using R = typename UnitQuotient<U1, decltype(typeHelper<α>())>::type;
		return Quantity<R, double>(autodiff_[α]);
	}
	template<index_type α>
	constexpr auto typeHelper() const {
		constexpr auto rc = nextLower(α);
		constexpr auto β = rc.second;
		if constexpr (α != β) {
			constexpr auto dim = rc.first;
			using ThisType = typename std::tuple_element<dim, I>::type;
			using RestType = decltype(typeHelper<β>());
			using NextType = typename UnitProduct<typename ThisType::units_type, RestType>::type;
			constexpr NextType rc { };
			return rc;
		} else {
			constexpr NullUnitType rc { };
			return rc;
		}
	}
	template<auto γ>
	auto get() const {
		if constexpr (std::is_same<decltype(γ), index_type>::value) {
			using R = typename UnitQuotient<U1, decltype(typeHelper<γ>())>::type;
			return Quantity<R, double>(autodiff_[γ]);
		} else {
			constexpr index_type α { γ };
			using R = typename UnitQuotient<U1, decltype(typeHelper<α>())>::type;
			return Quantity<R, double>(autodiff_[α]);
		}
	}
	static constexpr size_t size() {
		return pow<O>(D);
	}
	constexpr operator value_type() const {
		return value_type(double(autodiff_));
	}
	template<typename, int, typename >
	friend struct FwdAutoDiff;
private:
	FwdAutoDiffDouble<O, D> autodiff_;
};

