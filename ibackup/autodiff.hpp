/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "indices.hpp"
#include "quantity.hpp"
#include "util.hpp"
#include "vector.hpp"



template<typename Type, int order, typename IndependentVariableTypes = std::tuple<Type>>
struct FwdAutoDiff;

template<typename Type, int order, typename IndependentVariableTypes>
struct FwdAutoDiff {
	static constexpr int dimensionCount = std::tuple_size<IndependentVariableTypes>::value;
	using index_type = Indices<order, dimensionCount, int>;
	static constexpr size_t size() {
		return index_type::size();
	}
private:
	friend constexpr auto inverse(FwdAutoDiff const &x) {
		return compose([](Type x, size_t n) {
			using std::pow;
			return factorialPower(-1, n) / pow(x, n + 1);
		}, x);
	}
	template<class Function>
	friend constexpr auto compose(Function const &f, FwdAutoDiff const &g) {
		FwdAutoDiff h;
		Type const g0 = g[0];
		constexpr int Nmax = dimensionCount * order;
		std::array<Type, Nmax + 1> fn { };
		for (int n = 0; n <= Nmax; ++n) {
			fn[n] = f(g0, n) / factorial(n);
		}
		FwdAutoDiff δg = g - FwdAutoDiff(g0);
		h = FwdAutoDiff(fn[0]);
		FwdAutoDiff term(1.0);
		for (int n = 1; n <= Nmax; ++n) {
			term = term * δg;
			h += fn[n] * term;
		}
		return h;
	}
public:
	friend std::ostream& operator<<(std::ostream &os, FwdAutoDiff const &a) {
		os << "FwdAutoDiff<O=" << order << ", D=" << dimensionCount << ">{\n";
		index_type α = index_type::zero();
		while (α != index_type::end()) {
			os << "  α=(";
			for (int d = 0; d < dimensionCount; ++d) {
				os << α[d];
				if (d + 1 < dimensionCount)
					os << ",";
			}
			os << ") -> " << Type(a.C_[int(α)]) << "\n";
			α++;
		}
		os << "}";
		return os;
	}
	constexpr FwdAutoDiff() {
		C_.fill(0.0);
	}
	constexpr FwdAutoDiff(FwdAutoDiff const &other) {
		*this = other;
	}
	constexpr FwdAutoDiff(FwdAutoDiff &&other) {
		*this = std::move(other);
	}
	constexpr FwdAutoDiff(Type value) {
		*this = value;
	}
	constexpr FwdAutoDiff& operator=(FwdAutoDiff const &other) {
		C_ = other.C_;
		return *this;
	}
	constexpr FwdAutoDiff& operator=(FwdAutoDiff &&other) {
		C_ = std::move(other.C_);
		return *this;
	}
	constexpr FwdAutoDiff& operator=(Type value) {
		C_.fill(0.0);
		C_[0] = value;
		return *this;
	}
	constexpr Type operator()(std::array<Type, dimensionCount> const &x) const {
		Type y = 0.0;
		for (auto α = index_type::zero(); α != index_type::end(); α++) {
			y += C_[int(α)] * pow(x, α);
		}
		return y;
	}
	static constexpr FwdAutoDiff independent(Type value, size_t dim) {
		FwdAutoDiff diff;
		index_type indices = index_type::zero();
		indices[dim] = 1;
		diff.C_.fill(0.0);
		diff.C_[0] = value;
		diff.C_[int(indices)] = 1.0;
		return diff;
	}
	constexpr auto operator[](int a) const {
		return C_[a];
	}
	constexpr auto operator[](index_type α) const {
		return operator[](int(α));
	}
	constexpr auto operator+() const {
		return *this;
	}
	constexpr auto operator-() const {
		auto A = *this;
		for (size_t i = 0; i < size(); i++) {
			A.C_[i] = -this->C_[i];
		}
		return A;
	}
	constexpr auto operator+(FwdAutoDiff A) const {
		for (size_t i = 0; i < size(); i++) {
			A.C_[i] += this->C_[i];
		}
		return A;
	}
	constexpr auto operator-(FwdAutoDiff A) const {
		return *this + (-A);
	}
	constexpr auto operator*(FwdAutoDiff const &B) const {
		auto const &A = *this;
		FwdAutoDiff C;
		C.C_.fill(0.0);
		for (auto α = index_type::zero(); α != index_type::end(); α++) {
			for (auto β = index_type::zero(); β != index_type::end(); β++) {
				auto γ = α + β;
				if (γ.max() > order) {
					continue;
				}
				C.C_[int(γ)] += A.C_[int(α)] * B.C_[int(β)];
			}
		}
		return C;
	}
	constexpr auto operator/(FwdAutoDiff const &B) const {
		return operator*(inverse(B));
	}
	constexpr auto operator*(Type b) const {
		auto A = *this;
		for (size_t i = 0; i < size(); i++) {
			A.C_[i] *= b;
		}
		return A;
	}
	constexpr auto operator/(Type b) const {
		return operator*(1.0 / b);
	}
	constexpr FwdAutoDiff& operator*=(Type a) {
		*this = *this * a;
		return *this;
	}
	constexpr FwdAutoDiff& operator/=(Type a) {
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
	constexpr FwdAutoDiff& operator*=(FwdAutoDiff const &A) {
		*this = *this * A;
		return *this;
	}
	constexpr FwdAutoDiff& operator/=(FwdAutoDiff const &A) {
		*this = *this / A;
		return *this;
	}
	friend constexpr auto operator*(Type b, FwdAutoDiff A) {
		return A * b;
	}
	friend constexpr auto operator/(Type b, FwdAutoDiff A) {
		return b * inverse(A);
	}
	friend constexpr auto sqrt(FwdAutoDiff const &x) {
		constexpr Rational half(1, 2);
		return compose([](Type x0, size_t n) {
			using std::pow;
			using std::sqrt;
			if (n == 0) {
				return sqrt(x0);
			}
			return sqrt(x0) * factorialPower(half, n) / pow(x0, n);
		}, x);
	}

	friend constexpr auto cbrt(FwdAutoDiff const &x) {
		constexpr Rational third(1, 3);
		return compose([](Type x0, size_t n) {
			using std::pow;
			using std::cbrt;
			if (n == 0) {
				return cbrt(x0);
			}
			return cbrt(x0) * factorialPower(third, n) / pow(x0, n);
		}, x);
	}
	template<Rational power>
	friend constexpr auto pow(FwdAutoDiff const &x) {
		if constexpr (power.denominator() == 1) {
			if constexpr (power.numerator() < 0) {
				constexpr int m = -power.numerator();
				return inverse(pow<Rational(m)>(x));
			} else {
				return compose([](Type x, size_t n) {
					return pow<power.numerator()>(x) * factorialPower(power, n) / pow(x, n);
				}, x);
			}
		} else if constexpr (power.denominator() == 2) {
			return ::pow<power.numerator()>(x) * sqrt(x);
		} else if constexpr (power.denominator() == 3) {
			return ::pow<power.numerator()>(x) * cbrt(x);
		} else {
			return compose([](Type x, size_t n) {
				using std::pow;
				return pow(x, power - Rational(n)) * factorialPower(power, n);
			}, x);
		}
	}
	friend constexpr auto exp(FwdAutoDiff const &x) {
		using std::exp;
		auto const f = compose([](Type x, size_t) {
			return exp(x);
		}, x);
		return f;
	}
	friend constexpr auto log(FwdAutoDiff const &x) {
		using std::log;
		return compose([](Type x, size_t n) {
			if (n == 0) {
				return log(x);
			} else {
				return nonepow(n - 1) * factorial(n - 1) / std::pow(x, n);
			}
		}, x);
	}
	constexpr explicit operator Type() const {
		return C_[0];
	}
	auto constexpr gradient() const {
		std::array<FwdAutoDiff, dimensionCount> grad;
		for (int k = 0; k < dimensionCount; k++) {
			grad[k].C_.fill(Type(0));
			for (auto α = index_type::zero(); α != index_type::end(); α++) {
				if (α[k]) {
					const auto β = α - index_type::unit(k);
					grad[k].C_[β] = C_[α];
				}
			}
		}
		return grad;
	}
	template<typename, int, typename >
	friend struct FwdAutoDiff;
private:
	std::array<Type, size()> C_;
};




template<typename DependentUnitType, int order, typename IndependentVariableTypes>
struct FwdAutoDiff<Quantity<DependentUnitType, typename std::tuple_element_t<0, IndependentVariableTypes>::value_type>, order, IndependentVariableTypes> {
	static constexpr size_t dimensionCount = std::tuple_size<IndependentVariableTypes>::value;
	using Real = typename std::tuple_element_t<0, IndependentVariableTypes>::value_type;
	using index_type = typename FwdAutoDiff<Real, order, std::array<Real, dimensionCount>>::index_type;
private:
	using DependentVariableType = Quantity<DependentUnitType, Real>;
	using IndependentUnitTypes = typename UnwrapTuple<IndependentVariableTypes>::type;
	static constexpr auto nextLower(index_type α) {
		size_t dim = -1;
		for (size_t d = 0; d < dimensionCount; d++) {
			if (α[d]) {
				dim = d;
				α[d]--;
				break;
			}
		}
		return std::pair(dim, α);
	}
	template<size_t ... I>
	static auto makeGradientType(std::index_sequence<I...>) {
		return std::tuple<typename UnitQuotient<DependentUnitType, typename std::tuple_element<I, IndependentUnitTypes>::type>::type...> { };
	}
	template<std::array<int, dimensionCount> α>
	constexpr auto typeHelper() const {
		constexpr auto rc = nextLower(α);
		constexpr auto β = rc.second;
		if constexpr (α != β) {
			constexpr auto dim = rc.first;
			using ThisType = typename std::tuple_element<dim, IndependentUnitTypes>::type;
			using RestType = decltype(typeHelper<β>());
			using NextType = typename UnitProduct<ThisType, RestType>::type;
			constexpr NextType rc { };
			return rc;
		} else {
			constexpr NullUnitType rc { };
			return rc;
		}
	}
	template<size_t I>
	auto makeGradient() const {
		if constexpr (I < dimensionCount) {
			constexpr auto θ = index_type::unit(I);
			using ThisType = Quantity<typename UnitQuotient<DependentUnitType, typename std::tuple_element<I, IndependentUnitTypes>::type>::type, Real>;
			std::tuple<ThisType> thisTuple;
			std::get<0>(thisTuple) = ThisType(autodiff_[θ]);
			auto const returnTuple = std::tuple_cat(thisTuple, makeGradient<I + 1>());
			if constexpr (AllSame<IndependentUnitTypes>::value && (I == 0)) {
				return std::apply([](auto const &... values) {
					return std::array<ThisType, dimensionCount> { { values... } };
				}, returnTuple);
			} else {
				return returnTuple;
			}
		} else {
			return std::tuple<> { };
		}
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
	constexpr FwdAutoDiff(DependentVariableType const &qty) :
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
	constexpr FwdAutoDiff& operator=(DependentVariableType const &qty) {
		autodiff_ = qty.value();
		return *this;
	}
	template<size_t dim = 0>
	static constexpr auto independent(DependentVariableType const &qty) {
		FwdAutoDiff<DependentVariableType, order, IndependentVariableTypes> var;
		var.autodiff_ = FwdAutoDiff<Real, order, std::array<Real, dimensionCount>>::independent(qty.value(), dim);
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
	FwdAutoDiff& operator*=(DependentVariableType const &a) {
		static_assert(std::is_same<NullUnitType, DependentUnitType>::value, "Units mismatch");
		*this = *this * a;
		return *this;
	}
	template<typename OtherDependentUnitType>
	auto operator*(FwdAutoDiff<Quantity<OtherDependentUnitType, Real>, order, IndependentVariableTypes> const &A) const {
		using ReturnDependentType = Quantity<typename UnitProduct<DependentUnitType, OtherDependentUnitType>::type, Real>;
		FwdAutoDiff<ReturnDependentType, order, IndependentVariableTypes> B;
		B.autodiff_ = autodiff_ * A.autodiff_;
		return B;
	}
	template<typename OtherDependentUnitType>
	auto operator*(Quantity<OtherDependentUnitType, Real> const &a) const {
		using ReturnDependentType = Quantity<typename UnitProduct<DependentUnitType, OtherDependentUnitType>::type, Real>;
		FwdAutoDiff<ReturnDependentType, order, IndependentVariableTypes> B;
		B.autodiff_ = autodiff_ * a.value();
		return B;
	}
	auto operator*(Real a) const {
		FwdAutoDiff<DependentVariableType, order, IndependentVariableTypes> B;
		B.autodiff_ = autodiff_ * a;
		return B;
	}
	template<typename OtherDependentType>
	friend auto operator*(Quantity<OtherDependentType, Real> b, FwdAutoDiff A) {
		return A * b;
	}
	FwdAutoDiff& operator/=(FwdAutoDiff const &A) {
		*this = *this / A;
		return *this;
	}
	FwdAutoDiff& operator/=(DependentVariableType const &a) {
		static_assert(std::is_same<NullUnitType, DependentUnitType>::value, "Units mismatch");
		*this = *this / a;
		return *this;
	}
	template<typename OtherDependentUnitType>
	auto operator/(FwdAutoDiff<Quantity<OtherDependentUnitType, Real>, order, IndependentVariableTypes> const &A) const {
		using ReturnDependentType = Quantity<typename UnitQuotient<DependentUnitType, OtherDependentUnitType>::type, Real>;
		FwdAutoDiff<ReturnDependentType, order, IndependentVariableTypes> B;
		B.autodiff_ = autodiff_ / A.autodiff_;
		return B;
	}
	template<typename OtherDependentUnitType>
	auto operator/(Quantity<OtherDependentUnitType, Real> const &a) const {
		using ReturnDependentType = Quantity<typename UnitQuotient<DependentUnitType, OtherDependentUnitType>::type, Real>;
		FwdAutoDiff<ReturnDependentType, order, IndependentVariableTypes> B;
		B.autodiff_ = autodiff_ / a.value();
		return B;
	}
	auto operator/(Real a) const {
		FwdAutoDiff<DependentVariableType, order, IndependentVariableTypes> B;
		B.autodiff_ = autodiff_ / a;
		return B;
	}
	template<typename OtherDependentType>
	friend auto operator/(Quantity<OtherDependentType, Real> b, FwdAutoDiff A) {
		return b.value() / A;
	}
	friend auto operator/(Real b, FwdAutoDiff A) {
		return A / b;
	}
	friend auto exp(FwdAutoDiff x) {
		static_assert(std::is_same<NullUnitType, DependentUnitType>::value, "Units mismatch");
		x.autodiff_ = exp(x.autodiff_);
		return x;
	}
	friend auto log(FwdAutoDiff x) {
		static_assert(std::is_same<NullUnitType, DependentUnitType>::value, "Units mismatch");
		x.autodiff_ = log(x.autodiff_);
		return x;
	}
	template<Rational power>
	friend auto pow(FwdAutoDiff const &x) {
		using OtherDependentType = Quantity<typename UnitPower<DependentUnitType, power>::type, Real>;
		FwdAutoDiff<OtherDependentType, order, IndependentVariableTypes> B;
		B.autodiff_ = pow<power>(x.autodiff_);
		return B;
	}
	friend auto sqrt(FwdAutoDiff const &x) {
		using OtherDependentType = Quantity<typename UnitPower<DependentUnitType, Rational(1, 2)>::type, Real>;
		FwdAutoDiff<OtherDependentType, order, IndependentVariableTypes> B;
		B.autodiff_ = sqrt(x.autodiff_);
		return B;
	}
	friend auto cbrt(FwdAutoDiff const &x) {
		using OtherDependentType = typename UnitPower<DependentVariableType, Rational(1, 3)>::type;
		FwdAutoDiff<OtherDependentType, order, IndependentVariableTypes> B;
		B.autodiff_ = cbrt(x.autodiff_);
		return B;
	}
	template<std::array<int, dimensionCount> a>
	auto multiGet() const {
		using RUnitType = typename UnitQuotient<DependentUnitType, decltype(typeHelper<a>())>::type;
		return Quantity<RUnitType, Real>(autodiff_[a]);
	}
	template<auto a>
	auto get() const {
		constexpr index_type α(std::array<int, dimensionCount>( { a }));
		using IndependentUnitType = typename std::tuple_element<0, IndependentUnitTypes>::type;
		using ReturnType = typename UnitQuotient<DependentUnitType, IndependentUnitType>::type;
		return Quantity<ReturnType, Real>(autodiff_[α]);
	}
	auto gradient() const {
		std::array<FwdAutoDiff, dimensionCount> grad;
		for (int k = 0; k < dimensionCount; k++) {
			grad[k].autodiff_ = DependentVariableType(0);
			for (auto α = index_type::zero(); α != index_type::end(); α++) {
				if (α[k]) {
					const auto β = α - index_type::unit(k);
					grad[k].autodiff_[β] = autodiff_[α];
				}
			}
		}
		return grad;
	}
	static constexpr size_t size() {
		return pow<order>(dimensionCount);
	}
	constexpr operator DependentVariableType() const {
		return DependentVariableType(Real(autodiff_));
	}
	template<typename, int, typename >
	friend struct FwdAutoDiff;
private:
	FwdAutoDiff<Real, order, std::array<Real, dimensionCount>> autodiff_;
};

