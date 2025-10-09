#pragma once

#include <array>
#include <cmath>
#include <type_traits>

#include "numbers/rational.hpp"

template<int O, int D>
class FwdAutoDiffDouble {
	using index_type = std::array<size_t, D>;
	template<int n>
	static constexpr auto pow(double x) {
		if constexpr (n >= 0) {
			double xm = x;
			double xn = 1.0;
			int m = n;
			while (m) {
				if (m & 1) {
					xn *= xm;
				}
				m >>= 1;
				xm *= xm;
			}
			return xn;
		} else {
			return 1.0 / pow<-n>(x);
		}
	}
	static constexpr size_t size() {
		return pow<O>(D);
	}
	static constexpr double nonepow(int n) {
		using std::abs;
		return (abs(n) & 1) ? -1.0 : +1.0;
	}
	static constexpr double factorial(size_t n) {
		if (n > 0) {
			return n * factorial(n - 1);
		} else {
			return 1.0;
		}
	}
	static constexpr double factorialPower(Rational x, size_t n) {
		if (n) {
			n--;
			return (x - Rational(n)) * factorialPower(x, n);
		} else {
			return 1.0;
		}
	}
	static constexpr bool lte(index_type const &α, size_t b) {
		for (int d = 0; d < D; d++) {
			if (α[d] > b) {
				return false;
			}
		}
		return true;
	}
	static constexpr bool lt(index_type const &α, size_t b) {
		return !gte(α, b);
	}
	static constexpr bool gt(index_type const &α, size_t b) {
		return !lte(α, b);
	}
	static constexpr bool gte(index_type const &α, size_t b) {
		for (int d = 0; d < D; d++) {
			if (α[d] < b) {
				return false;
			}
		}
		return true;
	}
	static constexpr bool lte(index_type const &α, index_type const &β) {
		for (int d = 0; d < D; d++) {
			if (α[d] > β[d]) {
				return false;
			}
		}
		return true;
	}
	static constexpr bool lt(index_type const &α, index_type const &β) {
		return !gte(α, β);
	}
	static constexpr bool gt(index_type const &α, index_type const &β) {
		return !lte(α, β);
	}
	static constexpr bool gte(index_type const &α, index_type const &β) {
		return lte(β, α);
	}
	static constexpr bool eq(index_type const &α, index_type const &β) {
		for (int d = 0; d < D; d++) {
			if (α[d] != β[d]) {
				return false;
			}
		}
		return true;
	}
	static constexpr bool neq(index_type const &α, index_type const &β) {
		return !eq(α, β);
	}
	static constexpr size_t sum(index_type const &α) {
		size_t sum = 0;
		for (size_t d = 0; d < D; d++) {
			sum += α[d];
		}
		return sum;
	}
	static constexpr double factorial(index_type const &α) {
		double f = 1.0;
		for (size_t d = 0; d < D; d++) {
			f *= factorial(α[d]);
		}
		return f;
	}
	static constexpr size_t flatten(index_type const &α) {
		size_t a = 0;
		for (size_t i = 0; i < D; i++) {
			a = O * a + α[i];
		}
		return a;
	}
	static constexpr bool increment(index_type &α, index_type const &β) {
		int dim = D - 1;
		while (++α[dim] + β[dim] == O) {
			if (dim == 0) {
				return false;
			}
			α[dim--] = 0;
		}
		return true;
	}
	static constexpr index_type zero() {
		index_type ζ;
		ζ.fill(0);
		return ζ;
	}
	static constexpr index_type add(index_type const &α, index_type const &β) {
		index_type γ;
		for (size_t i = 0; i < D; i++) {
			γ[i] = α[i] + β[i];
		}
		return γ;
	}
	static constexpr index_type sub(index_type const &α, index_type const &β) {
		index_type γ;
		for (size_t i = 0; i < D; i++) {
			γ[i] = α[i] - β[i];
		}
		return γ;
	}
	static constexpr auto pow(std::array<double, D> const &x, index_type const &α) {
		double z = 1.0;
		for (size_t d = 0; d < D; d++) {
			for (size_t n = 0; n < α[d]; n++) {
				z *= x[d];
			}
		}
		return z;
	}
	double operator()(std::array<double, D> const &x) const {
		double y = 0.0;
		for (auto α = zero();; increment(α)) {
			y += C_[flatten(α)] * pow(x, α);
		}
		return y;
	}
	template<class Function>
	friend auto compose(Function const &foo, FwdAutoDiffDouble const &g) {
		FwdAutoDiffDouble h;
		std::array<double, O> f;
		auto const x = g.C_[0];
		for (size_t a = 0; a < O; a++) {
			f[a] = foo(x, a);
		}
		for (auto β = zero();; increment(β)) {
			auto βn = zero();
			for (size_t a = 0; gte(βn, a); a++) {
				h.C_[flatten(βn)] += f[a] * g.C_[flatten(β)] / (factorial(a) * factorial(β));
				βn = add(βn, β);
			}
		}
		return h;
	}
	constexpr auto inverse(FwdAutoDiffDouble const &x) const {
		return compose([](double x, size_t n) {
			using std::pow;
			return factorialPower(-1, n) / (pow(x, n + 1) * factorial(n));
		}, x);
	}
	static constexpr bool increment(index_type &α) {
		constexpr auto ζ = zero();
		return increment(α, ζ);
	}
public:
	constexpr FwdAutoDiffDouble() {
		C_.fill(0.0);
	}
	constexpr FwdAutoDiffDouble(FwdAutoDiffDouble const &other) {
		*this = other;
	}
	constexpr FwdAutoDiffDouble(FwdAutoDiffDouble &&other) {
		*this = std::move(other);
	}
	constexpr FwdAutoDiffDouble(double value) {
		*this = value;
	}
	constexpr FwdAutoDiffDouble& operator=(FwdAutoDiffDouble const &other) {
		C_ = other.C_;
		return *this;
	}
	constexpr FwdAutoDiffDouble& operator=(FwdAutoDiffDouble &&other) {
		C_ = std::move(other.C_);
		return *this;
	}
	constexpr FwdAutoDiffDouble& operator=(double value) {
		C_.fill(0.0);
		C_[0] = value;
		return *this;
	}
	static constexpr FwdAutoDiffDouble independent(double value, size_t dim) {
		FwdAutoDiffDouble diff;
		index_type indices;
		indices.fill(0);
		indices[dim] = 1;
		diff.C_.fill(0.0);
		diff.C_[0] = value;
		diff.C_[flatten(indices)] = 1.0;
		return diff;
	}
	auto operator[](auto α) const {
		if constexpr (D == 1) {
			return C_[α];
		} else {
			return C_[flatten(α)];
		}
	}
	auto operator+() const {
		return *this;
	}
	auto operator-() const {
		auto A = *this;
		for (size_t i = 0; i < size(); i++) {
			A.C_[i] = -this->C_[i];
		}
		return A;
	}
	auto operator+(FwdAutoDiffDouble A) const {
		for (size_t i = 0; i < size(); i++) {
			A.C_[i] += this->C_[i];
		}
		return A;
	}
	auto operator-(FwdAutoDiffDouble A) const {
		return *this + (-A);
	}
	auto operator*(FwdAutoDiffDouble const &B) const {
		auto const &A = *this;
		FwdAutoDiffDouble C;
		C.C_.fill(0.0);
		for (auto α = zero();; increment(α)) {
			for (auto β = zero();; increment(β, α)) {
				auto const γ = add(α, β);
				C.C_[flatten(γ)] += A.C_[flatten(α)] * B.C_[flatten(β)];
			}
		}
		return C;
	}
	auto operator/(FwdAutoDiffDouble const &B) const {
		return operator*(inverse(B));
	}
	auto operator*(double b) const {
		auto A = *this;
		for (size_t i = 0; i < size(); i++) {
			A.C_[i] *= b;
		}
		return A;
	}
	auto operator/(double b) const {
		return operator*(1.0 / b);
	}
	FwdAutoDiffDouble& operator*=(double a) {
		*this = *this * a;
		return *this;
	}
	FwdAutoDiffDouble& operator/=(double a) {
		*this = *this / a;
		return *this;
	}
	FwdAutoDiffDouble& operator+=(FwdAutoDiffDouble const &A) {
		*this = *this + A;
		return *this;
	}
	FwdAutoDiffDouble& operator-=(FwdAutoDiffDouble const &A) {
		*this = *this - A;
		return *this;
	}
	FwdAutoDiffDouble& operator*=(FwdAutoDiffDouble const &A) {
		*this = *this * A;
		return *this;
	}
	FwdAutoDiffDouble& operator/=(FwdAutoDiffDouble const &A) {
		*this = *this / A;
		return *this;
	}
	friend auto operator*(double b, FwdAutoDiffDouble A) {
		return A * b;
	}
	friend auto operator/(double b, FwdAutoDiffDouble A) {
		return b * inverse(A);
	}
	friend auto exp(FwdAutoDiffDouble const &x) {
		using std::exp;
		return compose([](double x, size_t) {
			return exp(x);
		}, x);
	}
	friend auto log(FwdAutoDiffDouble const &x) {
		using std::log;
		return compose([](double x, size_t n) {
			if (n) {
				return nonepow(n - 1) / (factorial(n) * factorial(n - 1) * pow<n>(x));
			} else {
				return log(x);
			}
		}, x);
	}
	friend auto sqrt(FwdAutoDiffDouble const &x) {
		using std::pow;
		using std::sqrt;
		constexpr Rational half(1, 2);
		double const sqrtx = sqrt(x.C_[0]);
		return compose([sqrtx](double x, size_t n) {
			return sqrtx * factorialPower(half, n) / (pow(x, n) * factorial(n));
		}, x);
	}
	friend auto cbrt(FwdAutoDiffDouble const &x) {
		using std::pow;
		using std::cbrt;
		constexpr Rational third(1, 3);
		double const cbrtx = cbrt(x.C_[0]);
		return compose([cbrtx](double x, size_t n) {
			return cbrtx * factorialPower(third, n) / (pow(x, n) * factorial(n));
		}, x);
	}
	template<Rational power>
	friend auto pow(FwdAutoDiffDouble const &x) {
		if constexpr (power.denominator() == 1) {
			return compose([](double x, size_t n) {
				return pow<power.numerator()>(x) * factorialPower(power, n) / (pow(x, n) * factorial(n));
			}, x);
		} else if constexpr (power.denominator() == 2) {
			return pow<power.numerator()>(x) * sqrt(x);
		} else if constexpr (power.denominator() == 3) {
			return pow<power.numerator()>(x) * cbrt(x);
		} else {
			using std::pow;
			return compose([](double x, size_t n) {
				return pow(x, power - Rational(n)) * factorialPower(power, n) / factorial(n);
			}, x);
		}
	}
	constexpr operator double() const {
		return C_[0];
	}
	template<typename, int, typename>
	friend struct FwdAutoDiff;
private:
	std::array<double, size()> C_;
};
