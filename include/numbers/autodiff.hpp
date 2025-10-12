#pragma once

#include <array>
#include <cmath>
#include <type_traits>
#include <unordered_map>

#include "BellPolynomial.hpp"
#include "numbers/rational.hpp"

template<int O, int D>
class FwdAutoDiffDouble {
	using index_type = std::array<int, D>;
	struct index_typeHash {
		size_t operator()(index_type const &α) const {
			constexpr size_t n = ((size_t(1) << size_t(31)) - size_t(1));
			constexpr size_t a = 48271;
			constexpr size_t b = 16807;
			std::hash<int> keyGen;
			size_t key = 42;
			for (int i = 0; i < D; i++) {
				key = ((a * key) % n) ^ ((b * keyGen(i)) % n);
			}
			return key;
		}
	};
	static constexpr auto pow(double x, int n) {
		if (n >= 0) {
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
			return 1.0 / pow(x, -n);
		}
	}
	template<int n>
	static constexpr auto pow(double x) {
		return pow(x, n);
	}
	static constexpr size_t size() {
		return pow<D>(O);
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
	static constexpr bool gte(index_type const &α, int b) {
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
	static constexpr auto sum(index_type const &α) {
		int sum = 0;
		for (size_t d = 0; d < D; d++) {
			sum += α[d];
		}
		return sum;
	}
	static constexpr auto max(index_type const &α) {
		auto max = std::numeric_limits<int>::min();
		for (size_t d = 0; d < D; d++) {
			auto const a = α[d];
			if (a > max) {
				max = a;
			}
		}
		return max;
	}
	static constexpr double factorial(index_type const &α) {
		double f = 1.0;
		for (size_t d = 0; d < D; d++) {
			f *= factorial(α[d]);
		}
		return f;
	}
	static constexpr index_type zero() {
		index_type ζ;
		ζ.fill(0);
		return ζ;
	}
	static constexpr index_type unit(size_t i) {
		index_type ζ = zero();
		ζ[i] = 1;
		return ζ;
	}
	static constexpr index_type end() {
	    index_type ζ;
	    ζ.fill(O);                      // sentinel beyond last valid
	    return ζ;
	}
	static constexpr index_type one() {
		index_type ζ;
		ζ.fill(+1);
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
	static constexpr index_type mul(index_type const &α, size_t b) {
		index_type γ;
		for (size_t i = 0; i < D; i++) {
			γ[i] = α[i] * b;
		}
		return γ;
	}
	static constexpr auto pow(std::array<double, D> const &x, index_type const &α) {
		double z = 1.0;
		for (size_t d = 0; d < D; d++) {
			for (int n = 0; n < α[d]; n++) {
				z *= x[d];
			}
		}
		return z;
	}
	friend constexpr auto inverse(FwdAutoDiffDouble const &x) {
		return compose([](double x, size_t n) {
			using std::pow;
			return factorialPower(-1, n) / pow(x, n + 1);
		}, x);
	}
	static constexpr void increment(index_type& α) {
	    for (int d = D - 1; d >= 0; --d) {
	        if (++α[d] < O) {
	            return;                 // normal step
	        }
	        α[d] = 0;                   // carry
	    }
	    α = end();                      // past-the-end sentinel
	}
	static constexpr size_t flatten(index_type const& α) {
	    size_t a = 0;
	    for (int d = 0; d < D; d++) {   //
	        a = O * a + α[d];
	    }
	    return a;
	}
	template<class Function>
	friend auto compose(Function const &f, FwdAutoDiffDouble const &g) {
	    FwdAutoDiffDouble h;
	    double const g0 = g[0];
	    constexpr int Nmax = D * (O - 1);
	    std::array<double, Nmax + 1> fn{};
	    for (int n = 0; n <= Nmax; ++n) {
	        fn[n] = f(g0, n) / factorial(n);
	    }
	    FwdAutoDiffDouble δg = g - FwdAutoDiffDouble(g0);
	    h = FwdAutoDiffDouble(fn[0]);
	    FwdAutoDiffDouble term(1.0);
	    for (int n = 1; n <= Nmax; ++n) {
	        term = term * δg;
	        h += fn[n] * term;
	    }
	    return h;
	}
public:
	friend std::ostream& operator<<(std::ostream& os, FwdAutoDiffDouble const& a) {
	    os << "FwdAutoDiffDouble<O=" << O << ", D=" << D << ">{\n";
	    index_type α = zero();
	    while (α != end()) {
	        os << "  α=(";
	        for (int d = 0; d < D; ++d) {
	            os << α[d];
	            if (d + 1 < D) os << ",";
	        }
	        os << ") -> " << a.C_[flatten(α)] << "\n";
	        increment(α);
	    }
	    os << "}";
	    return os;
	}
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
	constexpr double operator()(std::array<double, D> const &x) const {
		double y = 0.0;
		for (auto α = zero(); α != end(); increment(α)) {
			y += C_[flatten(α)] * pow(x, α);
		}
		return y;
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
	auto operator[](int a) const {
		return C_[a];
	}
	auto operator[](index_type α) const {
		return operator[](flatten(α));
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
	auto operator*(FwdAutoDiffDouble const& B) const {
	    auto const& A = *this;
	    FwdAutoDiffDouble C;
	    C.C_.fill(0.0);
	    for (auto α = zero(); α != end(); increment(α)) {
	        for (auto β = zero(); β != end(); increment(β)) {
	            auto γ = add(α, β);
	            if (max(γ) >= O) {
	                continue;
	            }
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
	friend auto sqrt(FwdAutoDiffDouble const &x) {
		constexpr Rational half(1, 2);
		return compose([](double x0, size_t n) {
			using std::pow;
			using std::sqrt;
			if (n == 0) {
				return sqrt(x0);
			}
			return sqrt(x0) * factorialPower(half, n) / pow(x0, n);
		}, x);
	}

	friend auto cbrt(FwdAutoDiffDouble const &x) {
		constexpr Rational third(1, 3);
		return compose([](double x0, size_t n) {
			using std::pow;
			using std::cbrt;
			if (n == 0) {
				return cbrt(x0);
			}
			return cbrt(x0) * factorialPower(third, n) / pow(x0, n);
		}, x);
	}
	template<Rational power>
	friend auto pow(FwdAutoDiffDouble const &x) {
		if constexpr (power.denominator() == 1) {
		    if constexpr (power.numerator() < 0) {
		        constexpr int m = -power.numerator();
		        return inverse(pow<Rational(m)>(x));
		    } else {
				return compose([](double x, size_t n) {
					return pow<power.numerator()>(x) * factorialPower(power, n) / pow(x, n);
				}, x);
		    }
		} else if constexpr (power.denominator() == 2) {
			return pow<power.numerator()>(x) * sqrt(x);
		} else if constexpr (power.denominator() == 3) {
			return pow<power.numerator()>(x) * cbrt(x);
		} else {
			return compose([](double x, size_t n) {
				using std::pow;
				return pow(x, power - Rational(n)) * factorialPower(power, n);
			}, x);
		}
	}
	friend auto exp(FwdAutoDiffDouble const &x) {
		using std::exp;
		auto const f = compose([](double x, size_t) {
			return exp(x);
		}, x);
		return f;
	}
	friend auto log(FwdAutoDiffDouble const &x) {
		using std::log;
		return compose([](double x, size_t n) {
			if (n == 0) {
				return log(x);
			} else {
				return nonepow(n - 1) * factorial(n - 1) / std::pow(x, n);
			}
		}, x);
	}
	constexpr operator double() const {
		return C_[0];
	}
	template<typename, int, typename >
	friend struct FwdAutoDiff;
private:
	std::array<double, size()> C_;
};
