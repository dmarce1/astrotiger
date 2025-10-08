#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <utility>
#include <vector>

template<int Dim, int Order>
constexpr auto makeMultiIndexTable() {
	// Enumerate non-decreasing index tuples of length k=0..Order.
	std::vector<std::array<int, Order>> indices;
	std::array<int, Order> cur { };

	auto recurse = [&](auto &&self, int depth, int start) -> void {
		if (depth == Order) {
			indices.push_back(cur);
			return;
		}
		for (int i = start; i < Dim; ++i) {
			cur[depth] = i;
			self(self, depth + 1, i); // allow repetition (mixed partials symmetrical)
		}
	};
	recurse(recurse, 0, 0);
	return indices;
}

// -----------------------------------------------------------------------------
// Compile-time factorial and combinatorics
// -----------------------------------------------------------------------------
constexpr int factorial(int n) {
	return (n <= 1) ? 1 : n * factorial(n - 1);
}

// number of combinations with repetition: C(Dim+Order,Order)
constexpr int numDerivatives(int Dim, int Order) {
	int num = 1;
	for (int i = 1; i <= Order; ++i)
		num = num * (Dim + i - 1) / i;
	return num;
}

// -----------------------------------------------------------------------------
// Multi-index enumeration (constexpr version)
// -----------------------------------------------------------------------------
template<int Dim, int Order>
struct MultiIndexTable {
	static constexpr int N = numDerivatives(Dim, Order);
	std::array<std::array<int, Order>, N> data { };

	constexpr MultiIndexTable() {
		int k = 0;
		std::array<int, Order> cur { };
		auto recurse = [&](auto &&self, int depth, int start) -> void {
			if (depth == Order) {
				data[k++] = cur;
				return;
			}
			for (int i = start; i < Dim; ++i) {
				cur[depth] = i;
				self(self, depth + 1, i);
			}
		};
		recurse(recurse, 0, 0);
	}
};

// -----------------------------------------------------------------------------
// FwdAutoDiffDouble
// -----------------------------------------------------------------------------
template<int Dim, int Order>
struct FwdAutoDiffDouble {
	static constexpr int dim() {
		return Dim;
	}
	static constexpr int order() {
		return Order;
	}
	static constexpr int count() {
		return numDerivatives(Dim, Order);
	}

	std::array<double, count()> coeffs { }; // value + derivatives

	// Constructors
	constexpr FwdAutoDiffDouble() :
			coeffs { } {
	}
	constexpr explicit FwdAutoDiffDouble(double v) {
		coeffs.fill(0.0);
		coeffs[0] = v;
	}

	static constexpr FwdAutoDiffDouble variable(double v, int i) {
		FwdAutoDiffDouble x;
		x.coeffs.fill(0.0);
		x.coeffs[0] = v;
		if constexpr (Order >= 1) {
			// index 1..Dim correspond to first derivatives
			x.coeffs[1 + i] = 1.0;
		}
		return x;
	}

	// Access
	constexpr double value() const {
		return coeffs[0];
	}
	constexpr double& value() {
		return coeffs[0];
	}

	// simple arithmetic (first-order correct; extendable)
	constexpr FwdAutoDiffDouble operator+(FwdAutoDiffDouble const &b) const {
		FwdAutoDiffDouble r;
		for (std::size_t n = 0; n < coeffs.size(); ++n) {
			r.coeffs[n] = coeffs[n] + b.coeffs[n];
		}
		return r;
	}

	constexpr FwdAutoDiffDouble operator-(FwdAutoDiffDouble const &b) const {
		FwdAutoDiffDouble r;
		for (std::size_t n = 0; n < coeffs.size(); ++n) {
			r.coeffs[n] = coeffs[n] - b.coeffs[n];
		}
		return r;
	}

	constexpr FwdAutoDiffDouble operator*(double s) const {
		FwdAutoDiffDouble r(*this);
		for (auto &c : r.coeffs) {
			c *= s;
		}
		return r;
	}

	friend constexpr FwdAutoDiffDouble operator*(double s, FwdAutoDiffDouble const &x) {
		return x * s;
	}

	// Product rule for first-order derivatives
	constexpr FwdAutoDiffDouble operator*(FwdAutoDiffDouble const &b) const {
		FwdAutoDiffDouble r;
		r.coeffs[0] = coeffs[0] * b.coeffs[0];
		if constexpr (Order >= 1) {
			for (int i = 0; i < Dim; ++i) {
				r.coeffs[1 + i] = coeffs[1 + i] * b.coeffs[0] + coeffs[0] * b.coeffs[1 + i];
			}
		}
		return r;
	}

	constexpr FwdAutoDiffDouble operator/(FwdAutoDiffDouble const &b) const {
		FwdAutoDiffDouble r;
		double inv = 1.0 / b.coeffs[0];
		r.coeffs[0] = coeffs[0] * inv;
		if constexpr (Order >= 1) {
			for (int i = 0; i < Dim; ++i) {
				r.coeffs[1 + i] = (coeffs[1 + i] * b.coeffs[0] - coeffs[0] * b.coeffs[1 + i]) * inv * inv;
			}
		}
		return r;
	}

	// Common unary functions (first-order)
	friend constexpr FwdAutoDiffDouble sin(FwdAutoDiffDouble const &x) {
		FwdAutoDiffDouble r;
		double s = std::sin(x.coeffs[0]);
		double c = std::cos(x.coeffs[0]);
		r.coeffs[0] = s;
		if constexpr (Order >= 1) {
			for (int i = 0; i < Dim; ++i)
				r.coeffs[1 + i] = x.coeffs[1 + i] * c;
		}
		return r;
	}

	friend constexpr FwdAutoDiffDouble cos(FwdAutoDiffDouble const &x) {
		FwdAutoDiffDouble r;
		double s = std::sin(x.coeffs[0]);
		double c = std::cos(x.coeffs[0]);
		r.coeffs[0] = c;
		if constexpr (Order >= 1) {
			for (int i = 0; i < Dim; ++i)
				r.coeffs[1 + i] = -x.coeffs[1 + i] * s;
		}
		return r;
	}

	friend constexpr FwdAutoDiffDouble exp(FwdAutoDiffDouble const &x) {
		FwdAutoDiffDouble r;
		double e = std::exp(x.coeffs[0]);
		r.coeffs[0] = e;
		if constexpr (Order >= 1) {
			for (int i = 0; i < Dim; ++i)
				r.coeffs[1 + i] = e * x.coeffs[1 + i];
		}
		return r;
	}
};

#include "units.hpp"

template<typename Unit, int Dim, int Order>
struct AutoDiffQuantity {
	using ValueType = FwdAutoDiffDouble<Dim, Order>;
	using UnitType = Unit;

	ValueType data { };

	constexpr AutoDiffQuantity() = default;
	constexpr explicit AutoDiffQuantity(Quantity<Unit> const &q) :
			data(q.value()) {
	}
	constexpr explicit AutoDiffQuantity(double v) :
			data(v) {
	}

	static constexpr AutoDiffQuantity variable(Quantity<Unit> const &q, int i) {
		AutoDiffQuantity r;
		r.data = ValueType::variable(q.value(), i);
		return r;
	}

	constexpr Quantity<Unit> value() const {
		return Quantity<Unit>(data.value());
	}

	// add/sub same unit
	constexpr AutoDiffQuantity operator+(AutoDiffQuantity const &b) const {
		AutoDiffQuantity r;
		r.data = data + b.data;
		return r;
	}

	constexpr AutoDiffQuantity operator-(AutoDiffQuantity const &b) const {
		AutoDiffQuantity r;
		r.data = data - b.data;
		return r;
	}

	// multiply/divide by scalar
	constexpr AutoDiffQuantity operator*(double s) const {
		AutoDiffQuantity r;
		r.data = data * s;
		return r;
	}
	friend constexpr AutoDiffQuantity operator*(double s, AutoDiffQuantity const &x) {
		return x * s;
	}

	// multiply/divide with other quantities (unit algebra)
	template<typename U2>
	constexpr auto operator*(AutoDiffQuantity<U2, Dim, Order> const &b) const {
		using OutUnit = typename UnitMultiply<Unit, U2>::type;
		AutoDiffQuantity<OutUnit, Dim, Order> r;
		r.data = data * b.data;
		return r;
	}

	template<typename U2>
	constexpr auto operator/(AutoDiffQuantity<U2, Dim, Order> const &b) const {
		using OutUnit = typename UnitDivide<Unit, U2>::type;
		AutoDiffQuantity<OutUnit, Dim, Order> r;
		r.data = data / b.data;
		return r;
	}

	// common unary
	friend constexpr AutoDiffQuantity sin(AutoDiffQuantity const &x) {
		static_assert(Unit::isDimensionless, "sin() requires dimensionless argument");
		AutoDiffQuantity r;
		r.data = sin(x.data);
		return r;
	}

	friend constexpr AutoDiffQuantity exp(AutoDiffQuantity const &x) {
		static_assert(Unit::isDimensionless, "exp() requires dimensionless argument");
		AutoDiffQuantity r;
		r.data = exp(x.data);
		return r;
	}
};
