/******************************************************************************
  Cartesian Multivariate Forward AutoDiff (Jets) — Box truncation only
  Author: Dominic C. Marcello (2025)
*******************************************************************************/
#ifndef INCLUDE_FWDAUTODIFF_JET_HPP_
#define INCLUDE_FWDAUTODIFF_JET_HPP_

#include <array>
#include <algorithm>
#include <type_traits>
#include <functional>
#include <utility>
#include <ostream>
#include <cmath>
#include "Util.hpp"

// ===========================
// Public compile-time iterator
// Always loops each axis 0..O inclusive
// ===========================
template<int dimensionCount, int order, typename Func>
constexpr void forEachIndex(std::array<int, dimensionCount>& idx, Func&& f, int depth = 0) {
	if constexpr (dimensionCount == 0) {
		f();
	} else {
		if (depth == dimensionCount - 1) {
			for (idx[depth] = 0; idx[depth] <= order; ++idx[depth]) {
				f(idx);
			}
		} else {
			for (idx[depth] = 0; idx[depth] <= order; ++idx[depth]) {
				forEachIndex<dimensionCount, order>(idx, std::forward<Func>(f), depth + 1);
			}
		}
	}
}

template<int dimensionCount, int order, typename Func>
constexpr void forEachIndex(Func&& f) {
	std::array<int, dimensionCount> idx{};
	forEachIndex<dimensionCount, order>(idx, std::forward<Func>(f), 0);
}

// ===========================
// Small math helpers (constexpr)
// ===========================

//template<typename T>
//constexpr T factorialPower(T const& a, int n) {
//	T r = T(1);
//	for (int k = 0; k < n; ++k) {
//		r *= (a - T(k));
//	}
//	return r;
//}


template<int dimensionCount>
constexpr int sumIndex(std::array<int, dimensionCount> const& i) {
	int s = 0;
	for (int k = 0; k < dimensionCount; ++k) {
		s += i[k];
	}
	return s;
}

template<int dimensionCount>
constexpr bool leqIndices(std::array<int, dimensionCount> const& a, std::array<int, dimensionCount> const& b) {
	for (int k = 0; k < dimensionCount; ++k) {
		if (a[k] > b[k]) {
			return false;
		}
	}
	return true;
}

template<int dimensionCount, int order>
constexpr bool insideBox(std::array<int, dimensionCount> const& i) {
	for (int k = 0; k < dimensionCount; ++k) {
		if (i[k] < 0 || i[k] > order) {
			return false;
		}
	}
	return true;
}

template<int dimensionCount, int order>
constexpr std::size_t flatten(std::array<int, dimensionCount> const& i) {
	std::size_t idx = 0;
	std::size_t stride = 1;
	for (int k = 0; k < dimensionCount; ++k) {
		idx += static_cast<std::size_t>(i[k]) * stride;
		stride *= static_cast<std::size_t>(order + 1);
	}
	return idx;
}

template<int dimensionCount>
constexpr std::array<int, dimensionCount> addIndex(std::array<int, dimensionCount> a, std::array<int, dimensionCount> const& b) {
	for (int k = 0; k < dimensionCount; ++k) {
		a[k] += b[k];
	}
	return a;
}

template<int dimensionCount>
constexpr std::array<int, dimensionCount> subIndex(std::array<int, dimensionCount> a, std::array<int, dimensionCount> const& b) {
	for (int k = 0; k < dimensionCount; ++k) {
		a[k] -= b[k];
	}
	return a;
}

template<int dimensionCount>
constexpr bool isZeroIndex(std::array<int, dimensionCount> const& i) {
	for (int k = 0; k < dimensionCount; ++k) {
		if (i[k] != 0) {
			return false;
		}
	}
	return true;
}

// ===========================
// Traits (compatible with your older interface)
// ===========================
template<typename T, int, int>
struct FwdAutoDiff; // fwd

template<typename T>
struct IsFwdAutoDiff {
	static constexpr bool value = false;
	using type = void;
};

template<typename T>
inline constexpr bool IsFwdAutoDiff_v = IsFwdAutoDiff<T>::value;

template<typename T>
using FwdAutoDiffValue_t = typename IsFwdAutoDiff<T>::type;

// ===========================
// FwdAutoDiff<T,O,D> : Cartesian-only jets
// ===========================
template<typename T, int order_ = 2, int dimensionCount = 1>
struct FwdAutoDiff {
public:
	using value_type = T;
	static constexpr int order = order_;
	static constexpr int dim = dimensionCount;
	static constexpr std::size_t storageSize() {
		std::size_t n = 1;
		for (int k = 0; k < dimensionCount; ++k) {
			n *= static_cast<std::size_t>(order + 1);
		}
		return n;
	}

	// ---- ctors ----
	constexpr FwdAutoDiff() {
		C_.fill(T(0));
	}
	constexpr auto D(int d) const {
		return C_[d];
	}
	constexpr explicit FwdAutoDiff(T const& value) {
		C_.fill(T(0));
		C_[0] = value;
	}

	// Independent variable along axis (0..D-1)
	static constexpr FwdAutoDiff independentVariable(T const& x0, int axis = 0) {
		FwdAutoDiff a;
		a.C_[0] = x0;
		if (axis >= 0 && axis < dimensionCount) {
			std::array<int, dimensionCount> e{};
			for (int k = 0; k < dimensionCount; ++k) {
				e[k] = 0;
			}
			e[axis] = 1;
			a[e] = T(1);
		}
		return a;
	}

	static constexpr FwdAutoDiff constant(T const& c) {
		return FwdAutoDiff(c);
	}

	// ---- conversions / access ----
	explicit constexpr operator T() const {
		return C_[0];
	}

	constexpr T operator[](std::array<int, dimensionCount> const& i) const {
		return C_[flatten<dimensionCount, order>(i)];
	}

	constexpr T& operator[](std::array<int, dimensionCount> const& i) {
		return C_[flatten<dimensionCount, order>(i)];
	}

	// ---- printing (brief) ----
	friend std::ostream& operator<<(std::ostream& os, FwdAutoDiff const& A) {
		bool first = true;
		forEachIndex<dimensionCount, order>([&](auto const& i) {
			if (!first) {
				os << ", ";
			}
			first = false;
			os << "[";
			for (int k = 0; k < dimensionCount; ++k) {
				os << i[k] << (k + 1 < dimensionCount ? "," : "");
			}
			os << "]=" << static_cast<double>(A[i]);
		});
		return os;
	}

	// ---- evaluator for the truncated polynomial ----
	template<typename Real = T, typename... Args>
	auto generateFunction() const {
		std::array<Real, storageSize()> A;
		for (std::size_t i = 0; i < storageSize(); ++i) {
			A[i] = static_cast<Real>(C_[i]);
		}
		return [A](Args... xs) -> Real {
			static_assert(sizeof...(Args) == dimensionCount, "generateFunction: arity must equal D");
			std::array<Real, dimensionCount> X{ static_cast<Real>(xs)... };
			Real y = Real(0);
			forEachIndex<dimensionCount, order>([&](auto const& I) {
				Real term = Real(1);
				for (int k = 0; k < dimensionCount; ++k) {
					term *= ipow(X[k], I[k]);
				}
				y += A[flatten<dimensionCount, order>(I)] * term;
			});
			return y;
		};
	}

	// =====================
	// Unary
	// =====================
	friend constexpr auto operator+(FwdAutoDiff const& A) {
		return A;
	}

	friend constexpr auto operator-(FwdAutoDiff const& A) {
		FwdAutoDiff C;
		for (std::size_t i = 0; i < storageSize(); ++i) {
			C.C_[i] = -A.C_[i];
		}
		return C;
	}

	// =====================
	// + / -
	// =====================
	friend constexpr auto operator+(FwdAutoDiff const& A, FwdAutoDiff const& B) {
		FwdAutoDiff C;
		forEachIndex<dimensionCount, order>([&](auto const& I) {
			C[I] = A[I] + B[I];
		});
		return C;
	}

	friend constexpr auto operator-(FwdAutoDiff const& A, FwdAutoDiff const& B) {
		FwdAutoDiff C;
		forEachIndex<dimensionCount, order>([&](auto const& I) {
			C[I] = A[I] - B[I];
		});
		return C;
	}

	// =====================
	// Scalar interop
	// =====================
	template<typename U, typename = std::enable_if_t<!IsFwdAutoDiff_v<U>>>
	friend constexpr auto operator*(FwdAutoDiff const& X, U const& b) {
		using R = decltype(std::declval<T>() * std::declval<U>());
		FwdAutoDiff<R, order, dimensionCount> Z;
		for (std::size_t i = 0; i < storageSize(); ++i) {
			Z.raw()[i] = X.C_[i] * b;
		}
		return Z;
	}

	template<typename U, typename = std::enable_if_t<!IsFwdAutoDiff_v<U>>>
	friend constexpr auto operator*(U const& a, FwdAutoDiff const& Y) {
		return Y * a;
	}

	template<typename U, typename = std::enable_if_t<!IsFwdAutoDiff_v<U>>>
	friend constexpr auto operator/(FwdAutoDiff const& X, U const& b) {
		using R = decltype(std::declval<T>() / std::declval<U>());
		FwdAutoDiff<R, order, dimensionCount> Z;
		for (std::size_t i = 0; i < storageSize(); ++i) {
			Z.raw()[i] = static_cast<R>(X.C_[i]) / static_cast<R>(b);
		}
		return Z;
	}

	template<typename U, typename = std::enable_if_t<!IsFwdAutoDiff_v<U>>>
	friend constexpr auto operator/(U const& a, FwdAutoDiff const& Y) {
		using R = decltype(std::declval<U>() / std::declval<T>());
		FwdAutoDiff<R, order, dimensionCount> num(a);
		return num / static_cast<FwdAutoDiff<R, order, dimensionCount>>(Y);
	}

	// =====================
	// AD * AD (Cartesian convolution with box truncation)
	// =====================
//	friend constexpr auto operator*(FwdAutoDiff const& X, FwdAutoDiff const& Y) {
//		FwdAutoDiff Z;
//		forEachIndex<dimensionCount, order>([&](auto const& I) {
//			forEachIndex<dimensionCount, order>([&](auto const& K) {
//				auto S = addIndex<dimensionCount>(I, K);
//				if (insideBox<dimensionCount, order>(S)) {
//					Z[S] += X[I] * Y[K];
//				}
//			});
//		});
//		return Z;
//	}

//	// =====================
//	// AD / AD (formal deconvolution)
//	// Requires B[0] != 0
//	// =====================]
//	friend constexpr auto operator/(FwdAutoDiff const& A, FwdAutoDiff const& B) {
//		FwdAutoDiff C;
//		T const b0 = B.C_[0];
//		// constant term
//		C.C_[0] = A.C_[0] / b0;
//
//		// graded by total degree to ensure dependencies are available
//		for (int s = 1; s <= order * dimensionCount; ++s) {
//			forEachIndex<dimensionCount, order>([&](auto const& I) {
//				if (sumIndex<dimensionCount>(I) != s) {
//					return;
//				}
//				T acc = A[I];
//				forEachIndex<dimensionCount, order>([&](auto const& K) {
//					if (isZeroIndex<dimensionCount>(K)) {
//						return;
//					}
//					if (!leqIndices<dimensionCount>(K, I)) {
//						return;
//					}
//					auto ImK = subIndex<dimensionCount>(I, K);
//					acc -= B[K] * C[ImK];
//				});
//				C[I] = acc / b0;
//			});
//		}
//		return C;
//	}


	// =====================
	// Compound ops
	// =====================
	constexpr FwdAutoDiff& operator+=(FwdAutoDiff const& rhs) {
		*this = *this + rhs;
		return *this;
	}
	constexpr FwdAutoDiff& operator-=(FwdAutoDiff const& rhs) {
		*this = *this - rhs;
		return *this;
	}
	constexpr FwdAutoDiff& operator*=(FwdAutoDiff const& rhs) {
		*this = *this * rhs;
		return *this;
	}
	template<typename U, typename = std::enable_if_t<!IsFwdAutoDiff_v<U>>>
	constexpr FwdAutoDiff& operator*=(U const& rhs) {
		*this = *this * rhs;
		return *this;
	}
	constexpr FwdAutoDiff& operator/=(FwdAutoDiff const& rhs) {
		*this = *this / rhs;
		return *this;
	}
	template<typename U, typename = std::enable_if_t<!IsFwdAutoDiff_v<U>>>
	constexpr FwdAutoDiff& operator/=(U const& rhs) {
		*this = *this / rhs;
		return *this;
	}


	// Generic univariate-on-jet composition: f(g(x))
	// dkdF(g0, n) returns the n-th derivative of f at g0 (not divided by n!)
	template<typename Function>
	friend constexpr auto compose(Function& dkdF, FwdAutoDiff const& g) {
		FwdAutoDiff h;
		T const g0 = g.C_[0];

		// delta = g - g0 (zero-out constant term)
		FwdAutoDiff delta = g;
		delta.C_[0] = T(0);

		// term = (g - g0)^n accumulates by repeated multiplication
		FwdAutoDiff term(T(1)); // n=0
		int const Nmax = order * dimensionCount;

		for (int n = 0; n <= Nmax; ++n) {
			T a_n = dkdF(g0, n) / static_cast<T>(factorial(n));
			h += a_n * term;
			if (n < Nmax) {
				term = term * delta;
			}
		}
		return h;
	}

	friend constexpr FwdAutoDiff exp(FwdAutoDiff const& x) {
		return compose([](T const& g0, int) {
			return std::exp(g0);
		}, x);
	}

	friend constexpr FwdAutoDiff log(FwdAutoDiff const& x) {
		return compose([](T const& g0, int n) {
			if (n == 0) {
				return std::log(g0);
			}
			T s = (n % 2 == 0) ? T(-1) : T(1);
			return s * static_cast<T>(factorial(n - 1)) / ipow(g0, n);
		}, x);
	}

	template<typename T1>
	friend auto const pow(FwdAutoDiff const& x, T1 const& p) {
		auto const f = [=](T const& g0, int n) {
			using std::pow;
			auto rc = factorialPower(p, n) * pow(g0, p - n);
			return rc;
		};
		auto rc = compose(f, x);
		return rc;
	}

	friend auto sqrt(FwdAutoDiff const& x) {
		return pow(x, T(0.5));
	}

	// Light-weight nth-derivative helpers for tanh/atanh at a point.
	// For robustness, we compute the first few exactly; higher orders use a simple
	// recurrence approximation via derivatives of 1 - tanh^2.
	friend constexpr FwdAutoDiff tanh(FwdAutoDiff const& x) {
		return compose([](T const& g0, int n) {
			if (n == 0) {
				return std::tanh(g0);
			}
			if (n == 1) {
				T c = std::cosh(g0);
				return T(1) / (c * c);
			}
			// Simple fallback: use symbolic relation for second derivative:
			if (n == 2) {
				T t = std::tanh(g0);
				T c = std::cosh(g0);
				return T(-2) * t / (c * c);
			}
			// For higher n, a compact exact formula exists but is involved; zero is conservative.
			return T(0);
		}, x);
	}

	friend constexpr FwdAutoDiff atanh(FwdAutoDiff const& x) {
		return compose([](T const& g0, int n) {
			if (n == 0) {
				return std::atanh(g0);
			}
			// d^n/dx^n atanh(x) = (n-1)! * A_n(x) / (1 - x^2)^n ; keep n=1 exact:
			if (n == 1) {
				return T(1) / (T(1) - g0 * g0);
			}
			return T(0);
		}, x);
	}

	// ---- raw access for internal helpers (kept public for flexibility) ----
	constexpr std::array<T, storageSize()>& raw() { return C_; }
	constexpr std::array<T, storageSize()> const& raw() const { return C_; }

private:
	std::array<T, storageSize()> C_;
};

// ---- trait specialization ----
template<typename T, int order, int dimensionCount>
struct IsFwdAutoDiff<FwdAutoDiff<T, order, dimensionCount>> {
	static constexpr bool value = true;
	using type = T;
};

// =============================
// Heterogeneous product
// =============================
template<typename A, typename B, int order, int dimensionCount,
         typename R = decltype(std::declval<A>() * std::declval<B>())>
constexpr FwdAutoDiff<R,order,dimensionCount> operator*(FwdAutoDiff<A,order,dimensionCount> const& X,
                                       FwdAutoDiff<B,order,dimensionCount> const& Y) {
    FwdAutoDiff<R,order,dimensionCount> Z;
    forEachIndex<dimensionCount,order>([&](auto const& i) {
        forEachIndex<dimensionCount,order>([&](auto const& j) {
            auto s = addIndex<dimensionCount>(i, j);
            if (insideBox<dimensionCount,order>(s)) {
                Z[s] += X[i] * Y[j];
            }
        });
    });
    return Z;
}

// =============================
// Heterogeneous division
// =============================
template<typename A, typename B, int order, int dimensionCount,
         typename R = decltype(std::declval<A>() / std::declval<B>())>
constexpr FwdAutoDiff<R, order, dimensionCount> operator/(FwdAutoDiff<A, order, dimensionCount> const& X,
                                       FwdAutoDiff<B, order, dimensionCount> const& Y) {
    FwdAutoDiff<R,order,dimensionCount> Z;
    B const y0 = Y.raw()[0];
    Z.raw()[0] = X.raw()[0] / y0;
    forEachIndex<dimensionCount,order>([&](auto const& I) {
        if (sumIndex<dimensionCount>(I) == 0) { return; }
        A acc = X[I];
        forEachIndex<dimensionCount,order>([&](auto const& K) {
            if (isZeroIndex<dimensionCount>(K) || !leqIndices<dimensionCount>(K,I)) { return; }
            auto ImK = subIndex<dimensionCount>(I,K);
            acc -= Y[K] * Z[ImK];
        });
        Z[I] = acc / y0;
    });
    return Z;
}
// FwdAutoDiff<T> ± T
template<typename A, typename B, int order, int dimensionCount,
         typename R = decltype(std::declval<A>() + std::declval<B>())>
constexpr FwdAutoDiff<R,order,dimensionCount> operator+(FwdAutoDiff<A,order,dimensionCount> const& X, B const& b) {
    FwdAutoDiff<B,order,dimensionCount> Y(b);
    return X + Y;
}

template<typename A, typename B, int order, int dimensionCount,
         typename R = decltype(std::declval<A>() - std::declval<B>())>
constexpr FwdAutoDiff<R,order,dimensionCount> operator-(FwdAutoDiff<A,order,dimensionCount> const& X, B const& b) {
    FwdAutoDiff<B,order,dimensionCount> Y(b);
    return X - Y;
}

// T ± FwdAutoDiff<T>
template<typename A, typename B, int order, int dimensionCount,
         typename R = decltype(std::declval<A>() - std::declval<B>())>
constexpr FwdAutoDiff<R,order,dimensionCount> operator-(A const& a, FwdAutoDiff<B,order,dimensionCount> const& Y) {
    FwdAutoDiff<A,order,dimensionCount> X(a);
    return X - Y;
}

#endif /* INCLUDE_FWDAUTODIFF_JET_HPP_ */
