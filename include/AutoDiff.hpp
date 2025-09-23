#pragma once

#include <algorithm>
#include "Multidices.hpp"
#include "BellPolynomial.hpp"

// Α,α,Β,β,Γ,γ,Δ,δ,Ε,ε,Ζ,ζ,Η,η,Θ,θ,ϑ,Ι,ι,Κ,κ,ϰ,Λ,λ,Μ,μ,Ν,μ,Ξ,ξ,Ο,ο,Π,π,ϖ,Ρ,ρ,ϱ,Σ,σ,ς,Τ,τ,Υ,υ,Φ,φ,ϕ,Χ,χ,Ψ,ψ,Ω,ω

//template<typename T, int N, int ...Ns>
//struct PolyN {
//	using value_type = typename std::conditional<(sizeof...(Ns) > 0), PolyN<T, Ns>, T>::type;
//	using reference = typename std::conditional<(sizeof...(Ns) > 0), PolyN<T, Ns>&, T&>::type;
//	using const_reference = typename std::conditional<(sizeof...(Ns) > 0), PolyN<T, Ns> const&, T const&>::type;
//
//	const_reference operator[](int i) const {
//		return coefficients_[i];
//	}
//	reference operator[](int i) {
//		return coefficients_[i];
//	}
//	PolyN operator+() const {
//		return *this;
//	}
//	PolyN operator-() const {
//		PolyN const &A = *this;
//		PolyN B;
//		for (int i = 0; i <= N; i++) {
//			B[i] = -A[i];
//		}
//		return B;
//	}
//	PolyN operator+(PolyN const &A) const {
//		PolyN C;
//		PolyN const &B = *this;
//		for (int i = 0; i <= N; i++) {
//			C[i] = A[i] + B[i];
//		}
//		return C;
//	}
//	PolyN operator-(PolyN const &A) const {
//		PolyN C;
//		PolyN const &B = *this;
//		for (int i = 0; i <= N; i++) {
//			C[i] = A[i] - B[i];
//		}
//		return C;
//	}
//	template<typename S, int M, int ...Ms>
//	auto operator*(PolyN<S, M, Ms...> const &A) const {
//		using ReturnT1 = decltype(T{} * S{});
//		PolyN<ReturnT1, N + M, (Ns + Ms)...> C;
//		PolyN const &B = *this;
//		for (int n = 0; n <= N + M; n++) {
//			C[n] = T(0);
//			for (int k = 0; k <= M; k++) {
//				C[n] += A[n - k] * B[k];
//			}
//		}
//		return C;
//	}
//	template<typename S>
//	auto operator*(const S &a) const {
//		using ReturnT1 = decltype(T{} * S{});
//		PolyN<ReturnT1, N, Ns...> C = *this;
//			return C;
//	}
//private:
//	std::array<value_T, N + 1> coefficients_;
//
//};
//
inline int stirling2(int n, int k) {
	static thread_local std::vector<std::vector<int>> S2;
	while (int(S2.size()) > n) {
		std::vector<int> s2(S2.size() + 2);
		int const n1 = S2.size();
		if (n1 == 0) {
			s2[0] = 1;
			s2[1] = 0;
		} else {
			for (int k1 = 0; k1 <= n1; k1++) {
				s2[k1] = k1 * S2[n1 - 1][k1] + S2[n1 - 1][k1 - 1];
			}
			s2.back() = 0;
		}
		S2.push_back(std::move(s2));
	}
	return S2[n][k];
}

template<typename T1, int O, int D>
struct AutoDiff {
	static constexpr size_t size() {
		return binco(O + D - 1, D);
	}
private:
	std::array<T1, size()> C_;
public:
	explicit operator T1() const {
		return C_[0];
	}
	T1 operator[](Multidices<D> const &mI) const {
		return C_[int(mI)];
	}
	T1& operator[](Multidices<D> const &mI) {
		return C_[int(mI)];
	}
	template<typename Real = T1>
	auto generateFunction() {
		std::array<Real, size()> A;
		std::copy(C_.begin(), C_.end(), A.begin());
		return [A](auto ...xs) {
			std::array<Real, D> X = { xs... };
			std::reverse(X.begin(), X.end());
			Real y = Real(0);
			for (Multidices<D> N = 0; N < abs(O); N++) {
				y += A[N] * pow(X, N);
			}
			return y;
		};
	}
	friend AutoDiff operator+(AutoDiff const &A) {
		return A;
	}
	friend AutoDiff operator-(AutoDiff const &A) {
		AutoDiff C;
		std::transform(A.C_.begin(), A.C_.end(), C.C_.begin(), std::negate<T1>());
		return C;
	}
	friend AutoDiff operator+(AutoDiff const &A, AutoDiff const &B) {
		AutoDiff C;
		std::transform(A.C_.begin(), A.C_.end(), B.C_.begin(), C.C_.begin(), std::plus<T1>());
		return C;
	}
	friend AutoDiff operator-(AutoDiff const &A, AutoDiff const &B) {
		AutoDiff C;
		std::transform(A.C_.begin(), A.C_.end(), B.C_.begin(), C.C_.begin(), std::minus<T1>());
		return C;
	}
	friend AutoDiff operator*(AutoDiff A, T1 const &b) {
		for (auto &a : A.C_) {
			a *= b;
		}
		return A;
	}
	friend AutoDiff operator*(T1 const &a, AutoDiff const &B) {
		return B * a;
	}
	AutoDiff& operator*=(AutoDiff const &A) {
		*this = *this * A;
		return *this;
	}
	AutoDiff& operator*=(T1 const &a) {
		*this = *this * a;
		return *this;
	}
	AutoDiff& operator/=(AutoDiff const &A) {
		*this = *this / A;
		return *this;
	}
	AutoDiff& operator/=(T1 const &a) {
		*this = *this / a;
		return *this;
	}
	AutoDiff& operator+=(AutoDiff const &A) {
		*this = *this + A;
		return *this;
	}
	AutoDiff& operator-=(AutoDiff const &A) {
		*this = *this - A;
		return *this;
	}
	friend std::ostream& operator<<(std::ostream &os, AutoDiff const &jet) {
		std::function<std::ostream& (std::ostream&, T1 const*, int, int)> print;
		print = [&print](std::ostream &os, T1 const *ptr, int n, int dim) -> std::ostream& {
			if (dim == 0) {
				os << print2string("{");
				for (int i = 0; i < n; i++) {
					os << print2string("%11.3e ", double(ptr[i]));
					if (i + 1 < n) {
						os << print2string(", ");
					}
				}
				os << print2string("}");
			} else {
				for (int i = 0; i < n; i++) {
					os << print2string("%c = %2i : ", 'j' + dim, i);
					print(os, ptr + binco(i + dim, 1 + dim), i + 1, dim - 1);
				}
			}
			os << print2string("\n");
			return os;
		};
		print(os, jet.C_.data(), O, D - 1);
		return os;
	}
	friend AutoDiff compose(std::function<T1(T1, int)> const &f, AutoDiff const &gj) {
		AutoDiff H;
		std::array<T1, O> dfdg;
		T1 const g0 = gj[0];
		for (int k = 0; k < O; k++) {
			dfdg[k] = f(g0, k);
		}
		H[0] = dfdg[0];
		for (Multidices<D> N = 1; abs(N) < O; N++) {
			T1 h = T1(0);
			for (Multidices<1> k = 1; k <= abs(N); k++) {
				auto const Bnk = multivariateBellPolynomial<D, 1>(N, k);
				T1 Bg = T1(0);
				for (auto const &term : Bnk) {
					T1 product = T1(1);
					for (auto const &factor : term) {
						auto const J = factor.first;
						auto const Kj = factor.second;
						product *= ipow(gj[J], Kj[0]) / T1(factorial(Kj[0]));
					}
					Bg += product;
				}
				h += dfdg[k] * Bg;
			}
			H[N] = h;
		}
		return H;
	}
	friend AutoDiff operator/(AutoDiff const &R, T1 const &A) {
		return R * (T1(1) / A);
	}
//	friend AutoDiff operator/(AutoDiff const &y, AutoDiff const &x) {
//		return y * (T1(1) / x);
//	}
	friend AutoDiff pow(AutoDiff const &x, T1 const &k) {
		return compose([k](T1 const &x, int n) {
			return factorialPower(k, n) * pow(x, k - n);
		}, x);
	}
	friend AutoDiff sqr(AutoDiff const &x) {
		return x * x;
	}
	friend AutoDiff sqrt(AutoDiff const &x) {
		return compose([](T1 const &x, int n) {
			return factorialPower(0.5, n) * std::sqrt(x) / ipow(x, n);
		},x);
	}
	friend AutoDiff exp(AutoDiff const &x) {
		return compose([](T1 const &x, int) {
			return std::exp(x);
		}, x);
	}
	friend AutoDiff log(AutoDiff const &x) {
		return compose([](T1 const &x, int) {
			return std::log(x);
		}, x);
	}
	friend AutoDiff tanh(AutoDiff const &x) {
		return compose([](T1 const &x, int) {
			return std::tanh(x);
		}, x);
	}
	friend AutoDiff atanh(AutoDiff const &x) {
		return compose([](T1 const &x, int) {
			return std::atanh(x);
		}, x);
	}
	AutoDiff(T1 value = T1(0)) {
		C_.fill(0);
		C_[0] = value;
	}
	AutoDiff(T1 const &value, int dim) {
		C_.fill(0);
		C_[0] = value;
		int const i = binco(1 + dim, dim);
		C_[i] = T1(1);
	}
};

template<typename T1, typename T2, int O, int D, typename T3 = decltype(T1()*T2())>
auto operator*(AutoDiff<T1, O, D> const &A, AutoDiff<T2, O, D> const &B) {
	AutoDiff<T3, O, D> C;
	for (Multidices<D> n = 0; abs(n) < O; n++) {
		for (Multidices<D> k = 0; abs(n + k) < O; k++) {
			C[n + k] += A[n] * B[k];
		}
	}
	return C;
}


template<typename T1, typename T2, int O, int D, typename T3 = decltype(T1()/T2())>
AutoDiff<T3, O, D> operator/(AutoDiff<T1, O, D> const &y, AutoDiff<T2, O, D> const &x) {
	return y * compose([](T1 const &x, int n) {
		return factorialPower(-1, n) * ipow(x, -(n + 1));
	}, x);
}
