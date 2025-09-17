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
//		using ReturnType = decltype(T{} * S{});
//		PolyN<ReturnType, N + M, (Ns + Ms)...> C;
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
//		using ReturnType = decltype(T{} * S{});
//		PolyN<ReturnType, N, Ns...> C = *this;
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

template<typename Type, int O, int D>
struct AutoDiff {
	static constexpr size_t size() {
		return binco(O + D - 1, D);
	}
private:
	std::array<Type, size()> C_;
public:
	explicit operator Type() const {
		return C_[0];
	}
	Type operator[](Multidices<D> const &mI) const {
		return C_[int(mI)];
	}
	Type& operator[](Multidices<D> const &mI) {
		return C_[int(mI)];
	}
	template<typename Real = Type>
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
		std::transform(A.C_.begin(), A.C_.end(), C.C_.begin(), std::negate<Type>());
		return C;
	}
	friend AutoDiff operator+(AutoDiff const &A, AutoDiff const &B) {
		AutoDiff C;
		std::transform(A.C_.begin(), A.C_.end(), B.C_.begin(), C.C_.begin(), std::plus<Type>());
		return C;
	}
	friend AutoDiff operator-(AutoDiff const &A, AutoDiff const &B) {
		AutoDiff C;
		std::transform(A.C_.begin(), A.C_.end(), B.C_.begin(), C.C_.begin(), std::minus<Type>());
		return C;
	}
	friend AutoDiff operator*(AutoDiff const &A, AutoDiff const &B) {
		AutoDiff C;
		for (Multidices<D> n = 0; abs(n) < O; n++) {
			for (Multidices<D> k = 0; abs(n + k) < O; k++) {
				C[n + k] += A[n] * B[k];
			}
		}
		return C;
	}
	friend AutoDiff operator*(AutoDiff A, Type const &b) {
		for (auto &a : A.C_) {
			a *= b;
		}
		return A;
	}
	friend AutoDiff operator*(Type const &a, AutoDiff const &B) {
		return B * a;
	}
	AutoDiff& operator*=(AutoDiff const &A) {
		*this = *this * A;
		return *this;
	}
	AutoDiff& operator*=(Type const &a) {
		*this = *this * a;
		return *this;
	}
	AutoDiff& operator/=(AutoDiff const &A) {
		*this = *this / A;
		return *this;
	}
	AutoDiff& operator/=(Type const &a) {
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
		std::function<std::ostream& (std::ostream&, Type const*, int, int)> print;
		print = [&print](std::ostream &os, Type const *ptr, int n, int dim) -> std::ostream& {
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
	friend AutoDiff compose(std::function<Type(Type, int)> const &f, AutoDiff const &gj) {
		AutoDiff H;
		std::array<Type, O> dfdg;
		Type const g0 = gj[0];
		for (int k = 0; k < O; k++) {
			dfdg[k] = f(g0, k);
		}
		H[0] = dfdg[0];
		for (Multidices<D> N = 1; abs(N) < O; N++) {
			Type h = Type(0);
			for (Multidices<1> k = 1; k <= abs(N); k++) {
				auto const Bnk = multivariateBellPolynomial<D, 1>(N, k);
				Type Bg = Type(0);
				for (auto const &term : Bnk) {
					Type product = Type(1);
					for (auto const &factor : term) {
						auto const J = factor.first;
						auto const Kj = factor.second;
						product *= ipow(gj[J], Kj[0]) / Type(factorial(Kj[0]));
					}
					Bg += product;
				}
				h += dfdg[k] * Bg;
			}
			H[N] = h;
		}
		return H;
	}
	//y^-1(y(x)) = x
//	AutoDiff invert() {
//		AutoDiff const& F = *this;
//		AutoDiff H;
//		H[0] = Type(0);
//		for (Multidices<D> N = 1; abs(N) <= 1; N++) {
//			H[N] = Type(1) / F[N];
//		}
//		for (Multidices<D> N = 2; abs(N) < O; N++) {
//			Type h = Type(0);
//			for (Multidices<1> k = 1; k <= abs(N); k++) {
//				auto const Bnk = multivariateBellPolynomial<D, 1>(N, k);
//				Type Bg = Type(0);
//				for (auto const &term : Bnk) {
//					Type product = Type(1);
//					for (auto const &factor : term) {
//						auto const J = factor.first;
//						auto const Kj = factor.second;
//						product *= ipow(F[J], Kj[0]) / factorial(Kj[0]);
//					}
//					Bg += product;
//				}
//				h += dfdg[k] * Bg;
//			}
//			H[N] = h;
//		}
//		return H;
//	}
	friend AutoDiff operator/(AutoDiff const &R, Type const &A) {
		return R * (Type(1) / A);
	}
	friend AutoDiff operator/(Type const &y, AutoDiff const &x) {
		return y * compose([](Type const &x, int n) {
			return factorialPower(-1, n) * pow(x, Type(-(n + 1)));
		}, x);
	}
	friend AutoDiff operator/(AutoDiff const &y, AutoDiff const &x) {
		return y * (Type(1) / x);
	}
	friend AutoDiff pow(AutoDiff const &x, Type const &k) {
		return compose([k](Type const &x, int n) {
			return factorialPower(k, n) * pow(x, k - n);
		}, x);
	}
	friend AutoDiff sqr(AutoDiff const &x) {
		return x * x;
	}
	friend AutoDiff sqrt(AutoDiff const &x) {
		return compose([](Type const &x, int n) {
			return factorialPower(0.5, n) * std::sqrt(x) / ipow(x, n);
		},x);
	}
	friend AutoDiff exp(AutoDiff const &x) {
		return compose([](Type const &x, int) {
			return std::exp(x);
		}, x);
	}
	friend AutoDiff log(AutoDiff const &x) {
		return compose([](Type const &x, int) {
			return std::log(x);
		}, x);
	}
	friend AutoDiff tanh(AutoDiff const &x) {
		return compose([](Type const &x, int) {
			return std::tanh(x);
		}, x);
	}
	friend AutoDiff atanh(AutoDiff const &x) {
		return compose([](Type const &x, int) {
			return std::atanh(x);
		}, x);
	}
	AutoDiff(Type value = Type(0)) {
		C_.fill(0);
		C_[0] = value;
	}
	AutoDiff(Type const &value, int dim) {
		C_.fill(0);
		C_[0] = value;
		int const i = binco(1 + dim, dim);
		C_[i] = Type(1);
	}
};

//template<typename T, typename U, int N>
//struct PowerSeries {
//	template<typename Function>
//	PowerSeries(Function const &f) {
//		for (int k = 0; k < N; k++) {
//			Cn[k] = f(0.0, k) / factorial(k);
//		}
//	}
//	PowerSeries(std::array<T, N> const &cc) {
//		for (int k = 0; k < N; k++) {
//			Cn[k] = cc[k] / factorial(k);
//		}
//	}
//	U operator()(U x) const {
//		U sum = U(Cn[0]);
//		U xn = x;
//		for (int k = 1; k < N - 1; k++) {
//			sum += Cn[k] * xn;
//			xn *= x;
//		}
//		sum += Cn.back() * xn;
//		return sum;
//	}
//	template<typename, int, int>
//	friend struct Jet;
//private:
//	std::array<T, N> Cn;
//};
///
//template<int N, int D1, int D2>
//struct MultivariateBellPolynomials {
//	using JIndex = Multidices<D1>;
//	using KIndex = Multidices<D2>;
//	MultivariateBellPolynomials() {
//		for (int n = 1; n < N; n++) {
//			std::vector<std::vector<int>> Kjs;
//			std::vector<int> Kj(n, 0);
//			bool done = false;
//			while (!done) {
//				int sum = 0;
//				for (int j = 0; j < n; j++) {
//					sum += (j + 1) * Kj[j];
//				}
//				if (sum == n) {
//					Kjs.push_back(Kj);
//				}
//				int j = 0;
//				while (Kj[j]++ >= (n / (j + 1))) {
//					Kj[j] = 0;
//					j++;
//					if (j == n) {
//						done = true;
//						break;
//					}
//				}
//			}
//			for (auto const &Kj : Kjs) {
//				for (int j = 0; j <= n; j++) {
//					KIndex K0 = 0;
//					std::vector<KIndex> Ks;
//					for (KIndex K = k; abs(K) == k; K++) {
//						K0 += K;
//						JIndex J0 = 0;
//						for (JIndex J = j; abs(J) == j; J++) {
//							J0 += J;
//							if ((K0 <= n) && (J0 <= n)) {
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//}
//;

//template<typename T, int O, int D>
//struct Jet {
//	using IndexType = Multidices<D>;
//	template<typename U>
//	PowerSeries<T, U, O> generatePowerSeries() {
//		PowerSeries<T, U, O> P(data_.begin(), data_.end());
//	}
//	T operator[](IndexType const &I) const {
//		return data_[int(I)];
//	}
//	T& operator[](int i) {
//		return data_[i];
//	}
//	T operator[](int i) const {
//		return data_[i];
//	}
//	T& operator[](IndexType const &I) {
//		return data_[int(I)];
//	}
//	Jet& operator*=(Jet const &other) {
//		*this = *this * other;
//		return *this;
//	}
//	Jet& operator*=(T const &other) {
//		*this = *this * other;
//		return *this;
//	}
//	Jet& operator/=(T const &other) {
//		*this = *this / other;
//		return *this;
//	}
//	friend Jet sqr(Jet const &B) {
//		Jet A;
//		for (IndexType N = 0; abs(N) < O; N++) {
//			for (IndexType M = 0; abs(M) < abs(N); M++) {
//				auto const NpM = N + M;
//				if (abs(NpM) < O) {
//					A[NpM] += T(2) * B[N] * B[M] * binco(N + M, N);
//				}
//			}
//			auto const NpN = 2 * N;
//			if (abs(NpN) < O) {
//				A[NpN] += sqr(B[N]) * factorial(2 * N) / sqr(factorial(N));
//			}
//		}
//		return A;
//	}
//	friend Jet operator*(Jet const &B, Jet const &C) {
//		Jet A;
//		for (IndexType N = 0; abs(N) < O; N++) {
//			for (IndexType M = 0; abs(M) < O; M++) {
//				auto const NpM = N + M;
//				if (abs(NpM) < O) {
//					A[NpM] += B[N] * C[M] * binco(N + M, N);
//				}
//			}
//		}
//		return A;
//	}
//	static constexpr Jet zeroJet() {
//		Jet Xn;
//		return Xn;
//	}
//	static constexpr Jet oneJet() {
//		Jet Xn;
//		IndexType I;
//		I[0] = 0;
//		Xn[I] = 1;
//		return Xn;
//	}
//	friend Jet operator/(Jet const &G, Jet const &H) {
//		return G * (1.0 / H);
//	}
//	friend Jet operator/(T const &G, Jet const &H) {
//		auto const inverse = [](Jet const &X) {
//			Jet Y;
//			auto const invX = 1.0 / X[0];
//			Y[0] = invX;
//			for (IndexType K = 0; abs(K) < O; K++) {
//				for (IndexType N = 0; abs(N) < abs(K); N++) {
//					Y[K] -= Y[N] * X[N - K] * binco(K, N) * invX;
//				}
//			}
//			return Y;
//		};
//		return G * inverse(H);
//	}
//	friend Jet operator*(Jet A, T const &v) {
//		for (auto &x : A.data_) {
//			x *= v;
//		}
//		return A;
//	}
//	friend Jet operator*(T const &v, Jet A) {
//		for (auto &x : A.data_) {
//			x *= v;
//		}
//		return A;
//	}
//	friend Jet operator/(Jet A, T const &v) {
//		for (auto &x : A.data_) {
//			x /= v;
//		}
//		return A;
//	}
//	Jet& operator+=(Jet const &C) {
//		*this = *this + C;
//		return *this;
//	}
//	Jet& operator-=(Jet const &C) {
//		*this = *this - C;
//		return *this;
//	}
//	friend Jet operator+(Jet const &B, Jet const &C) {
//		Jet A;
//		for (int i = 0; i < size_; i++) {
//			A.data_[i] = B.data_[i] + C.data_[i];
//		}
//		return A;
//	}
//	friend Jet operator-(Jet const &B, Jet const &C) {
//		Jet A;
//		for (int i = 0; i < size_; i++) {
//			A.data_[i] = B.data_[i] - C.data_[i];
//		}
//		return A;
//	}
//	Jet() {
//		data_.fill(0);
//	}
//	Jet(T value) {
//		data_.fill(0);
//		data_[0] = value;
//	}
//	friend Jet upShift(Jet const &A, int d = 0) {
//		Jet Deriv;
//		auto const I = IndexType::unit(d);
//		for (IndexType K = 0; abs(K) + 1 < O; K++) {
//			IndexType J = K + I;
//			Deriv[J] = A[K];
//		}
//		return Deriv;
//	}
//	friend Jet downShift(Jet const &A, int d = 0) {
//		Jet Deriv;
//		auto const I = IndexType::unit(d);
//		for (IndexType K = 0; abs(K) + 1 < O; K++) {
//			IndexType J = K + I;
//			Deriv[K] = A[J];
//		}
//		return Deriv;
//	}
//	// Α,α,Β,β,Γ,γ,Δ,δ,Ε,ε,Ζ,ζ,Η,η,Θ,θ,ϑ,Ι,ι,Κ,κ,ϰ,Λ,λ,Μ,μ,Ν,ν,Ξ,ξ,Ο,ο,Π,π,ϖ,Ρ,ρ,ϱ,Σ,σ,ς,Τ,τ,Υ,υ,Φ,φ,ϕ,Χ,χ,Ψ,ψ,Ω,ω
//	template<typename Function>
//	static Jet composite(Function const &f, Jet const &G) {
//		Jet H;
//		Jet<T, O + 1, 1> F;
//		std::array<Jet, O + 1> Q;
//		T const g = G[0];
//		for (int n = 0; n <= O; n++) {
//			F[n] = f(g, n);
//		}
//		Q[0][0] = 1;
//		auto dG = G;
//		dG[0] = 0;
//		Q[1] = dG;
//		for (int n = 1; n <= O; n++) {
//			Q[n + 1] = upShift(dG * Q[n]);
//			dG = upShift(dG);
//		}
//		for (int n = 0; n <= O; n++) {
//			std::cout << Q[n] << "\n";
//			H += Q[n] * F[n];
//		}
//		std::cout << "---------\n";
//		std::cout << H << "\n";
//		std::cout << "---------\n";
//		return H;
//	}
//	template<typename Function>
//	static Jet faàDiBruno(Function const &F, Jet const &G) {
//		Jet Q;
//		auto const fH = PowerSeries<T, Jet, 2 * O>(F);
//		Q = fH(G);
//		return Q;
//	}
//	friend Jet exp(Jet const &g) {
//		return composite([](T const &x, int n) {
//			return std::exp(x);
//		}, g);
//	}
//	static Jet genVar(T const &value, int dim = 0) {
//		Jet X;
//		X.data_[0] = value;
//		IndexType I;
//		I[dim] = 1;
//		X.data_[int(I)] = 1;
//		return X;
//	}
//	friend std::ostream& operator<<(std::ostream &os, Jet const &A) {
//		std::string str;
//		for (int α = 0; α < O; α++) {
//			os << print2string("   α={%i}-> ", α);
//			bool first = true;
//			for (IndexType η = 0; abs(η) < O; η++) {
//				if (abs(η) != α) {
//					continue;
//				}
//				if (!first) {
//					os << ", ";
//				}
//				os << print2string("{(");
//				for (int d = 0; d < D; d++) {
//					os << print2string("%i", η[d]);
//					if (d + 1 < D) {
//						os << print2string(",");
//					}
//				}
//				os << print2string("), %10.3e}", A[η]);
//				first = false;
//			}
//		}
//		os << "\n";
//		return os;
//	}
//private:
//	static constexpr int size_ = binco(O + D - 1, D);
//	std::array<T, size_> data_;
//}
//;

