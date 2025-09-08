#pragma once

#include "Multidices.hpp"



template<typename T, typename U, int N>
struct PowerSeries {
	template<typename Function>
	PowerSeries(Function const &f) {
		for (int k = 0; k < N; k++) {
			Cn[k] = f(0.0, k) / factorial(k);
		}
	}
	PowerSeries(std::array<T, N> const &cc) {
		for (int k = 0; k < N; k++) {
			Cn[k] = cc[k] / factorial(k);
		}
	}
	U operator()(U x) const {
		U sum = U(Cn[0]);
		U xn = x;
		for (int k = 1; k < N - 1; k++) {
			sum += Cn[k] * xn;
			xn *= x;
		}
		sum += Cn.back() * xn;
		return sum;
	}
	template<typename, int, int>
	friend struct Jet;
private:
	std::array<T, N> Cn;
};
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

template<typename T, int O, int D>
struct Jet {
	using IndexType = Multidices<D>;
	template<typename U>
	PowerSeries<T, U, O> generatePowerSeries() {
		PowerSeries<T, U, O> P(data_.begin(), data_.end());
	}
	T operator[](IndexType const &I) const {
		return data_[int(I)];
	}
	T& operator[](int i) {
		return data_[i];
	}
	T operator[](int i) const {
		return data_[i];
	}
	T& operator[](IndexType const &I) {
		return data_[int(I)];
	}
	Jet& operator*=(Jet const &other) {
		*this = *this * other;
		return *this;
	}
	Jet& operator*=(T const &other) {
		*this = *this * other;
		return *this;
	}
	Jet& operator/=(T const &other) {
		*this = *this / other;
		return *this;
	}
	friend Jet sqr(Jet const &B) {
		Jet A;
		for (IndexType N = 0; abs(N) < O; N++) {
			for (IndexType M = 0; abs(M) < abs(N); M++) {
				auto const NpM = N + M;
				if (abs(NpM) < O) {
					A[NpM] += T(2) * B[N] * B[M] * binco(N + M, N);
				}
			}
			auto const NpN = 2 * N;
			if (abs(NpN) < O) {
				A[NpN] += sqr(B[N]) * factorial(2 * N) / sqr(factorial(N));
			}
		}
		return A;
	}
	friend Jet operator*(Jet const &B, Jet const &C) {
		Jet A;
		for (IndexType N = 0; abs(N) < O; N++) {
			for (IndexType M = 0; abs(M) < O; M++) {
				auto const NpM = N + M;
				if (abs(NpM) < O) {
					A[NpM] += B[N] * C[M] * binco(N + M, N);
				}
			}
		}
		return A;
	}
	static constexpr Jet zeroJet() {
		Jet Xn;
		return Xn;
	}
	static constexpr Jet oneJet() {
		Jet Xn;
		IndexType I;
		I[0] = 0;
		Xn[I] = 1;
		return Xn;
	}
	friend Jet operator/(Jet const &G, Jet const &H) {
		return G * (1.0 / H);
	}
	friend Jet operator/(T const &G, Jet const &H) {
		auto const inverse = [](Jet const &X) {
			Jet Y;
			auto const invX = 1.0 / X[0];
			Y[0] = invX;
			for (IndexType K = 0; abs(K) < O; K++) {
				for (IndexType N = 0; abs(N) < abs(K); N++) {
					Y[K] -= Y[N] * X[N - K] * binco(K, N) * invX;
				}
			}
			return Y;
		};
		return G * inverse(H);
	}
	friend Jet operator*(Jet A, T const &v) {
		for (auto &x : A.data_) {
			x *= v;
		}
		return A;
	}
	friend Jet operator*(T const &v, Jet A) {
		for (auto &x : A.data_) {
			x *= v;
		}
		return A;
	}
	friend Jet operator/(Jet A, T const &v) {
		for (auto &x : A.data_) {
			x /= v;
		}
		return A;
	}
	Jet& operator+=(Jet const &C) {
		*this = *this + C;
		return *this;
	}
	Jet& operator-=(Jet const &C) {
		*this = *this - C;
		return *this;
	}
	friend Jet operator+(Jet const &B, Jet const &C) {
		Jet A;
		for (int i = 0; i < size_; i++) {
			A.data_[i] = B.data_[i] + C.data_[i];
		}
		return A;
	}
	friend Jet operator-(Jet const &B, Jet const &C) {
		Jet A;
		for (int i = 0; i < size_; i++) {
			A.data_[i] = B.data_[i] - C.data_[i];
		}
		return A;
	}
	Jet() {
		data_.fill(0);
	}
	Jet(T value) {
		data_.fill(0);
		data_[0] = value;
	}
	friend Jet upShift(Jet const &A, int d = 0) {
		Jet Deriv;
		auto const I = IndexType::unit(d);
		for (IndexType K = 0; abs(K) + 1 < O; K++) {
			IndexType J = K + I;
			Deriv[J] = A[K];
		}
		return Deriv;
	}
	friend Jet downShift(Jet const &A, int d = 0) {
		Jet Deriv;
		auto const I = IndexType::unit(d);
		for (IndexType K = 0; abs(K) + 1 < O; K++) {
			IndexType J = K + I;
			Deriv[K] = A[J];
		}
		return Deriv;
	}
	// Α,α,Β,β,Γ,γ,Δ,δ,Ε,ε,Ζ,ζ,Η,η,Θ,θ,ϑ,Ι,ι,Κ,κ,ϰ,Λ,λ,Μ,μ,Ν,ν,Ξ,ξ,Ο,ο,Π,π,ϖ,Ρ,ρ,ϱ,Σ,σ,ς,Τ,τ,Υ,υ,Φ,φ,ϕ,Χ,χ,Ψ,ψ,Ω,ω
	template<typename Function>
	static Jet composite(Function const &f, Jet const &G) {
		Jet H;
		Jet<T, O + 1, 1> F;
		std::array<Jet, O + 1> Q;
		T const g = G[0];
		for (int n = 0; n <= O; n++) {
			F[n] = f(g, n);
		}
		Q[0][0] = 1;
		auto dG = G;
		dG[0] = 0;
		Q[1] = dG;
		for (int n = 1; n <= O; n++) {
			Q[n + 1] = upShift(dG * Q[n]);
			dG = upShift(dG);
		}
		for (int n = 0; n <= O; n++) {
			std::cout << Q[n] << "\n";
			H += Q[n] * F[n];
		}
		std::cout << "---------\n";
		std::cout << H << "\n";
		std::cout << "---------\n";
		return H;
	}
	template<typename Function>
	static Jet faàDiBruno(Function const &F, Jet const &G) {
		Jet Q;
		auto const fH = PowerSeries<T, Jet, 2 * O>(F);
		Q = fH(G);
		return Q;
	}
	friend Jet exp(Jet const &g) {
		return composite([](T const &x, int n) {
			return std::exp(x);
		}, g);
	}
	static Jet genVar(T const &value, int dim = 0) {
		Jet X;
		X.data_[0] = value;
		IndexType I;
		I[dim] = 1;
		X.data_[int(I)] = 1;
		return X;
	}
	friend std::ostream& operator<<(std::ostream &os, Jet const &A) {
		std::string str;
		for (int α = 0; α < O; α++) {
			os << print2string("   α={%i}-> ", α);
			bool first = true;
			for (IndexType η = 0; abs(η) < O; η++) {
				if (abs(η) != α) {
					continue;
				}
				if (!first) {
					os << ", ";
				}
				os << print2string("{(");
				for (int d = 0; d < D; d++) {
					os << print2string("%i", η[d]);
					if (d + 1 < D) {
						os << print2string(",");
					}
				}
				os << print2string("), %10.3e}", A[η]);
				first = false;
			}
		}
		os << "\n";
		return os;
	}
private:
	static constexpr int size_ = binco(O + D - 1, D);
	std::array<T, size_> data_;
}
;

