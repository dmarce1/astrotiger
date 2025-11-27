#pragma once

#include <array>
#include <cmath>
#include <limits>

#include "math.hpp"
#include "matrix.hpp"

template<typename T, int N>
struct FourierLegendreTables {
	static constexpr T zero = T(0), one = T(1), two = T(2);
	static constexpr T half = one / two;
	static constexpr T π = std::numbers::pi_v<T>;
	static constexpr T legendre(int n, T x) {
		T p0, p1, p2;
		p0 = zero;
		p1 = one;
		for (int i = 0; i < n; i++) {
			T const η = i;
			p2 = ((two * η + one) * x * p1 - η * p0) / (η + one);
			p0 = p1;
			p1 = p2;
		}
		return p1;
	}
	static constexpr T legendreDerivative(int n, T x) {
		using std::abs;
		T const η = n;
		if (n == 0) {
			return zero;
		} else if (abs(x) == one) {
			T const sgn = x * nonepow(n + 1);
			return sgn * half * η * (η + one);
		} else {
			return (x * legendre(n, x) - legendre(n - 1, x)) * η / (sqr(x) - one);
		}
	}
	template<int N0 = N>
	static constexpr SquareMatrix<T, N0> decompose(SquareMatrix<T, N0> A) {
		for (int i = 0; i < N0; i++) {
			for (int n = i + 1; n < N0; n++) {
				A(n, i) /= A(i, i);
			}
			for (int n = i + 1; n < N0; n++) {
				for (int m = i + 1; m < N0; m++) {
					A(n, m) -= A(n, i) * A(i, m);
				}
			}
		}
		return A;
	}
	static constexpr auto evenSplit(SquareMatrix<T, N> A) {
		constexpr int Ne = (N + 1) / 2;
		SquareMatrix<T, Ne> Ae;
		for (int n = 0; n < Ne; n++) {
			for (int k = 0; k < Ne; k++) {
				Ae(n, k) = A(2 * n, N - k - 1);
			}
		}
		return Ae;
	}
	static constexpr auto oddSplit(SquareMatrix<T, N> A) {
		constexpr int No = N / 2;
		SquareMatrix<T, No> Ao;
		for (int n = 0; n < No; n++) {
			for (int k = 0; k < No; k++) {
				Ao(n, k) = A(2 * n + 1, k);
			}
		}
		return Ao;
	}
	constexpr FourierLegendreTables() {
		using std::abs;
		using std::cos;
		Vector<T, N> X(zero);
		Vector<T, N> W(zero);
		Vector<T, N> norm(zero);
		SquareMatrix<T, N> Λ(zero);
		SquareMatrix<T, N> Π(zero);
		for (int km = 0, kp = N - 1; km <= kp; km++, kp--) {
			T xmin = cos((4.0 * km - 1.0) / (4.0 * N - 2.0) * π);
			T xmax = cos((4.0 * km + 1.0) / (4.0 * N - 2.0) * π);
			T xmid;
			do {
				xmid = half * (xmin + xmax);
				T const pmid = legendre(N, xmid);
				T const pmin = legendre(N, xmin);
				if (pmin * pmid < zero) {
					xmax = xmid;
				} else {
					xmin = xmid;
				}
			} while (xmax > xmin);
			T const x = abs(xmid);
			T const dpdx = legendreDerivative(N, x);
			W[km] = W[kp] = two / ((one - sqr(x)) * sqr(dpdx));
			X[km] = -x;
			X[kp] = +x;
		}
		for (int n = 0; n < N; n++) {
			norm[n] = two / (two * n + one);
			for (int k = 0; k < N; k++) {
				Λ(n, k) = legendre(n, X[k]);
			}
		}
		for (int km = 0, kp = N - 1; km < kp; km++, kp--) {
			Π(km, km) = +one;
			Π(km, kp) = +one;
			Π(kp, km) = -one;
			Π(kp, kp) = +one;
		}
		for (int k = N / 2; k != (N + 1) / 2; k++) {
			Π(k, k) = one;
		}
		Vector<int, N> p;
		int i = 0;
		for (int n = 0; n < N; n += 2) {
			p[i++] = n;
		}
		i = N;
		for (int n = 1; n < N; n += 2) {
			p[--i] = n;
		}
		auto const P = permutationMatrix<T, N>(p);
		Λs_ = Π * transpose(Λ) * inverse(P);
		Λa_ = P * inverse(transpose(Λ)) * inverse(Π);
		X_ = X;
		W_ = W;
		LUa_ = decompose(Λa_);
		LUs_ = decompose(Λs_);
	}
	SquareMatrix<T, N> LUa_;
	SquareMatrix<T, N> LUs_;
	SquareMatrix<T, N> Λa_;
	SquareMatrix<T, N> Λs_;
	Vector<T, N> X_;
	Vector<T, N> W_;
};

template<typename T, int P>
auto gaussLegendreQuadrature() {
	constexpr FourierLegendreTables<T, P> tables;
}

template<typename T, int D, int P>
auto fourierLegendre(Vector<T, pow<D>(P)> Z, int dir) {
	using std::fma;
	using std::swap;
	FourierLegendreTables<T, P> tables { };
	auto LU = dir > 0 ? tables.LUa_ : tables.LUs_;
	constexpr int N = pow<D - 1>(P);
	auto const transpose = [](int d1, int d2, auto x) {
		if (d1 > d2) {
			std::swap(d1, d2);
		} else {
			if (d1 == d2) {
				return;
			}
		}
		int const Na = pow(P, D - (d2 + 1));
		int const Nb = pow(P, d2 - (d1 + 1));
		int const Nc = pow(P, d1);
		for (int a = 0; a < Na; a++) {
			for (int b = 0; b < Nb; b++) {
				for (int c = 0; c < Nc; c++) {
					for (int p1 = 0; p1 < P; p1++) {
						for (int p2 = 0; p2 < p1; p2++) {
							int const n1 = c + Nc * (p1 + P * (b + Nb * (p2 + P * a)));
							int const n2 = c + Nc * (p2 + P * (b + Nb * (p1 + P * a)));
							swap(x[n1], x[n2]);
						}
					}
				}
			}
		}
	};
	for (int d = 0; d < D; d++) {
		transpose(D - 1, d, Z);
		if (dir > 0) {
			for (int km = 0, kp = N * (P - 1); km < kp; km += N, kp -= N) {
				for (int i = 0; i < N; i++) {
			//		Z[km + i] += Z[kp + i];
				}
				for (int i = 0; i < N; i++) {
			//		Z[kp + i] = 2 * Z[kp + i] - Z[km + i];
				}
			}
		}
		for (int n = 0; n < N * ((P + 1) / 2); n += N) {
			for (int i = 0; i < N; i++) {
			//	Z[n + i] *= LU(n, n);
			}
			for (int i = 0; i < N; i++) {
				for (int k = n + N; k < N * ((P + 1) / 2); k += N) {
			//		Z[n + i] += LU(n, k) * Z[k + i];
				}
			}
		}
		for (int n = N * ((P + 1) / 2); n >= 0; n -= N) {
			for (int k = n - N; k >= 0; k -= N) {
				for (int i = 0; i < N; i++) {
				//	Z[n + i] += LU(n, k) * Z[k + i];
				}
			}
		}
		for (int n = N * ((P + 1) / 2); n < N * P; n += N) {
			for (int i = 0; i < N; i++) {
			//	Z[n + i] *= LU(n, n);
			}
			for (int k = n + N; k < N * P; k += N) {
				for (int i = 0; i < N; i++) {
				//	Z[n + i] += LU(n, k) * Z[k + i];
				}
			}
		}
		for (int n = N * (P - 1); n >= N * P / 2; n -= N) {
			for (int k = n - N; k >= N * (P / 2); k -= N) {
				for (int i = 0; i < N; i++) {
					//Z[n + i] += LU(n, k) * Z[k + i];
				}
			}
		}
		if (dir < 0) {
			for (int km = 0, kp = N * (P - 1); km < kp; km += N, kp -= N) {
				for (int i = 0; i < N; i++) {
				//	Z[km + i] = T(0.5) * (Z[km + i] - Z[kp + i]);
				}
				for (int i = 0; i < N; i++) {
				//	Z[kp + i] += Z[km + i];
				}
			}
		}
		transpose(D - 1, d, Z);
	}
	return Z;
}

//template<typename T, int D, int P>
//auto fourierLegendreTransform2(Vector<T, pow<D>(P)> U) {
//	using std::fma;
//	using std::swap;
//	constexpr FourierLegendreTables<T, P> tables { };
//	constexpr int P2D = pow<D>(P);
//	constexpr auto Λ = inverse(tables.Mnk_) * tables.Λnk_ * diagonal<T, P>(tables.Wk_);
//	constexpr auto LU = FourierLegendreTables<T, P>::decompose(Λ);
//	int Nlo = 1;
//	 int Nhi = pow(P, D - 1);
//	auto *X = &(U[0]);
//	for (int d = 0; d < D; d++) {
//		for (int base = 0; base < Nhi; base++) {
//			auto *Z = X + base * Nlo * P;
//			for (int n = 0; n < P; n++) {
//				for (int i = 0; i < Nlo; i++) {
//					Z[Nlo * n + i] = LU(n, n) * Z[Nlo * n + i];
//				}
//				for (int k = n + 1; k < P; k++) {
//					for (int i = 0; i < Nlo; i++) {
//						Z[Nlo * n + i] += LU(n, k) * Z[Nlo * k + i];
//					}
//				}
//			}
//			for (int n = P - 1; n >= 0; n--) {
//				for (int k = n - 1; k >= 0; k--) {
//					for (int i = 0; i < Nlo; i++) {
//						Z[Nlo * n + i] += LU(n, k) * Z[Nlo * k + i];
//					}
//				}
//			}
//		}
//		Nlo *= P;
//		Nhi /= P;
//	}
//	return U;
//}

//#include <algorithm>
//#include <cmath>
//#include <array>
//#include <limits>
//#include <numbers>
//#include <utility>
//
//#include "math.hpp"
//#include "matrix.hpp"
//
//template<int N, int D = 0>
//constexpr auto legendreP(auto x) {
//	using T = std::remove_cvref_t<decltype(x)>;
//	std::array<T, D + 1> Pn { }, Pnp1 { }, Pnm1 { };
//	Pn[0] = T(1);
//	for (int n = 0; n < N; n++) {
//		Pnp1[0] = (T(2 * n + 1) * (x * Pn[0]) - T(n) * Pnm1[0]) / T(n + 1);
//		for (int d = 1; d <= D; d++) {
//			Pnp1[d] = (T(2 * n + 1) * (Pn[d - 1] * T(d) + x * Pn[d]) - T(n) * Pnm1[d]) / T(n + 1);
//		}
//		for (int d = 0; d <= D; d++) {
//			Pnm1[d] = Pn[d];
//			Pn[d] = Pnp1[d];
//		}
//	}
//	return Pn;
//}
//
//enum class QuadratureType : int {
//	legendre, lobatto
//};
//
//template<typename T, QuadratureType Q, int N, int K = (Q == QuadratureType::lobatto) ? N : N + 1>
//struct QuadraturePoints {
//	std::array<T, K> x;
//	std::array<T, K> w;
//};
//
//template<typename T, int P>
//constexpr auto gaussLegendrePoints() {
//	using std::abs;
//	using std::atan;
//	using std::cos;
//	using std::sin;
//	using std::numeric_limits;
//	static constexpr T ε = 4 * numeric_limits < T > ::epsilon();
//	static constexpr T one = T(1);
//	static constexpr T two = T(2);
//	static constexpr T half = one / two;
//	static constexpr T π = std::numbers::pi_v<T>;
//	QuadraturePoints<T, QuadratureType::legendre, P> pts;
//	for (int n = 0; n < (P + 1) / 2; n++) {
//		T θ, dθ, w, x;
//		θ = π * (one - half * (T(2 * n + 1)) / T(P));
//		do {
//			x = cos(θ);
//			auto const Pn = legendreP<P, 1>(x);
//			w = one / (sin(θ) * Pn[1]);
//			dθ = w * Pn[0];
//			θ += dθ;
//		} while (abs(dθ / θ) > ε);
//		int const m = P - 1 - n;
//		w = two * sqr(w);
//		pts.x[n] = -x;
//		pts.w[n] = w;
//		pts.x[m] = +x;
//		pts.w[m] = w;
//	}
//	return pts;
//}
//
//template<typename T, int P>
//constexpr auto gaussLobattoPoints() {
//	using std::abs;
//	using std::atan;
//	using std::cos;
//	using std::sin;
//	using std::numeric_limits;
//	static constexpr T ε = 4 * numeric_limits < T > ::epsilon();
//	static constexpr T zero = T(0);
//	static constexpr T one = T(1);
//	static constexpr T two = T(2);
//	static constexpr T half = one / two;
//	static constexpr T π = std::numbers::pi_v<T>;
//	QuadraturePoints<T, QuadratureType::lobatto, P> pts;
//	T wtot = zero;
//	for (int n = 1; n < (P + 1) / 2; n++) {
//		T θ, dθ, x;
//		θ = π * (one - T(n) / T(P));
//		do {
//			x = cos(θ);
//			auto const Pn = legendreP<P, 2>(x);
//			dθ = Pn[1] / (sin(θ) * Pn[2]);
//			θ += dθ;
//		} while (abs(dθ) > ε);
//		int const m = P - 1 - n;
//		auto const Pn = legendreP<P, 2>(x);
//		auto const w = two / (P * (P + 1) * sqr(Pn[0]));
//		pts.x[n] = -x;
//		pts.w[n] = w;
//		pts.x[m] = +x;
//		pts.w[m] = w;
//		auto const α = (2 * n + 1 < P) ? two : one;
//		wtot += α * w;
//	}
//	T const we = one - half * wtot;
//	pts.x[0] = -one;
//	pts.w[0] = we;
//	pts.x[P] = +one;
//	pts.w[P] = we;
//	return pts;
//}
//
//template<typename T, int N, QuadratureType W, QuadratureType ...Ws>
//using QuadratureTuple = std::remove_cvref<decltype(std::tuple_cat(
//				std::tuple<QuadraturePoints<T, W, N>>(),
//				std::conditional_t<sizeof...(Ws) == 0, std::tuple<>, QuadratureTuple<T, N, Ws...>>()))>;
//
//template<typename T, int N>
//constexpr auto genQuadratureTuple() {
//	return std::tuple();
//}
//
//template<typename T, int N, QuadratureType Q, QuadratureType ...Qs>
//constexpr auto genQuadratureTuple() {
//	auto const qnp1 = genQuadratureTuple<T, N, Qs...>();
//	if constexpr (qTypes[dim] == QuadratureType::legendre) {
//		return std::tuple_cat(std::tuple(gaussLegendrePoints<T, N>()), qnp1);
//	} else if constexpr (qTypes[dim] == QuadratureType::lobatto) {
//		auto const qn = gaussLobattoPoints<T, N>();
//		return std::tuple_cat(std::tuple(gaussLobattoPoints<T, N>()), qnp1);
//	}
//}
//
//template<typename T, int N>
//constexpr auto genAnalysisTransforms() {
//	return std::tuple();
//}
//
//template<typename T, int N, QuadratureType Q, QuadratureType ...Qs>
//constexpr auto genAnalysisTransforms() {
//	constexpr int dim = sizeof...(Qs);
//	constexpr int K = qSizes[dim];
//	constexpr int Ko2 = (K + 1) / 2;
//	constexpr auto Q = std::get<dim>(qPoints);
//	Matrix<T, N, K> Λ(zero);
//	for (int k = 0; k < K; k++) {
//		T const x = Q.x[k];
//		T const w = Q.w[k];
//		T Pnm1 = zero;
//		T Pn = one;
//		for (int n = 0; n < N; n++) {
//			T const η = T(n);
//			if (((even(n) && (k < Ko2)) || (odd(n) && (k >= Ko2)))) {
//				Λ(n, k) = (η + half) * w * Pn;
//			}
//			if (n + 1 < N) {
//				T const Pnp1 = ((two * η + one) * x * Pn - η * Pnm1) / (η + one);
//				Pnm1 = Pn;
//				Pn = Pnp1;
//			}
//		}
//	}
//	return std::tuple_cat(std::tuple(Λ), genAnalysisTransforms<T, N, Qs...>());
//}
//
//template<typename T, int D, int N, QuadratureType Q, QuadratureType ...Qs>
//struct FourierLegendreTables {
//	std::array<QuadratureType, D> qTypes;
//	std::array<int, D> qSizes;
//	QuadratureTuple<T, N, Q, Qs...> qPoints;
////	auto qPoints = []() {
////	}();
////	template<int dim>
////	constexpr auto FourierLegendre<T, D, N, Qs...>::genQuadratureTuple() {
////		constexpr auto qLegendre = gaussLegendrePoints<T, N>();
////		constexpr auto qLobatto = gaussLobattoPoints<T, N>();
////		if constexpr (dim == D) {
////			return std::tuple();
////		} else {
////			constexpr auto qTypes = std::array<QuadratureType, D> { Qs... };
////			auto const qnp1 = genQuadratureTuple<dim + 1>();
////			if constexpr (qTypes[dim] == QuadratureType::legendre) {
////				return std::tuple_cat(std::tuple(qLegendre), qnp1);
////			} else if constexpr (qTypes[dim] == QuadratureType::lobatto) {
////				auto const qn = gaussLobattoPoints<T, N>();
////				return std::tuple_cat(std::tuple(qLobatto), qnp1);
////			}
////		}
////	}
//	constexpr FourierLegendreTables() {
//		constexpr T zero = T(0);
//		constexpr T one = T(1);
//		constexpr T two = T(2);
//		constexpr T half = one / two;
//		qTypes = std::array<QuadratureType, D> { { Q, Qs... } };
//		qPoints = genQuadratureTuple<T, N, Q, Qs...>();\
//		for (int d = 0; d < D; d++) {
//			qSizes[dim] = (qTypes[dim] == QuadratureType::lobatto) ? (N + 1) : N;
//		}
//	}
//};
//
//template<typename T, int D, int N, QuadratureType ...Qs>
//struct FourierLegendre {
//	constexpr FourierLegendre() {
//		auto const Λ = std::get<0>(aTransforms);
//		auto const lu = luDecomposition(Λ);
//		std::cout << Λ;
//		std::cout << lu.U;
//		std::cout << lu.L;
//		std::cout << lu.P;
//	}
//private:
//};
//
//template<typename T, int D, int N, QuadratureType ...Qs>
//template<int dim>
//constexpr auto FourierLegendre<T, D, N, Qs...>::genQuadratureTuple() {
//	constexpr auto qLegendre = gaussLegendrePoints<T, N>();
//	constexpr auto qLobatto = gaussLobattoPoints<T, N>();
//	if constexpr (dim == D) {
//		return std::tuple();
//	} else {
//		constexpr auto qTypes = std::array<QuadratureType, D> { Qs... };
//		auto const qnp1 = genQuadratureTuple<dim + 1>();
//		if constexpr (qTypes[dim] == QuadratureType::legendre) {
//			return std::tuple_cat(std::tuple(qLegendre), qnp1);
//		} else if constexpr (qTypes[dim] == QuadratureType::lobatto) {
//			auto const qn = gaussLobattoPoints<T, N>();
//			return std::tuple_cat(std::tuple(qLobatto), qnp1);
//		}
//	}
//}
//
//template<typename T, int D, int N, QuadratureType ...Qs>
//template<int dim>
//constexpr auto FourierLegendre<T, D, N, Qs...>::genQuadratureSizes() {
//	if constexpr (dim == D) {
//		return std::array<int, D>();
//	} else {
//		constexpr auto qTypes = std::array<QuadratureType, D> { Qs... };
//		auto qSizes = genQuadratureSizes<dim + 1>();
//		if constexpr (qTypes[dim] == QuadratureType::legendre) {
//			qSizes[dim] = N;
//		} else if constexpr (qTypes[dim] == QuadratureType::lobatto) {
//			qSizes[dim] = N + 1;
//		}
//		return qSizes;
//	}
//}
//
//template<typename T, int D, int N, QuadratureType ...Qs>
//template<int dim>
//constexpr auto FourierLegendre<T, D, N, Qs...>::genAnalysisTransforms() {
//}
//
