/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/
#include "ValArray.hpp"

#include <hpx/hpx_init.hpp>
#include "Hydrodynamics.hpp"
#include "LegendreP.hpp"
#include "Polynomial.hpp"
#include "Real.hpp"
#include "HydroGrid.hpp"
#include "Options.hpp"
#include "Radiation.hpp"
#include "Root.hpp"
#include "Vector.hpp"
#include <unordered_map>
#include <numeric>
#include <stack>
#include <type_traits>
#include <debug/vector>
#include <span>
using namespace std;
#include "TriangularArray.hpp"
#include "Real.hpp"
#include "Interval.hpp"
#include "Limiters.hpp"
#include "Quadrature.hpp"
#include "Sod.hpp"
#include "Zobrist.hpp"
#include <random>
#include <vector>

template<int M, int P>
void compute() {
	using namespace Hydrodynamics;
	using state_type = ConservedState<Real, 1>;
	using prim_type = PrimitiveState<Real, 1>;
	constexpr int BW = 2;
	constexpr int N = M + 2 * BW;
	constexpr int NF = state_type::NFields;
	Real const dx = Real(1) / Real(M);
	Real const half = Real(0.5);
	Real const zero = Real(0.0);
	Real const one = Real(1.0);
	std::vector<state_type> U(N);
	std::vector<state_type> S(N);
	std::vector<state_type> U0(N);
	std::vector<state_type> dUdt(N);
	std::vector<state_type> UL(N);
	std::vector<state_type> UR(N);
	std::vector<state_type> F(N);
	Real tMax = Real(.01);
	Real t = Real(0);
	Real dt;
	for (int n = 0; n < N; n++) {
		state_type v;
		if (n > N / 2) {
			v.D = Real(sod_init.rhor);
			v.S = zero;
			v.E = Real(sod_init.pr / (sod_init.gamma - 1.0));
			v.tau = Real(energy2Entropy(sod_init.rhor, sod_init.pr / (sod_init.gamma - 1.0)));
		} else {
			v.D = Real(sod_init.rhol);
			v.S = zero;
			v.E = Real(sod_init.pl / (sod_init.gamma - 1.0));
			v.tau = Real(energy2Entropy(sod_init.rhol, sod_init.pl / (sod_init.gamma - 1.0)));
		}
		U[n] = state_type(v);
	}
	int n = 0;
	printf("Starting...\n");
	while (t < tMax) {
		U0 = U;
		for (int n = BW - 1; n < BW + 1; n++) {
			const auto maxE = max(max(U0[n].E, U0[n + 1].E), U0[n - 1].E);
			U[n].dualEnergyUpdate(maxE);
		}
		U = U0;
		for (int n = BW; n < N - BW; n++) {

		}
		for (int rk = 0; rk < BW; rk++) {
			for (int n = 0; n < BW; n++) {
				U[n] = U[BW];
				U[N - n - 1] = U[N - BW - 1];
			}
			for (int n = 1; n < N - 1; n++) {
				state_type const up = U[n + 1];
				state_type const u0 = U[n];
				state_type const um = U[n - 1];
				for (int f = 0; f < NF; f++) {
					S[n][f] = minmod(u0[f] - um[f], up[f] - u0[f]);
				}
			}
			for (int n = 1; n < N - 1; n++) {
				UL[n + 1] = U[n] + half * S[n];
				UR[n] = U[n] + half * S[n];
			}
			Real lambdaMax = Real(0);
			for (int n = BW; n < N - BW + 1; n++) {
				auto const rc = riemannSolver(UL[n], UR[n]);
				F[n] = rc.flux;
				lambdaMax = std::max(lambdaMax, rc.signalSpeed);
			}
			if (rk == 0) {
				dt = Real(0.4) * dx / lambdaMax;
				dt = std::min(dt, tMax - t);
				printf("%i %e %e %e\n", n, t, dt, lambdaMax);
			}
			for (int n = BW; n < N - BW; n++) {
				dUdt[n] = -(F[n + 1] - F[n]) * (one / dx);
			}
			if (rk == 0) {
				for (int n = BW; n < N - BW; n++) {
					U[n] += dt * dUdt[n];
				}
			} else {
				for (int n = BW; n < N - BW; n++) {
					U[n] = half * (U0[n] + U[n] + dt * dUdt[n]);
				}
			}
		}
		t += dt;
		n++;
	}
	{
		FILE *fp = fopen("test.txt", "wt");
		double l1 = 0.0;
		constexpr sod_init_t sod_init = { 1.0, 0.125, 1.0, 0.1, 5.0 / 3.0 };
		for (int n = BW; n < N - BW; n++) {
			sod_state_t sod;
			prim_type V(U[n]);
			double const x = double((n - N / 2) + 0.5) * double(dx);
			exact_sod(&sod, &sod_init, x, t, dx);
			double const err = abs(double(U[n].D) - sod.rho);
			l1 += err * dx;
			fprintf(fp, "%e %e %e %e %e\n", x, U[n].D, sod.rho, U[n].tau, energy2Entropy(V.rho, V.rho * V.eps));
		}
		fclose(fp);
		printf("Error = %12.4e\n", l1);
	}
}

void testIntegers();
void testRadiation();

template<typename T, int N>
Math::SquareMatrix<T, N> interpolationCoefficients() {
	static constexpr T half = T(0.5), one = T(1);
	Math::SquareMatrix<T, N> C;
	int factorial = T(1);
	for (int k = 0; k < N; k++) {
		T x = -half * Real(N - 1);
		for (int n = 0; n < N; n++) {
			T const numP = integerPower(x + half, k + 1);
			T const numM = integerPower(x - half, k + 1);
			T const den = T(k + 1);
			C[n, k] = (numP - numM) / den;
			x = x + one;
		}
		factorial *= k + 1;
	}
	C = matrixInverse(C);
	return C;
}

template<typename T>
Polynomial<T> legendrePolynomial(int n) {
	Polynomial<T> Pn;
	if (n == 0) {
		Pn[0] = Real(1);
	} else {
		if (n % 2 == 0) {
			Pn[0] = Real(1);
			Pn[1] = Real(0);
		} else {
			Pn[0] = Real(0);
			Pn[1] = Real(1);
		}
		//(n - k) * (n + k + 1) = n(n + 1) - k(k+1)
		// (k + 1) * (k + 2) = k^2 + 3*k + 2;
		for (int k = 0; k <= n - 2; k++) {
			Pn[k + 2] = -Real((n - k) * (n + k + 1)) / Real((k + 1) * (k + 2)) * Pn[k];
		}
		Real norm = Real(0);
		for (int k = 0; k <= n; k++) {
			norm += Pn[k];
		}
		norm = Real(1) / norm;
		for (int k = 0; k <= n; k++) {
			Pn[k] *= norm;
		}
	}
	return Pn;
}

int hpx_main(int argc, char *argv[]) {
	static constexpr Real one = Real(1), zero = Real(0), half = Real(0.5);
	printf("BEGIN\n");
	processOptions(argc, argv);
	static constexpr int N = 3;
	auto C = interpolationCoefficients<Real, N>();
	SquareMatrix<Real, N> a;
	for (int n = 0; n < N; n++) {
		if (n % 2 == 0) {
			a[n, 0] = Real(1);
			a[n, 1] = Real(0);
		} else {
			a[n, 1] = Real(1);
			a[n, 0] = Real(0);
		}
		for (int k = 0; k <= n - 2; k++) {
			a[n, k + 2] = -Real((n - k) * (n + k + 1)) / Real((k + 1) * (k + 2)) * a[n, k];
		}
		for (int k = n + 1; k <= N; k++) {
			a[n, k] = Real(0);
		}
		Real norm = Real(0);
		for (int k = 0; k <= N; k++) {
			norm += a[n, k];
		}
		norm = Real(1) / norm;
		for (int k = 0; k <= N; k++) {
			a[n, k] *= norm;
			a[n, k] *= integerPower(Real(2), k);
		}
	}
	std::cout << toString(C);
	Real const f[N] = { one, -one -one, one };
	Polynomial<Real> P;
	Vector<Real, N> Q;
	for (int n = 0; n < N; n++) {
		P[n] = zero;
		for (int k = 0; k < N; k++) {
			P[n] += C[n, k] * f[k];
		}
	}
	auto const A = matrixInverse(a);
	for (int n = 0; n < N; n++) {
		Q[n] = zero;
		for (int k = 0; k < N; k++) {
			Q[n] += A[n, k] * P[k];
		}
	}
	for (int n = 0; n < N; n++) {
		Real const dx = one;
		Real const x = (Real(n) + half) - Real(N) * half;
		Real const xa = x - half * dx;
		Real const xb = x + half * dx;
		printf("%e %e\n", x, polynomialIntegrate(P, xa, xb));

	}
	Polynomial<Real> R(Real(0));
	for (int n = 0; n < N; n++) {
		auto const Pn = legendrePolynomial<Real>(n);
		for (int k = 0; k <= n; k++) {
			R[n] += P[n] * Pn[k];
		}
	}
	std::cout << Math::toString(P) << "\n";
	printf("\n");
	std::cout << Math::toString(R) << "\n";
//	InterpolationCoefficients<Real, 3> c;
	//	testRadiation();

	HydroGrid<double, NDIM, PORDER> hydroTest;
//	TriangularIndices<3, 6> test;
	/*	using namespace Math;
	 constexpr int Ndim = 3;
	 TriangularArray<Real, 3, 4> test;
	 Vector<int, Ndim> dims { 5, 5, 5 };
	 MultiArray<Real, Ndim> A(dims);
	 auto view1 = A.getDefaultView();
	 auto view2 = A.getDefaultView();
	 view1 = view2;
	 processOptions(argc, argv);
	 compute<1024, 2>();*/
	printf("END\n");
	return hpx::local::finalize();
}

namespace Math {
Real normalDistribution() {
	static std::default_random_engine gen(42);
	static std::normal_distribution d { 0.0, 1.0 };
	return Real(d(gen));
}

template<int N>
Vector<Real, N> leastSquares(std::vector<Real> y, std::vector<Real> x) {
	using namespace Math;
	Real const zero(0), one(1), two(2);
	int const M = y.size();
	std::vector<Vector<Real, N>> F(M);
	SquareMatrix<Real, N> J;
	Vector<Real, N> beta;
	for (int m = 0; m < M; m++) {
		Real const xm = x[m];
		Real xn;
		for (int n = 0; n < N; n++) {
			xn = (n == 0) ? one : (xn * xm);
			F[m][n] = -xn;
		}
	}
	for (int n = 0; n < N; n++) {
		for (int m = 0; m < N; m++) {
			J[n, m] = zero;
			for (int k = 0; k < M; k++) {
				J[n, m] += F[k][m] * F[k][n];
			}
		}
	}
	auto const iJ = matrixInverse(J);
	for (int n = 0; n < N; n++) {
		beta[n] = zero;
	}
	for (int m = 0; m < M; m++) {
		beta -= y[m] * toVector(iJ * toColumnVector(F[m]));
	}
	Real r2 = zero;
	for (int m = 0; m < M; m++) {
		Real const xm = x[m];
		Real xn;
		Real r = y[m];
		for (int n = 0; n < N; n++) {
			xn = (n == 0) ? one : (xn * xm);
			r -= beta[n] * xn;
		}
		r2 += abs(r);
	}
	printf("%e\n\n", r2 / M);

	return beta;
}

Real derivative(std::function<Real(Real)> const &f, Real x) {
	static Real const eps = pow(Real(std::numeric_limits<double>::epsilon()), Real(1.0 / 6.0));
	static Real const mindx = Real(std::numeric_limits<double>::min());
	Real const third(1.0 / 3.0), half(0.5), zero(0), one(1), four(4);
	Real const dx = max(eps * abs(x), mindx);
	Real const dxinv = one / dx;
	Real const dfdx1 = half * (f(x + dx) - f(x - dx)) * dxinv;
	Real const dfdx2 = (f(x + half * dx) - f(x - half * dx)) * dxinv;
	return (four * dfdx2 - dfdx1) * third;
}
}

int main(int argc, char *argv[]) {
	/*Real err(0);
	 for (Real x = Real(std::numeric_limits<double>::min()); x < Real(100.0001); x *= Real(2)) {
	 Real const y = exp(x) * expm1(x);
	 Real const dydx1 = exp(x) * expm1(x) + exp(x) * exp(x);
	 Real const dydx2 = derivative([](Real x) {
	 return exp(x) * expm1(x);
	 }, x);
	 err += abs(dydx2 / dydx1 - Real(1));
	 printf("%e %e %e %e\n", x, dydx1, dydx2, dydx2 / dydx1 - Real(1));
	 }
	 printf("%e\n", err);*/
	/*	constexpr int N = 4;
	 constexpr int M = 100000;
	 std::vector<Real> x(M), y(M);
	 for (int m = 0; m < M; m++) {
	 x[m] = (Real(m) / Real(M));
	 y[m] = exp(x[m]);
	 }
	 auto const beta = leastSquares<N>(y, x);
	 for (int n = 0; n < N; n += 1) {
	 printf("%i %e\n", n, Real(1) / beta[n]);
	 }*/

	/*	constexpr int N = 4;
	 using namespace Math;
	 SquareMatrix<double, N> L;
	 SquareMatrix<double, N> U;
	 SquareMatrix<double, N> A = std::array<std::array<double, N>, N>( { { { 6, 1, 1, 1 }, { 1, 4, 1, 1 },
	 { 1, 1, 1, 1 }, { -1, 5, -1, 5 } } });
	 auto const P = matrixLUDecomposition(A, L, U);

	 std::cout << to_string(L) << "\n";
	 std::cout << to_string(U) << "\n";
	 std::cout << to_string(A) << "\n";
	 std::cout << to_string(L * U) << "\n";
	 std::cout << to_string(L * U - P * A) << "\n";*/

	/*	auto constexpr N = 256;
	 printf( "closed Newton-Cotes\n");
	 auto const rules1 = closedNewtonCotesRules<N>();
	 printf( "open Newton-Cotes\n");
	 auto const rules2 = openNewtonCotesRules<N>();
	 printf( "Gauss-Legendre\n");
	 auto const rules3 = gaussLegendreRules<N>();
	 printf( "Gauss-Lobatto\n");
	 auto const rules4 = gaussLobattoRules<N>();
	 printf( "Clenshaw-Curtis\n");
	 auto const rules5 = clenshawCurtisRules<N>();
	 Real w[5] = { Real(0), Real(0), Real(0), Real(0) };
	 for (int n = 0; n < N; n++) {
	 w[0] += rules1.w[n];
	 w[1] += rules2.w[n];
	 w[2] += rules3.w[n];
	 w[3] += rules4.w[n];
	 w[4] += rules5.w[n];
	 printf("%i (%11.4e, %11.4e) (%11.4e, %11.4e) (%11.4e, %11.4e) (%11.4e, %11.4e) (%11.4e, %11.4e)\n", n,
	 rules1.x[n], rules1.w[n], rules2.x[n], rules2.w[n], rules3.x[n], rules3.w[n], rules4.x[n], rules4.w[n],
	 rules5.x[n], rules5.w[n]);
	 }
	 for (int n = 0; n < 5; n++) {
	 printf("%e ", w[n] - Real(2));
	 }*/
	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };
	cfg.push_back("hpx.stacks.small_size=524288");
	hpx::init_params init_params;
	init_params.cfg = std::move(cfg);
	auto rc = hpx::init(argc, argv, init_params);
	return rc;
}
