#include <cmath>
#include <stdio.h>
#include <vector>
#include <array>
#include <cassert>
#include <astrotiger/chemistry.hpp>

#include <functional>
#include <limits>

#define KperEv (11604.525)

#define nH  0
#define nHP  1
#define nHE 2
#define nHEP   3
#define nHEPP  4
#define nE  5
#define NS  6
#define NF NS

const double arad = 7.5646e-15;
const auto tiny = 1.0e+3 * std::numeric_limits<double>::min();
const auto dhuge = std::numeric_limits<double>::max() / 1.0e+3;

using real = double;

template<int N>
std::array<std::array<real, N>, N> matrix_inverse(std::array<std::array<real, N>, N> A) {
	std::array<std::array<real, N>, N> Ainv;
	for (int n = 0; n < N; n++) {
		for (int m = 0; m < N; m++) {
			Ainv[n][m] = n == m ? 1.0 : 0.0;
			if (std::abs(A[n][m]) < 1.0e-50) {
				A[n][m] = 0.0;
			}
		}
	}
	for (int q = 0; q < N; q++) {
		for (int n = q; n < N; n++) {
			if (A[n][q] != 0.0) {
				std::swap(A[n], A[q]);
				std::swap(Ainv[n], Ainv[q]);
				break;
			}
		}
		for (int n = q; n < N; n++) {
			if (A[n][q] != 0.0) {
				const auto inv = 1.0 / A[n][q];
				for (int k = 0; k < q; k++) {
					A[n][k] *= inv;
				}
				A[n][q] = 1.0;
				for (int k = q + 1; k < N; k++) {
					A[n][k] *= inv;
				}
				for (int k = 0; k < N; k++) {
					Ainv[n][k] *= inv;
				}
			}
		}
		for (int n = q + 1; n < N; n++) {
			if (A[n][q] != 0.0) {
				A[n][q] -= A[q][q];
				for (int l = q + 1; l < N; l++) {
					A[n][l] -= A[q][l];
				}
				for (int l = 0; l < N; l++) {
					Ainv[n][l] -= Ainv[q][l];
				}
			}
		}
	}
	for (int q = N - 1; q >= 0; q--) {
		for (int n = q - 1; n >= 0; n--) {
			for (int k = 0; k < N; k++) {
				Ainv[n][k] -= Ainv[q][k] * A[n][q];
			}
			A[n][q] = 0.0;
		}
	}
	return Ainv;
}

void chemical_rates(double &k1, double &k2, double &k3, double &k4, double &k5, double &k6, double T) {

	const auto tev = T / KperEv;
	const auto logtev = std::log(tev);
	const auto logtev2 = logtev * logtev;
	const auto logtev3 = logtev2 * logtev;
	const auto logtev4 = logtev2 * logtev2;
	const auto logtev5 = logtev2 * logtev3;
	const auto logtev6 = logtev3 * logtev3;
	const auto logtev7 = logtev3 * logtev4;
	const auto logtev8 = logtev4 * logtev4;
	const auto logtev9 = logtev4 * logtev5;
	k1 = exp(
			-32.71396786375 + 13.53655609057 * logtev - 5.739328757388 * logtev2 + 1.563154982022 * logtev3 - 0.2877056004391 * logtev4
					+ 0.03482559773736999 * logtev5 - 0.00263197617559 * logtev6 + 0.0001119543953861 * logtev7 - 2.039149852002e-6 * logtev8);

	k2 = exp(
			-28.61303380689232 - 0.7241125657826851 * logtev - 0.02026044731984691 * logtev2 - 0.002380861877349834 * logtev3 - 0.0003212605213188796 * logtev4
					- 0.00001421502914054107 * logtev5 + 4.989108920299513e-6 * logtev6 + 5.755614137575758e-7 * logtev7 - 1.856767039775261e-8 * logtev8
					- 3.071135243196595e-9 * logtev9);

	k3 = exp(
			-44.09864886561001 + 23.91596563469 * logtev - 10.75323019821 * logtev2 + 3.058038757198 * logtev3 - 0.5685118909884001 * logtev4
					+ 0.06795391233790001 * logtev5 - 0.005009056101857001 * logtev6 + 0.0002067236157507 * logtev7 - 3.649161410833e-6 * logtev8);

	k4 = 3.92e-13 * std::pow(tev, -0.6353);

	if (tev > 0.1) {
		k4 += 1.54e-9 * (1. + 0.3 * std::exp(-8.099328789667 / tev)) / (exp(40.49664394833662 / tev) * std::pow(tev, 1.5));
	}

	k5 = exp(
			-68.71040990212001 + 43.93347632635 * logtev - 18.48066993568 * logtev2 + 4.701626486759002 * logtev3 - 0.7692466334492 * logtev4
					+ 0.08113042097303 * logtev5 - 0.005324020628287001 * logtev6 + 0.0001975705312221 * logtev7 - 3.165581065665e-6 * logtev8);
	const auto Tp7 = std::pow(T, 0.7);
	k6 = 1.33764e-9 / (Tp7 * (1.0 + 6.309573e-5 * Tp7));

}

const double hplanck = 6.6261e-27;
const double clight = 2.99792458e10;
const double kb = 1.3807e-16;
const double ergtoev = 6.242e+11;
const double evtoerg = 1.0 / ergtoev;

double homo_rate(double z, double alpha, double beta, double gamma, double z0, double z1) {
	auto c0 = beta * std::pow(z - z0, 2) / (1.0 + gamma * std::pow(z + z1, 2));
	c0 = std::min(c0, 100.0);
	c0 = std::max(c0, -100.0);
	return std::pow(1 + z, alpha) * std::exp(c0);
}

void radiation_rates(double &I20, double &I21, double &I22, double z) {
	I20 = 1.04e-12 * homo_rate(z, 0.231, -0.6818, 0.1646, 1.855, 0.3097);
	I21 = 1.84e-14 * homo_rate(z, -1.038, -1.1640, 0.1940, 1.973, -0.6561);
	I22 = 5.79e-13 * homo_rate(z, 0.278, -0.8260, 0.1730, 1.973, 0.2880);
}

void heating_rates(double &I20, double &I21, double &I22, double z) {
	I20 = 8.86e-25 * homo_rate(z, 0.231, -0.6818, 0.1646, 1.855, 0.3097);
	I21 = 5.86e-24 * homo_rate(z, -1.038, -1.1640, 0.1940, 1.973, -0.6561);
	I22 = 2.17e-25 * homo_rate(z, 0.278, -0.8260, 0.1730, 1.973, 0.2880);
}

double species_T(species s0, double egas0) {
	return egas0 / ((1.5 * s0.H + 3 * s0.Hp + 1.5 * s0.He + 3 * s0.Hep + 4.5 * s0.Hepp) * kb);
}

double species_T(std::array<double, NS> U, double egas0) {
	return egas0 / ((1.5 * U[nH] + 3 * U[nHP] + 1.5 * U[nHE] + 3 * U[nHEP] + 4.5 * U[nHEPP]) * kb);
}
double species_energy(species s0, double T) {
	return ((1.5 * s0.H + 3 * s0.Hp + 1.5 * s0.He + 3 * s0.Hep + 4.5 * s0.Hepp) * kb) * T;
}

double ion_energy(species s) {
	const auto evtoerg = 1.60218e-12;
	const auto Hion = -13.6 * evtoerg;
	const auto Heion = -(24.6 + 13.6 * 4) * evtoerg;
	const auto Hepion = -24.6 * evtoerg;
	return s.H * Hion + s.He * Heion + s.Hep * Hepion;
}

void cooling_rate(species s, double &dEdt, double &dEdtdT, double z, double T) {
	double C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13, dC1dT, dC2dT, dC3dT, dC4dT, dC5dT, dC6dT, dC7dT, dC8dT, dC9dT, dC10dT, dC11dT, dC12dT,
			dC13dT;

	const auto ne = s.Hp + s.Hep + 2 * s.Hepp;
	const auto T3 = T / 1e3;
	const auto T5 = T / 1e5;
	const auto T6 = T / 1e6;
	const auto tev = T / KperEv;
	const auto logtev = std::log(T);
	const auto tiny = std::numeric_limits<double>::min();

	const auto k1 = std::exp(
			-32.71396786375 + 13.53655609057 * logtev - 5.739328757388 * std::pow(logtev, 2) + 1.563154982022 * std::pow(logtev, 3)
					- 0.2877056004391 * std::pow(logtev, 4) + 0.03482559773736999 * std::pow(logtev, 5) - 0.00263197617559 * std::pow(logtev, 6)
					+ 0.0001119543953861 * std::pow(logtev, 7) - 2.039149852002e-6 * std::pow(logtev, 8));

	const auto dk1dT = (1 * 13.53655609057 - 2 * 5.739328757388 * logtev + 3 * 1.563154982022 * std::pow(logtev, 2) - 4 * 0.2877056004391 * std::pow(logtev, 3)
			+ 5 * 0.03482559773736999 * std::pow(logtev, 4) - 6 * 0.00263197617559 * std::pow(logtev, 5) + 7 * 0.0001119543953861 * std::pow(logtev, 6)
			- 8 * 2.039149852002e-6 * std::pow(logtev, 7)) * k1 / tev / KperEv;

	const auto k3 = std::exp(
			-44.09864886561001 + 23.91596563469 * logtev - 10.75323019821 * std::pow(logtev, 2) + 3.058038757198 * std::pow(logtev, 3)
					- 0.5685118909884001 * std::pow(logtev, 4) + 0.06795391233790001 * std::pow(logtev, 5) - 0.005009056101857001 * std::pow(logtev, 6)
					+ 0.0002067236157507 * std::pow(logtev, 7) - 3.649161410833e-6 * std::pow(logtev, 8));

	const auto dk3dT = (1 * 23.91596563469 - 2 * 10.75323019821 * logtev + 3 * 3.058038757198 * std::pow(logtev, 2)
			- 4 * 0.5685118909884001 * std::pow(logtev, 3) + 5 * 0.06795391233790001 * std::pow(logtev, 4) - 6 * 0.005009056101857001 * std::pow(logtev, 5)
			+ 7 * 0.0002067236157507 * std::pow(logtev, 6) - 8 * 3.649161410833e-6 * std::pow(logtev, 7)) * k3 / tev / KperEv;

	const auto k5 = std::exp(
			-68.71040990212001 + 43.93347632635 * logtev - 18.48066993568 * std::pow(logtev, 2) + 4.701626486759002 * std::pow(logtev, 3)
					- 0.7692466334492 * std::pow(logtev, 4) + 0.08113042097303 * std::pow(logtev, 5) - 0.005324020628287001 * std::pow(logtev, 6)
					+ 0.0001975705312221 * std::pow(logtev, 7) - 3.165581065665e-6 * std::pow(logtev, 8));

	const auto dk5dT = (1 * 43.93347632635 - 2 * 18.48066993568 * logtev + 3 * 4.701626486759002 * std::pow(logtev, 2)
			- 4 * 0.7692466334492 * std::pow(logtev, 3) + 5 * 0.08113042097303 * std::pow(logtev, 4) - 6 * 0.005324020628287001 * std::pow(logtev, 5)
			+ 7 * 0.0001975705312221 * std::pow(logtev, 6) - 8 * 3.165581065665e-6 * std::pow(logtev, 7)) * k5 / tev / KperEv;

	if (T > 250) {
		C1 = 7.50e-19 * std::exp(-118348 / T) / (1.0 + std::sqrt(T / 100000));
		dC1dT = 7.50e-19 * (-500 * (-236696000 + std::sqrt(10) * (-236696 + T) * std::sqrt(T)))
				/ (std::exp(118348 / T) * std::pow(1000 + std::sqrt(10) * std::sqrt(T), 2) * std::pow(T, 2));
	} else {
		C1 = dC1dT = 0.0;
	}
	if (T > 25) {
		C2 = 9.10e-27 * std::pow(T, -.1687) * std::exp(-13179 / T) / (1.0 + std::sqrt(T / 100000));
		dC2dT = 9.10e-27
				* (1.3179e10 * std::pow(T, 1.8374000000000001) + 4.167565728335908e7 * std::pow(T, 2.3373999999999997)
						- 168700.00000000003 * std::pow(T, 2.8373999999999997) - 2114.6150713545953 * std::pow(T, 3.3373999999999997))
				/ (std::exp(13179 / T) * std::pow(1000. + 3.1622776601683795 * std::sqrt(T), 2) * std::pow(T, 4.0061));
	} else {
		C2 = dC2dT = 0.0;
	}
	if (T > 1000) {
		C3 = 5.54e-17 * std::pow(T, -.397) * std::exp(-473638 / T) / (1.0 + std::sqrt(T / 100000));
		dC3dT = 5.54e-17
				* (4.73638e11 * std::pow(T, 2.294) + 1.4977748664068308e9 * std::pow(T, 2.794) - 397000. * std::pow(T, 3.2940000000000005)
						- 2836.5630611710367 * std::pow(T, 3.7940000000000005))
				/ (std::exp(473638 / T) * std::pow(1000. + 3.1622776601683795 * std::sqrt(T), 2) * std::pow(T, 4.691000000000001));
	} else {
		C3 = dC3dT = 0.0;
	}
	C4 = 2.18e-11 * k1;
	dC4dT = 2.18e-11 * dk1dT;
	C5 = 3.94e-11 * k3;
	dC5dT = 3.94e-11 * dk3dT;
	C6 = 8.72e-11 * k5;
	dC6dT = 8.72e-11 * dk5dT;
	if (T > 100) {
		C7 = 5.01e-27 * std::pow(T, -.1687) * std::exp(-55338 / T) / (1.0 + std::sqrt(T / 100000));
		dC7dT = 5.01e-27
				* (5.5338e10 * std::pow(T, 1.8374000000000001) + 1.749941211583978e8 * std::pow(T, 2.3373999999999997)
						- 168700. * std::pow(T, 2.8373999999999997) - 2114.6150713545953 * std::pow(T, 3.3373999999999997))
				/ (std::exp(55338 / T) * std::pow(1000. + 3.1622776601683795 * std::sqrt(T), 2) * std::pow(T, 4.0061));
	} else {
		C7 = dC7dT = 0.0;
	}
	C8 = 8.70e-27 * std::sqrt(T) * std::pow(T3, -0.2) / (1.0 + std::pow(T6, 0.7));
	dC8dT = 8.70e-27
			* -((25238.293779207717 - 2.9999999999999964e8 / std::pow(T, 0.7))
					/ (2.5118864315095773e8 + 31697.863849222253 * std::pow(T, 0.7) + std::pow(T, 1.4)));
	C9 = 1.55e-26 * std::pow(T, 0.3647);
	dC9dT = 0.3647 * C9 / T;
	if (T > 1000) {
		C10 = 1.24e-13 * std::pow(T, -1.5) * (1 + 0.3 * std::exp(-94000 / T)) * std::exp(-470000 / T);
		dC10dT = 1.24e-13
				* (((169200 + 470000 * std::exp(94000 / T)) * std::pow(T, 2.5) + (-0.45 - 1.5 * std::exp(94000 / T)) * std::pow(T, 3.5))
						/ (std::exp(564000 / T) * std::pow(T, 6.)));
	} else {
		C10 = dC10dT = 0.0;
	}
	C11 = 3.48e-26 * std::sqrt(T) * std::pow(T / 1000, -0.2) / (1.0 + std::pow(T / 1000000, 0.7));
	dC11dT = 3.48e-26
			* (-((25238.293779207717 - 3e8 / std::pow(T, 0.7)) / (2.5118864315095773e8 + 31697.863849222253 * std::pow(T, 0.7) + 1. * std::pow(T, 1.4))));
	C12 = 1.43e-27 * std::sqrt(T) * (1.1 + 0.34 * std::exp(-std::pow(5.50 - std::log10(T), 2) / 3.0));
	dC12dT = 1.43e-27
			* (0.55 / std::pow(T, 0.5)
					+ (std::pow(T, 1.0924131003119228) * (0.00002971599807934954 - 1.7857483385425343e-6 * std::log(T)))
							/ std::exp(0.06287056567053795 * std::pow(std::log(T), 2)));
	C13 = 5.64e-36 * std::pow(1 + z, 4) * (T - 2.73 * (1 + z));
	dC13dT = 5.64e-36 * std::pow(1 + z, 4);

	double J20, J21, J22;
	heating_rates(J20, J21, J22, z);

//	C1 = 0.0;
//	C2 = 0.0;
//	C3 = 0.0;
//	C4 = 0.0;
//	C5 = 0.0;
//	C6 = 0.0;
//	C7 = 0.0;
////	C8 = 0.0;
//	C9 = 0.0;
//	C10 = 0.0;
//	C11 = 0.0;
//	C12 = 0.0;
//	C13 = 0.0;
//	dC1dT = 0.0;
//	dC2dT = 0.0;
//	dC3dT = 0.0;
//	dC4dT = 0.0;
//	dC5dT = 0.0;
//	dC6dT = 0.0;
//	dC7dT = 0.0;
////	dC8dT = 0.0;
//	dC9dT = 0.0;
//	dC10dT = 0.0;
//	dC11dT = 0.0;
//	dC12dT = 0.0;
//	dC13dT = 0.0;
//	J20 = J21 = J22 = 0.0;

	dEdt = -ne
			* (C1 * s.H + C2 * ne * s.He + C3 * s.Hep + C4 * s.H + C5 * s.He + C6 * s.Hep + C7 * ne * s.Hep + C8 * s.Hp + C9 * s.Hep + C10 * s.Hep
					+ C11 * s.Hepp + C12 * (s.Hp + s.Hep + s.Hepp) + C13);
	dEdtdT = -ne
			* (dC1dT * s.H + dC2dT * ne * s.He + dC3dT * s.Hep + dC4dT * s.H + dC5dT * s.He + dC6dT * s.Hep + dC7dT * ne * s.Hep + dC8dT * s.Hp + dC9dT * s.Hep
					+ dC10dT * s.Hep + dC11dT * s.Hepp + dC12dT * (s.Hp + s.Hep + s.Hepp) + dC13dT);

	dEdt += J20 * s.H + J21 * s.He + J22 * s.Hep;
//	printf( "%e\n", dEdt);
}

double compute_next_energy(species s, double egas0, double z, double tmax) {
//
//	const auto T0 = species_T(s, egas0);
//	const auto ne = s.Hp + s.Hep + 2 * s.Hepp;
//	const auto cv = 1.5 * kb * (s.H + 2 * s.Hp + s.He + 2 * s.Hep + 3 * s.Hepp + ne);
//	double T = T0;
//	double egas = egas0;
//	double err;
//	double last_dT = 0.0;
//	double dT = 0.0;
//	double dEdt, dEdtdT;
//	constexpr double toler = 1e-12;
//	double t = 0.0;
//	while (t < tmax) {
//		//	printf( "!\n");
//		cooling_rate(s, dEdt, dEdtdT, z, T);
//		if (std::abs(dEdt) > 0.0) {
//			const auto dt_lim = std::abs(egas0 / dEdt * 0.1);
//			const auto dt = std::min(dt_lim, tmax - t);
//			do {
//				const auto f = egas - egas0 - dt * dEdt;
//				const auto dfdT = cv - dt * dEdtdT;
//				dT = -f / dfdT;
//				dT = std::min(std::max(dT, -0.1 * T), 0.1 * T);
//				T += dT;
//				egas = cv * T;
//				err = std::abs(dT / T);
//				if (err > toler) {
//					cooling_rate(s, dEdt, dEdtdT, z, T);
//				}
//
////			printf("---%e %e %e %e %e %e\n", T0, T, dT, f, err);
//			} while (err > toler);
//	//		printf("%e %e %e %e %e %e\n", t / tmax, (t + dt) / tmax, dt / tmax, T, egas0, dEdt);
//			t += dt;
//			egas0 = egas;
//		} else {
//			t = tmax;
//		}
//	}
//	return egas;
//
}

double compute_next_ne(double ne, species s0, species &s, double e0, double &e, double z, double dt) {
	e = e0;
	double K1, K2, K3, K4, K5, K6;
	double J20, J21, J22;
	const auto hinv = 1.0 / (1 + dt * (J20 + (K1 + K2) * ne));
	const auto hsum = (s0.H + s0.Hp);
	const auto hesum = (s0.He + s0.Hep + s0.Hepp);
	const auto cv = (ne + s0.H + s0.Hp + s0.He + s0.Hep + s0.Hepp) * 1.5 * kb;
	const auto T = e0 / cv;
	chemical_rates(K1, K2, K3, K4, K5, K6, T);
	radiation_rates(J20, J21, J20, z);
	s.H = (s0.H + dt * hsum * K2 * ne) * hinv;
	s.Hp = hsum - s.H;
	const auto heinv = 1.0
			/ (1 + dt * (J21 + J22 + dt * J22 * K3 * ne + dt * J21 * (J22 + (K5 + K6) * ne) + ne * (K3 + K4 + K5 + K6 + dt * (K3 * K5 + (K3 + K4) * K6) * ne)));
	s.He = (s0.He + dt * K4 * ne * (s0.Hep + dt * (s0.Hep + s0.Hepp) * K6 * ne) + dt * s0.He * (J22 + ne * (K4 + K5 + K6 + dt * K4 * K6 * ne))) * heinv;
	s.Hep = (s0.Hep + dt * s0.He * J21 + dt * s0.Hep * J21 + dt * ((s0.He + s0.Hep) * K3 + (s0.Hep + s0.Hepp + dt * hesum * J21) * K6) * ne
			+ dt * dt * hesum * K3 * K6 * ne * ne) * heinv;
	s.Hepp = hesum - s.He - s.Hep;
	return s.Hp + s.Hep + 2 * s.Hepp;
}

double chemistry_update(species s0, species &s, double egas0, double &egas, double z, double dt) {
	double nemin = 0.0;
	double nemax = (s0.H + s0.Hp) + 2 * (s0.He + s0.Hep + s0.Hepp);
	egas = egas0;
	double err;
	do {
		const auto nemid = 0.5 * (nemax + nemin);
		const auto f1 = nemax - compute_next_ne(nemax, s0, s, egas0, egas, z, dt);
		const auto f2 = nemid - compute_next_ne(nemid, s0, s, egas0, egas, z, dt);
		if( f1 * f2 < 0.0 ) {
			nemin = nemid;
		} else {
			nemax = nemid;
		}
		err = 1.0 - nemin / nemax;
	} while (err > 1.0e-7);
	return species_T(s,egas);
}

void chemistry_test() {
	species s0;
	const auto n = 1.0;
	s0.H = 0.91 * n;
	s0.Hp = 0.01 * n;
	s0.He = 0.08 * n;
	s0.Hep = 0.005 * n;
	s0.Hepp = 0.005 * n;
	double T0 = 1e+7;
	auto s = s0;
	printf("%14s %14s %14s %14s %14s %14s %14s %14s %14s %14s\n", "t", "T0", "T", "na", "ne", "H", "Hp", "He", "Hep", "Hepp");
	for (double dt = 1e1; dt < 1e18; dt *= 10.0) {
		const auto egas0 = species_energy(s0, T0);
		double egas;
		auto T = chemistry_update(s0, s, egas0, egas, 0.0, dt);
		const auto ne = s.Hp + s.Hep + 2 * s.Hepp;
		const auto na = s.Hp + s.H + 4 * s.He + 4 * s.Hep + 4 * s.Hepp;
		printf("%14.4e %14.4e %14.4e %14.4e %14.4e %14.4e %14.4e %14.4e %14.4e %14.4e\n", dt / (3600 * 24 * 365), T0, T, na, ne, s.H, s.Hp, s.He, s.Hep,
				s.Hepp);
	}

}
