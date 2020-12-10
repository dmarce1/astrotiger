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

double ion_energy(species s) {
	const auto evtoerg = 1.60218e-12;
	const auto Hion = -13.6 * evtoerg;
	const auto Heion = -(24.6 + 13.6 * 4) * evtoerg;
	const auto Hepion = -24.6 * evtoerg;
	return s.H * Hion + s.He * Heion + s.Hep * Hepion;
}

void cooling_rate(double &C1, double &C2, double &C3, double &C4, double &C5, double &C6, double &C7, double &C8, double &C9, double &C10, double &C11,
		double &C12, double &C13, double &dC1dT, double &dC2dT, double &dC3dT, double &dC4dT, double &dC5dT, double &dC6dT, double &dC7dT, double &dC8dT,
		double &dC9dT, double &dC10dT, double &dC11dT, double &dC12dT, double &dC13dT, double T, double z) {
	const auto T3 = T / 1e3;
	const auto T5 = T / 1e5;
	const auto T6 = T / 1e6;
	const auto tev = T / KperEv;
	const auto logtev = std::log(T);
	const auto tiny = std::numeric_limits<double>::min();
//
//	const auto k1 = exp(
//			-32.71396786375 + 13.53655609057 * logtev - 5.739328757388 * std::pow(logtev, 2) + 1.563154982022 * std::pow(logtev, 3)
//					- 0.2877056004391 * std::pow(logtev, 4) + 0.03482559773736999 * std::pow(logtev, 5) - 0.00263197617559 * std::pow(logtev, 6)
//					+ 0.0001119543953861 * std::pow(logtev, 7) - 2.039149852002e-6 * std::pow(logtev, 8));
//
//	const auto dk1dT = (1 * 13.53655609057 - 2 * 5.739328757388 * logtev + 3 * 1.563154982022 * std::pow(logtev, 2) - 4 * 0.2877056004391 * std::pow(logtev, 3)
//			+ 5 * 0.03482559773736999 * std::pow(logtev, 4) - 6 * 0.00263197617559 * std::pow(logtev, 5) + 7 * 0.0001119543953861 * std::pow(logtev, 6)
//			- 8 * 2.039149852002e-6 * std::pow(logtev, 7)) * k1 / tev / KperEv;
//
//	const auto k3 = exp(
//			-44.09864886561001 + 23.91596563469 * logtev - 10.75323019821 * std::pow(logtev, 2) + 3.058038757198 * std::pow(logtev, 3)
//					- 0.5685118909884001 * std::pow(logtev, 4) + 0.06795391233790001 * std::pow(logtev, 5) - 0.005009056101857001 * std::pow(logtev, 6)
//					+ 0.0002067236157507 * std::pow(logtev, 7) - 3.649161410833e-6 * std::pow(logtev, 8));
//
//	const auto dk3dT = (1 * 23.91596563469 - 2 * 10.75323019821 * logtev + 3 * 3.058038757198 * std::pow(logtev, 2)
//			- 4 * 0.5685118909884001 * std::pow(logtev, 3) + 5 * 0.06795391233790001 * std::pow(logtev, 4) - 6 * 0.005009056101857001 * std::pow(logtev, 5)
//			+ 7 * 0.0002067236157507 * std::pow(logtev, 6) - 8 * 3.649161410833e-6 * std::pow(logtev, 7)) * k3 / tev / KperEv;
//
//	const auto k5 = exp(
//			-68.71040990212001 + 43.93347632635 * logtev - 18.48066993568 * std::pow(logtev, 2) + 4.701626486759002 * std::pow(logtev, 3)
//					- 0.7692466334492 * std::pow(logtev, 4) + 0.08113042097303 * std::pow(logtev, 5) - 0.005324020628287001 * std::pow(logtev, 6)
//					+ 0.0001975705312221 * std::pow(logtev, 7) - 3.165581065665e-6 * std::pow(logtev, 8));
//
//	const auto dk5dT = (1 * 43.93347632635 - 2 * 18.48066993568 * logtev + 3 * 4.701626486759002 * std::pow(logtev, 2)
//			- 4 * 0.7692466334492 * std::pow(logtev, 3) + 5 * 0.08113042097303 * std::pow(logtev, 4) - 6 * 0.005324020628287001 * std::pow(logtev, 5)
//			+ 7 * 0.0001975705312221 * std::pow(logtev, 6) - 8 * 3.165581065665e-6 * std::pow(logtev, 7)) * k5 / tev / KperEv;
//
//	C1 = 7.50e-19 * (1.0 + std::sqrt(T / 100000)) * std::exp(-118348 / T);
//	dC1dT = 7.50e-19 * ((118348. + (374.2492365256074 + 0.0015811388300841897 * T) * std::sqrt(T)) / (std::exp(118348 / T) * T * T));
//	C2 = 9.10e-27 * (1.0 + std::sqrt(T / 100000)) * std::pow(T, -.1687) * std::exp(-13179 / T);
//	dC2dT = 9.10e-27
//			* ((13179 * std::pow(T, 1.8374) + 41.67565728335908 * std::pow(T, 2.3374) - 0.1687 * std::pow(T, 2.8374)
//					+ 0.0010476625888137844 * std::pow(T, 3.3374)) / (std::exp(13179 / T) * std::pow(T, 4.0061)));
//	C3 = 5.54e-17 * (1.0 + std::sqrt(T / 100000)) * std::pow(T, -.397) * std::exp(-473638 / T);
//	dC3dT =
//			5.54e-17
//					* ((473638. * std::pow(T, 2.294) + 1497.774866406831 * std::pow(T, 2.794) - 0.397 * std::pow(T, 3.294)
//							+ 0.00032571459899734327 * std::pow(T, 3.794)) / (std::exp(473638 / T) * std::pow(T, 4.691)));
//	C4 = 2.18e-11 * k1;
//	dC4dT = 2.18e-11 * dk1dT;
//	C5 = 3.94e-11 * k3;
//	dC5dT = 3.94e-11 * dk3dT;
//	C6 = 8.72e-11 * k5;
//	dC6dT = 8.72e-11 * dk5dT;
//	C7 = 5.01e-27 * (1.0 + std::sqrt(T / 100000)) * std::pow(T, -.1687) * std::exp(-55338 / T);
//	dC7dT = 5.01e-27
//			* ((55338 * std::pow(T, 1.8374) + 174.9941211583978 * std::pow(T, 2.3374) - 0.1687 * std::pow(T, 2.8374)
//					+ 0.0010476625888137844 * std::pow(T, 3.3374)) / (std::exp(55338 / T) * std::pow(T, 4.0061)));
//	C8 = 8.70e-27 * std::sqrt(T) * std::pow(T3, -0.2) / (1.0 + std::pow(T6, 0.7));
//	dC8dT = 8.70e-27 * (-((25238.293779207717 - 3e8 / std::pow(T, 0.7)) / (2.5118864315095773e8 + 31697.863849222253 * std::pow(T, 0.7) + std::pow(T, 1.4))));
//	C9 = 1.55e-26 * std::pow(T, 0.3647);
//	dC9dT = 0.3647 * C9 / T;
//	C10 = 1.24e-13 * std::pow(T, -1.5) * (1 + 0.3 * std::exp(-94000 / T)) * std::exp(-470000 / T);
//	dC10dT = 1.24e-13
//			* (((169200 + 470000 * std::exp(94000 / T)) * std::pow(T, 2.5) + (-0.45 - 1.5 * std::exp(94000 / T)) * std::pow(T, 3.5))
//					/ (std::exp(564000 / T) * std::pow(T, 6.)));
//	C11 = 3.48e-26 * std::sqrt(T) * std::pow(T / 1000, -0.2) / (1.0 + std::pow(T / 1000000, 0.7));
//	dC11dT = 3.48e-26
//			* (-((25238.293779207717 - 3e8 / std::pow(T, 0.7)) / (2.5118864315095773e8 + 31697.863849222253 * std::pow(T, 0.7) + 1. * std::pow(T, 1.4))));
//	C12 = 1.43e-27 * std::sqrt(T) * (1.1 + 0.34 * std::exp(-std::pow(5.50 - std::log10(T), 2) / 3.0));
//	dC12dT = 1.43e-27
//			* (0.55 / std::pow(T, 0.5)
//					+ (std::pow(T, 1.0924131003119228) * (0.00002971599807934954 - 1.7857483385425343e-6 * std::log(T)))
//							/ std::exp(0.06287056567053795 * std::pow(std::log(T), 2)));
//	C13 = 5.64e-36 * std::pow(1 + z, 4) * (T - 2.73 * (1 + z));
//	dC13dT = 5.64e-36 * std::pow(1 + z, 4);
//
}

bool compute_next_state(std::array<double, NS> U0, std::array<double, NS> &U, double T, double z, double dt) {
	double K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19;
	double J20, J21, J22, J23, J24, J25, J26, J27, J28;

	const auto H0 = U0[nH];
	const auto Hp0 = U0[nHP];
	const auto He0 = U0[nHE];
	const auto Hep0 = U0[nHEP];
	const auto Hepp0 = U0[nHEPP];
	const auto ne0 = U0[nE];
	auto &ne = U[nE];
	auto &H = U[nH];
	auto &Hp = U[nHP];
	auto &He = U[nHE];
	auto &Hep = U[nHEP];
	auto &Hepp = U[nHEPP];
	U = U0;

	const auto Htot = H0 + Hp0;
	const auto Hetot = Hep0 + He0 + Hepp0;
	chemical_rates(K1, K2, K3, K4, K5, K6, T);
	radiation_rates(J20, J21, J22, z);
//	K1 = K2 = K3 = K4 = K5 = K6 = 0.0;
	double err;
#define List(a0,a1,a2,a3,a4,a5) {a0,a1,a2,a3,a4,a5}
#define Power std::pow
	std::array<double, NF> f = { -H + H0 - Hp + Hp0, Hp - Hp0 + dt * Hp * K2 * ne - dt * H * (J20 + K1 * ne), -He + He0 - Hep + Hep0 - Hepp + Hepp0, Hep - Hep0
			+ dt * Hep * (J22 + (K4 + K5) * ne) - dt * (He * J21 + He * K3 * ne + Hepp * K6 * ne), Hepp - Hepp0 + dt * Hepp * K6 * ne
			- dt * Hep * (J22 + K5 * ne), -Hep - 2 * Hepp - Hp + ne };
	do {

		const double a10 = -(dt * (J20 + K1 * ne));
		const double a11 = 1 + dt * K2 * ne;
		const double a15 = -(dt * (H * K1 - Hp * K2));
		const double a32 = -(dt * (J21 + K3 * ne));
		const double a33 = 1 - dt * (-J22 - K4 * ne - K5 * ne);
		const double a34 = -(dt * K6 * ne);
		const double a35 = -(dt * (He * K3 - Hep * K4 - Hep * K5 + Hepp * K6));
		const double a43 = -(dt * (J22 + K5 * ne));
		const double a44 = 1 + dt * K6 * ne;
		const double a45 = -(dt * (Hep * K5 - Hepp * K6));
		const double detA = a15 * (-(a34 * a43) + a32 * (a43 - a44) + a33 * a44)
				+ a11 * (-2 * a35 * a43 + a33 * a44 + a35 * a44 + a32 * (a43 - a44 - a45) + 2 * a33 * a45 - a34 * (a43 + a45))
				+ a10 * (2 * a35 * a43 - a33 * a44 - a35 * a44 - 2 * a33 * a45 + a34 * (a43 + a45) + a32 * (-a43 + a44 + a45));
		const std::array<std::array<double, NF>, NF> AinvdetA = { {
		List(
				a15 * (a34 * a43 - a33 * a44 + a32 * (-a43 + a44))
				+ a11 * (a34 * a43 + 2 * a35 * a43 - a35 * a44 + a34 * a45 + a32 * (-a43 + a44 + a45) - a33 * (a44 + 2 * a45)),
				2 * a35 * a43 - a33 * a44 - a35 * a44 - 2 * a33 * a45 + a34 * (a43 + a45) + a32 * (-a43 + a44 + a45), a15 * a32 * (-2 * a43 + a44),
				a15 * (-2 * a43 + a44), -(a15 * (a32 - 2 * a33 + a34)), a15 * (-(a34 * a43) + a32 * (a43 - a44) + a33 * a44)),
		List(a10 * (-2 * a35 * a43 + a33 * a44 + a35 * a44 + a32 * (a43 - a44 - a45) + 2 * a33 * a45 - a34 * (a43 + a45)),
				-2 * a35 * a43 + a33 * a44 + a35 * a44 + a32 * (a43 - a44 - a45) + 2 * a33 * a45 - a34 * (a43 + a45), a15 * a32 * (2 * a43 - a44),
				a15 * (2 * a43 - a44), a15 * (a32 - 2 * a33 + a34), a15 * (a34 * a43 - a33 * a44 + a32 * (-a43 + a44))),
		List(a10 * (a35 * (-a43 + a44) + (a33 - a34) * a45), a35 * (-a43 + a44) + (a33 - a34) * a45,
				a15 * (a34 * a43 - a33 * a44) + a10 * (-2 * a35 * a43 + a33 * a44 + a35 * a44 + 2 * a33 * a45 - a34 * (a43 + a45))
				+ a11 * (a34 * a43 + 2 * a35 * a43 - a35 * a44 + a34 * a45 - a33 * (a44 + 2 * a45)),
				(-a10 + a11 + a15) * (a43 - a44) + (a10 - a11) * a45, -((-a10 + a11 + a15) * (a33 - a34)) + (-a10 + a11) * a35,
				(a10 - a11) * (a35 * (a43 - a44) + (-a33 + a34) * a45)),
		List(-(a10 * (a35 * a44 + (a32 - a34) * a45)), -(a35 * a44) + (-a32 + a34) * a45,
				a32 * ((a11 + a15) * a44 + 2 * a11 * a45 - a10 * (a44 + 2 * a45)), (a11 + a15) * a44 + 2 * a11 * a45 - a10 * (a44 + 2 * a45),
				(-a10 + a11 + a15) * (a32 - a34) + 2 * (a10 - a11) * a35, (a10 - a11) * (a35 * a44 + (a32 - a34) * a45)),
		List(a10 * (a35 * a43 + (a32 - a33) * a45), a35 * a43 + (a32 - a33) * a45, a32 * (-(a15 * a43) + a10 * (a43 + a45) - a11 * (a43 + a45)),
				-(a15 * a43) + a10 * (a43 + a45) - a11 * (a43 + a45), -((-a10 + a11 + a15) * (a32 - a33)) + (-a10 + a11) * a35,
				-((a10 - a11) * (a35 * a43 + (a32 - a33) * a45))),
		List(a10 * (-(a34 * a43) + a32 * (a43 - a44) + a33 * a44), -(a34 * a43) + a32 * (a43 - a44) + a33 * a44, (a10 - a11) * a32 * (2 * a43 - a44),
				(a10 - a11) * (2 * a43 - a44), (a10 - a11) * (a32 - 2 * a33 + a34), (a10 - a11) * (a34 * a43 - a33 * a44 + a32 * (-a43 + a44))) } };
		for (int n = 0; n < NF; n++) {
			auto dU = 0.0;
			for (int m = 0; m < NF; m++) {
				dU += AinvdetA[n][m] * f[m] / detA;
			}
			U[n] = U[n] - dU;
			if (U[n] < 0.0) {
				return false;
			}
		}
		err = 0.0;
		f = { -H + H0 - Hp + Hp0, Hp - Hp0 + dt * Hp * K2 * ne - dt * H * (J20 + K1 * ne), -He + He0 - Hep + Hep0 - Hepp + Hepp0, Hep - Hep0
				+ dt * Hep * (J22 + (K4 + K5) * ne) - dt * (He * J21 + He * K3 * ne + Hepp * K6 * ne), Hepp - Hepp0 + dt * Hepp * K6 * ne
				- dt * Hep * (J22 + K5 * ne), -Hep - 2 * Hepp - Hp + ne };
		for (int i = 0; i < NF; i++) {
			err += f[i] * f[i];
		}
		err = std::sqrt(err);
		err /= (H + Hp + He + Hep + Hepp + ne);
	} while (err > 1.0e-9);
	return true;
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

double chemistry_update(const species sn, species &snp1, double e, double scale0, double scale1, double dt) {
	std::array<double, NS> Un;
	Un[nH] = sn.H;
	Un[nHP] = sn.Hp;
	Un[nHE] = sn.He;
	Un[nHEP] = sn.Hep;
	Un[nHEPP] = sn.Hepp;
	Un[nE] = sn.Hp + sn.Hep + 2 * sn.Hepp;

	std::array<double, NS> Unp1 = Un;
//	compute_next_state(Un, Unp1, species_T(Un, e), scale1, dt);
	double this_dt = dt;
	double tmax = dt;
	double t = 0.0;
	while (t < tmax) {
//		printf("stepping %e %e\n", t / tmax, this_dt / tmax);
		const auto w = (t + this_dt) / tmax;
		const auto scale = scale0 * (1.0 - w) + scale1 * w;
		if (compute_next_state(Un, Unp1, species_T(Un, e), scale, this_dt)) {
			t += this_dt;
			this_dt = std::min(tmax - t, this_dt * 2.0);
			Un = Unp1;
		} else {
			this_dt /= 2.0;
			Unp1 = Un;
		}
	}
	snp1.H = Unp1[nH];
	snp1.Hp = Unp1[nHP];
	snp1.He = Unp1[nHE];
	snp1.Hep = Unp1[nHEP];
	snp1.Hepp = Unp1[nHEPP];
	return species_T(snp1, e);
}

void chemistry_test() {
	species s0;
	const auto n = 1.0;
	s0.H = 0.92 * n;
	s0.Hp = 0.001 * n;
	s0.He = 0.08 * n;
	s0.Hep = 0.00 * n;
	s0.Hepp = 0.00 * n;
	double T0 = 1e+7;
	auto s = s0;
	printf("%14s %14s %14s %14s %14s %14s %14s %14s %14s %14s\n", "t", "T0", "T", "na", "ne", "H", "Hp", "He", "Hep", "Hepp");
	for (double dt = 1e1; dt < 1e18; dt *= 10.0) {
		auto T = chemistry_update(s0, s, species_energy(s0, T0), 0.0, 0.0, dt);
		const auto ne = s.Hp + s.Hep + 2 * s.Hepp;
		const auto na = s.Hp + s.H + 4 * s.He + 4 * s.Hep + 4 * s.Hepp;
		printf("%14.4e %14.4e %14.4e %14.4e %14.4e %14.4e %14.4e %14.4e %14.4e %14.4e\n", dt / (3600 * 24 * 365), T0, T, na, ne, s.H, s.Hp, s.He, s.Hep,
				s.Hepp);
	}

}
