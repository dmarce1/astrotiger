#include <cmath>
#include <stdio.h>
#include <vector>
#include <array>
#include <cassert>

#define KperEv (11604.525)

#define nH   0
#define nHP  1
#define nHN  2
#define nH2  3
#define nH2P  4
#define nHE  5
#define nHEP  6
#define nHEPP 7
#define NS  8

void rates(double &k1, double &k2, double &k3, double &k4, double &k5, double &k6, double &k7, double &k8,
		double &k9, double &k10, double &k11, double &k12, double &k13, double &k14, double &k15, double &k16,
		double &k17, double &k18, double &k19, double T, bool caseb) {
	const auto tev = T / KperEv;
	const auto logtev = std::log(T);
	const auto tiny = std::numeric_limits<double>::min();
	if (tev > 0.8) {
		k1 = exp(
				-32.71396786375 + 13.53655609057 * logtev - 5.739328757388 * std::pow(logtev, 2) + 1.563154982022 * std::pow(logtev, 3)
						- 0.2877056004391 * std::pow(logtev, 4) + 0.03482559773736999 * std::pow(logtev, 5) - 0.00263197617559 * std::pow(logtev, 6)
						+ 0.0001119543953861 * std::pow(logtev, 7) - 2.039149852002e-6 * std::pow(logtev, 8));

		k3 = exp(
				-44.09864886561001 + 23.91596563469 * logtev - 10.75323019821 * std::pow(logtev, 2) + 3.058038757198 * std::pow(logtev, 3)
						- 0.5685118909884001 * std::pow(logtev, 4) + 0.06795391233790001 * std::pow(logtev, 5) - 0.005009056101857001 * std::pow(logtev, 6)
						+ 0.0002067236157507 * std::pow(logtev, 7) - 3.649161410833d - 6 * std::pow(logtev, 8));

		k4 = (1.54e-9 * (1. + 0.3 / exp(8.099328789667 / tev)) / (exp(40.49664394833662 / tev) * std::pow(tev, 1.5)) + 3.92d - 13 / std::pow(tev, 0.6353));

		k5 = exp(
				-68.71040990212001 + 43.93347632635 * logtev - 18.48066993568 * std::pow(logtev, 2) + 4.701626486759002 * std::pow(logtev, 3)
						- 0.7692466334492 * std::pow(logtev, 4) + 0.08113042097303 * std::pow(logtev, 5) - 0.005324020628287001 * std::pow(logtev, 6)
						+ 0.0001975705312221 * std::pow(logtev, 7) - 3.165581065665e-6 * std::pow(logtev, 8));
	} else {
		k1 = tiny;
		k3 = tiny;
		k4 = 3.92e-13 / std::pow(tev, 0.6353);
		k5 = tiny;
	}

	if (caseb) {
		k4 = 1.26e-14 * std::pow((5.7067e5 / T), (0.75));
	}

	if (caseb) {
		if (T < 1.0e9) {
			k2 = 4.881357e-6 * std::pow(T, -1.5) * std::pow((1. + 1.14813e2 * std::pow(T, (-0.407))), (-2.242));
		} else {
			k2 = tiny;
		}
	} else {
		if (T > 5500.) {
			k2 = exp(
					-28.61303380689232 - 0.7241125657826851 * logtev - 0.02026044731984691 * std::pow(logtev, 2) - 0.002380861877349834 * std::pow(logtev, 3)
							- 0.0003212605213188796 * std::pow(logtev, 4) - 0.00001421502914054107 * std::pow(logtev, 5)
							+ 4.989108920299513e-6 * std::pow(logtev, 6) + 5.755614137575758d - 7 * std::pow(logtev, 7)
							- 1.856767039775261e-8 * std::pow(logtev, 8) - 3.071135243196595e-9 * std::pow(logtev, 9));
		} else {
			k2 = k4;
		}
	}
	if (caseb) {
		if (T < 1.0e9) {
			k6 = 7.8155e-5 * std::pow(T, -1.5) * std::pow((1. + 2.0189e2 * std::pow(T, -0.407)), (-2.242));
		} else {
			k6 = tiny;
		}
	} else {
		k6 = 3.36e-10 / sqrt(T) / std::pow((T / 1.e3), 0.2) / (1. + std::pow((T / 1.e6), 0.7));
	}
	k7 = 6.77e-15 * std::pow(tev, 0.8779);

	if (tev > 0.1) {
		k8 = exp(
				-20.06913897587003 + 0.2289800603272916 * logtev + 0.03599837721023835 * std::pow(logtev, 2) - 0.004555120027032095 * std::pow(logtev, 3)
						- 0.0003105115447124016 * std::pow(logtev, 4) + 0.0001073294010367247 * std::pow(logtev, 5) - 8.36671960467864e-6 * std::pow(logtev, 6)
						+ 2.238306228891639d - 7 * std::pow(logtev, 7));
	} else {
		k8 = 1.43e-9;
	}

	if (T > 6.7e3) {
		k9 = 5.81e-16 * std::pow(T / 56200., (-0.6657 * std::log10(T / 56200.)));
	} else {
		k9 = 1.85e-23 * std::pow(T, 1.8);
	}

	k10 = 6.0e-10;

	if (tev > 0.3) {
		k13 = 1.0670825e-10 * std::pow(tev, 2.012) / (exp(4.463 / tev) * std::pow((1. + 0.2472 * tev), 3.512));

		k11 = exp(
				-24.24914687731536 + 3.400824447095291 * logtev - 3.898003964650152 * std::pow(logtev, 2) + 2.045587822403071 * std::pow(logtev, 3)
						- 0.5416182856220388 * std::pow(logtev, 4) + 0.0841077503763412 * std::pow(logtev, 5) - 0.007879026154483455 * std::pow(logtev, 6)
						+ 0.0004138398421504563 * std::pow(logtev, 7) - 9.36345888928611e-6 * std::pow(logtev, 8));

		k12 = 4.38e-10 * exp(-102000. / T) * std::pow(T, 0.35);
	} else {
		k13 = tiny;
		k11 = tiny;
		k12 = tiny;
	}

	if (tev > 0.04) {
		k14 = exp(
				-18.01849334273 + 2.360852208681 * logtev - 0.2827443061704 * std::pow(logtev, 2) + 0.01623316639567 * std::pow(logtev, 3)
						- 0.03365012031362999 * std::pow(logtev, 4) + 0.01178329782711 * std::pow(logtev, 5) - 0.001656194699504 * std::pow(logtev, 6)
						+ 0.0001068275202678 * std::pow(logtev, 7) - 2.631285809207e-6 * std::pow(logtev, 8));
	} else {
		k14 = tiny;
	}

	if (tev > 0.1) {
		k15 = exp(
				-20.37260896533324 + 1.139449335841631 * logtev - 0.1421013521554148 * std::pow(logtev, 2) + 0.00846445538663 * std::pow(logtev, 3)
						- 0.0014327641212992 * std::pow(logtev, 4) + 0.0002012250284791 * std::pow(logtev, 5) + 0.0000866396324309 * std::pow(logtev, 6)
						- 0.00002585009680264 * std::pow(logtev, 7) + 2.4555011970392e-6 * std::pow(logtev, 8) - 8.06838246118e-8 * std::pow(logtev, 9));
	} else {
		k15 = 2.56e-9 * std::pow(tev, 1.78186);
	}

	k16 = 6.5e-9 / sqrt(tev);

	if (T > 1.0e4) {
		k17 = 4.0e-4 * std::pow(T, (-1.4)) * exp(-15100. / T);
	} else {
		k17 = 1.0e-8 * std::pow(T, (-0.4));

	}

	k18 = 5.56396e-8 / std::pow(tev, 0.6035);
	k19 = 4.64e-8 / sqrt(tev);

}

double compute_next_ne(std::array<double, NS> U0, std::array<double, NS> &U, double ne, double T, double dt) {
	double K1;
	double K2;
	double K3;
	double K4;
	double K5;
	double K6;
	double K7;
	double K8;
	double K9;
	double K10;
	double K11;
	double K12;
	double K13;
	double K14;
	double K15;
	double K16;
	double K17;
	double K18;
	double K19;
	rates(K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, T, false);
	auto &H0 = U0[nH];
	auto &Hp0 = U0[nHP];
	auto &Hn0 = U0[nHN];
	auto &H20 = U0[nH2];
	auto &H2p0 = U0[nH2P];
	auto &He0 = U0[nHE];
	auto &Hep0 = U0[nHEP];
	auto &Hepp0 = U0[nHEPP];
	auto &H = U[nH];
	auto &Hp = U[nHP];
	auto &Hn = U[nHN];
	auto &H2 = U[nH2];
	auto &H2p = U[nH2P];
	auto &He = U[nHE];
	auto &Hep = U[nHEP];
	auto &Hepp = U[nHEPP];
	const auto den = (1 + dt * K3 * ne + dt * K4 * ne + dt * K5 * ne + dt * K6 * ne + std::pow(dt, 2) * K3 * K5 * std::pow(ne, 2)
			+ std::pow(dt, 2) * K3 * K6 * std::pow(ne, 2) + std::pow(dt, 2) * K4 * K6 * std::pow(ne, 2));
	He = -((-(dt * K4 * ne * (dt * Hepp0 * K6 * ne + Hep0 * (1 + dt * K6 * ne)))
			+ He0 * (std::pow(dt, 2) * K5 * K6 * std::pow(ne, 2) - (1 + dt * K4 * ne + dt * K5 * ne) * (1 + dt * K6 * ne))) / den);
	Hep = -((-Hep0 - dt * He0 * K3 * ne - dt * Hep0 * K3 * ne - dt * Hep0 * K6 * ne - dt * Hepp0 * K6 * ne - std::pow(dt, 2) * He0 * K3 * K6 * std::pow(ne, 2)
			- std::pow(dt, 2) * Hep0 * K3 * K6 * std::pow(ne, 2) - std::pow(dt, 2) * Hepp0 * K3 * K6 * std::pow(ne, 2)) / den);
	Hepp = -((-Hepp0 - dt * Hepp0 * K3 * ne - dt * Hepp0 * K4 * ne - dt * Hep0 * K5 * ne - dt * Hepp0 * K5 * ne
			- std::pow(dt, 2) * He0 * K3 * K5 * std::pow(ne, 2) - std::pow(dt, 2) * Hep0 * K3 * K5 * std::pow(ne, 2)
			- std::pow(dt, 2) * Hepp0 * K3 * K5 * std::pow(ne, 2)) / den);
	double err = 0.0;
#define List(a,b,c,d,e) {a,b,c,d,e}
	do {

		const std::array<std::array<double, 5>, 5> A = { {
		List(1 + dt * (H2p * K10 - 2 * H2 * K13 - Hn * K15 + Hn * K8 + Hp * K9 + K1 * ne + K7 * ne),
				dt * (-(H2 * K11) - 2 * Hn * K16 + H * K9 - K2 * ne), dt * (-(H * K15) - 2 * Hp * K16 - H2p * K19 + H * K8 - K14 * ne),
				dt * (-(Hp * K11) - 2 * H * K13 - 2 * K12 * ne), dt * (H * K10 - Hn * K19 - 2 * K18 * ne)),
		List(dt * (-(H2p * K10) + Hp * K9 - K1 * ne), 1 + dt * (H2 * K11 + Hn * K16 + Hn * K17 + H * K9 + K2 * ne), dt * (Hp * K16 + Hp * K17),
				dt * Hp * K11, -(dt * H * K10)),
		List(dt*(Hn*K15 + Hn*K8 - K7*ne),dt*(Hn*K16 + Hn*K17),1 + dt*(H*K15 + Hp*K16 + Hp*K17 + H2p*K19 + H*K8 + K14*ne),0,dt*Hn*K19),
		List(dt*(-(H2p*K10) + H2*K13 - Hn*K8),dt*H2*K11,dt*(-(H2p*K19) - H*K8),1 + dt*(Hp*K11 + H*K13 + K12*ne),dt*(-(H*K10) - Hn*K19)),
		List(dt*(H2p*K10 - Hp*K9),dt*(-(H2*K11) - Hn*K17 - H*K9),dt*(-(Hp*K17) + H2p*K19),-(dt*Hp*K11),1 + dt*(H*K10 + Hn*K19 + K18*ne)) } };

		const std::array<double, 5> f = { H - H0
				+ dt
						* (H * H2p * K10 - H2 * Hp * K11 - 2 * H * H2 * K13 - H * Hn * K15 - 2 * Hn * Hp * K16 - H2p * Hn * K19 + H * Hn * K8 + H * Hp * K9
								+ H * K1 * ne - 2 * H2 * K12 * ne - Hn * K14 * ne - 2 * H2p * K18 * ne - Hp * K2 * ne + H * K7 * ne), Hp - Hp0
				+ dt * (-(H * H2p * K10) + H2 * Hp * K11 + Hn * Hp * K16 + Hn * Hp * K17 + H * Hp * K9 - H * K1 * ne + Hp * K2 * ne), Hn - Hn0
				+ dt * (H * Hn * K15 + Hn * Hp * K16 + Hn * Hp * K17 + H2p * Hn * K19 + H * Hn * K8 + Hn * K14 * ne - H * K7 * ne), H2 - H20
				+ dt * (-(H * H2p * K10) + H2 * Hp * K11 + H * H2 * K13 - H2p * Hn * K19 - H * Hn * K8 + H2 * K12 * ne), H2p - H2p0
				+ dt * (H * H2p * K10 - H2 * Hp * K11 - Hn * Hp * K17 + H2p * Hn * K19 - H * Hp * K9 + H2p * K18 * ne) };

		const auto a00 = A[0][0];
		const auto a01 = A[0][1];
		const auto a02 = A[0][2];
		const auto a03 = A[0][3];
		const auto a04 = A[0][4];
		const auto a10 = A[1][0];
		const auto a11 = A[1][1];
		const auto a12 = A[1][2];
		const auto a13 = A[1][3];
		const auto a14 = A[1][4];
		const auto a20 = A[2][0];
		const auto a21 = A[2][1];
		const auto a22 = A[2][2];
		const auto a23 = A[2][3];
		const auto a24 = A[2][4];
		const auto a30 = A[3][0];
		const auto a31 = A[3][1];
		const auto a32 = A[3][2];
		const auto a33 = A[3][3];
		const auto a34 = A[3][4];
		const auto a40 = A[4][0];
		const auto a41 = A[4][1];
		const auto a42 = A[4][2];
		const auto a43 = A[4][3];
		const auto a44 = A[4][4];

		const std::array<std::array<double, 5>, 5> invAdetA = { {
		List(
				a12 * a24 * a33 * a41 - a11 * a24 * a33 * a42 - a12 * a24 * a31 * a43 + a11 * a24 * a32 * a43 + a12 * a21 * a34 * a43
				- a11 * a22 * a34 * a43 + a14 * (-(a22 * a33 * a41) + a21 * a33 * a42 + a22 * a31 * a43 - a21 * a32 * a43)
				- a12 * a21 * a33 * a44 + a11 * a22 * a33 * a44
				+ a13 * (-(a24 * a32 * a41) + a22 * a34 * a41 + a24 * a31 * a42 - a21 * a34 * a42 - a22 * a31 * a44 + a21 * a32 * a44),
				-(a02 * a24 * a33 * a41) + a01 * a24 * a33 * a42 + a02 * a24 * a31 * a43 - a01 * a24 * a32 * a43 - a02 * a21 * a34 * a43
				+ a01 * a22 * a34 * a43 + a04 * (a22 * a33 * a41 - a21 * a33 * a42 - a22 * a31 * a43 + a21 * a32 * a43) + a02 * a21 * a33 * a44
				- a01 * a22 * a33 * a44
				+ a03 * (a24 * a32 * a41 - a22 * a34 * a41 - a24 * a31 * a42 + a21 * a34 * a42 + a22 * a31 * a44 - a21 * a32 * a44),
				a02 * a14 * a33 * a41 - a02 * a13 * a34 * a41 - a01 * a14 * a33 * a42 + a01 * a13 * a34 * a42 - a02 * a14 * a31 * a43
				+ a01 * a14 * a32 * a43 + a02 * a11 * a34 * a43 - a01 * a12 * a34 * a43
				+ a04 * (a13 * a32 * a41 - a12 * a33 * a41 - a13 * a31 * a42 + a11 * a33 * a42 + a12 * a31 * a43 - a11 * a32 * a43)
				+ (a02 * a13 * a31 - a01 * a13 * a32 - a02 * a11 * a33 + a01 * a12 * a33) * a44
				+ a03 * (-(a14 * a32 * a41) + a12 * a34 * a41 + a14 * a31 * a42 - a11 * a34 * a42 - a12 * a31 * a44 + a11 * a32 * a44),
				a02 * a13 * a24 * a41 - a01 * a13 * a24 * a42 + a02 * a14 * a21 * a43 - a01 * a14 * a22 * a43 - a02 * a11 * a24 * a43
				+ a01 * a12 * a24 * a43 + a04 * (-(a13 * a22 * a41) + a13 * a21 * a42 - a12 * a21 * a43 + a11 * a22 * a43)
				- a02 * a13 * a21 * a44 + a01 * a13 * a22 * a44
				+ a03 * (a14 * a22 * a41 - a12 * a24 * a41 - a14 * a21 * a42 + a11 * a24 * a42 + a12 * a21 * a44 - a11 * a22 * a44),
				-(a02 * a13 * a24 * a31) + a01 * a13 * a24 * a32 - a02 * a14 * a21 * a33 + a01 * a14 * a22 * a33 + a02 * a11 * a24 * a33
				- a01 * a12 * a24 * a33 + a04 * (a13 * a22 * a31 - a13 * a21 * a32 + a12 * a21 * a33 - a11 * a22 * a33) + a02 * a13 * a21 * a34
				- a01 * a13 * a22 * a34
				+ a03 * (-(a14 * a22 * a31) + a12 * a24 * a31 + a14 * a21 * a32 - a11 * a24 * a32 - a12 * a21 * a34 + a11 * a22 * a34)),
		List(
				-(a12 * a24 * a33 * a40) + a10 * a24 * a33 * a42 + a12 * a24 * a30 * a43 - a10 * a24 * a32 * a43 - a12 * a20 * a34 * a43
				+ a10 * a22 * a34 * a43 + a14 * (a22 * a33 * a40 - a20 * a33 * a42 - a22 * a30 * a43 + a20 * a32 * a43) + a12 * a20 * a33 * a44
				- a10 * a22 * a33 * a44
				+ a13 * (a24 * a32 * a40 - a22 * a34 * a40 - a24 * a30 * a42 + a20 * a34 * a42 + a22 * a30 * a44 - a20 * a32 * a44),
				a02 * a24 * a33 * a40 - a00 * a24 * a33 * a42 - a02 * a24 * a30 * a43 + a00 * a24 * a32 * a43 + a02 * a20 * a34 * a43
				- a00 * a22 * a34 * a43 + a04 * (-(a22 * a33 * a40) + a20 * a33 * a42 + a22 * a30 * a43 - a20 * a32 * a43)
				- a02 * a20 * a33 * a44 + a00 * a22 * a33 * a44
				+ a03 * (-(a24 * a32 * a40) + a22 * a34 * a40 + a24 * a30 * a42 - a20 * a34 * a42 - a22 * a30 * a44 + a20 * a32 * a44),
				-(a02 * a14 * a33 * a40) + a02 * a13 * a34 * a40 + a00 * a14 * a33 * a42 - a00 * a13 * a34 * a42 + a02 * a14 * a30 * a43
				- a00 * a14 * a32 * a43 - a02 * a10 * a34 * a43 + a00 * a12 * a34 * a43
				+ a04 * (-(a13 * a32 * a40) + a12 * a33 * a40 + a13 * a30 * a42 - a10 * a33 * a42 - a12 * a30 * a43 + a10 * a32 * a43)
				- a02 * a13 * a30 * a44 + a00 * a13 * a32 * a44 + a02 * a10 * a33 * a44 - a00 * a12 * a33 * a44
				+ a03 * (a14 * a32 * a40 - a12 * a34 * a40 - a14 * a30 * a42 + a10 * a34 * a42 + a12 * a30 * a44 - a10 * a32 * a44),
				-(a02 * a13 * a24 * a40) + a00 * a13 * a24 * a42 - a02 * a14 * a20 * a43 + a00 * a14 * a22 * a43 + a02 * a10 * a24 * a43
				- a00 * a12 * a24 * a43 + a04 * (a13 * a22 * a40 - a13 * a20 * a42 + a12 * a20 * a43 - a10 * a22 * a43) + a02 * a13 * a20 * a44
				- a00 * a13 * a22 * a44
				+ a03 * (-(a14 * a22 * a40) + a12 * a24 * a40 + a14 * a20 * a42 - a10 * a24 * a42 - a12 * a20 * a44 + a10 * a22 * a44),
				a02 * a13 * a24 * a30 - a00 * a13 * a24 * a32 + a02 * a14 * a20 * a33 - a00 * a14 * a22 * a33 - a02 * a10 * a24 * a33
				+ a00 * a12 * a24 * a33 + a04 * (-(a13 * a22 * a30) + a13 * a20 * a32 - a12 * a20 * a33 + a10 * a22 * a33)
				- a02 * a13 * a20 * a34 + a00 * a13 * a22 * a34
				+ a03 * (a14 * a22 * a30 - a12 * a24 * a30 - a14 * a20 * a32 + a10 * a24 * a32 + a12 * a20 * a34 - a10 * a22 * a34)),
		List(
				a11 * a24 * a33 * a40 - a10 * a24 * a33 * a41 - a11 * a24 * a30 * a43 + a10 * a24 * a31 * a43 + a11 * a20 * a34 * a43
				- a10 * a21 * a34 * a43 + a14 * (-(a21 * a33 * a40) + a20 * a33 * a41 + a21 * a30 * a43 - a20 * a31 * a43)
				- a11 * a20 * a33 * a44 + a10 * a21 * a33 * a44
				+ a13 * (-(a24 * a31 * a40) + a21 * a34 * a40 + a24 * a30 * a41 - a20 * a34 * a41 - a21 * a30 * a44 + a20 * a31 * a44),
				-(a01 * a24 * a33 * a40) + a00 * a24 * a33 * a41 + a01 * a24 * a30 * a43 - a00 * a24 * a31 * a43 - a01 * a20 * a34 * a43
				+ a00 * a21 * a34 * a43 + a04 * (a21 * a33 * a40 - a20 * a33 * a41 - a21 * a30 * a43 + a20 * a31 * a43) + a01 * a20 * a33 * a44
				- a00 * a21 * a33 * a44
				+ a03 * (a24 * a31 * a40 - a21 * a34 * a40 - a24 * a30 * a41 + a20 * a34 * a41 + a21 * a30 * a44 - a20 * a31 * a44),
				a01 * a14 * a33 * a40 - a01 * a13 * a34 * a40 - a00 * a14 * a33 * a41 + a00 * a13 * a34 * a41 - a01 * a14 * a30 * a43
				+ a00 * a14 * a31 * a43 + a01 * a10 * a34 * a43 - a00 * a11 * a34 * a43
				+ a04 * (a13 * a31 * a40 - a11 * a33 * a40 - a13 * a30 * a41 + a10 * a33 * a41 + a11 * a30 * a43 - a10 * a31 * a43)
				+ (a01 * a13 * a30 - a00 * a13 * a31 - a01 * a10 * a33 + a00 * a11 * a33) * a44
				+ a03 * (-(a14 * a31 * a40) + a11 * a34 * a40 + a14 * a30 * a41 - a10 * a34 * a41 - a11 * a30 * a44 + a10 * a31 * a44),
				a01 * a13 * a24 * a40 - a00 * a13 * a24 * a41 + a01 * a14 * a20 * a43 - a00 * a14 * a21 * a43 - a01 * a10 * a24 * a43
				+ a00 * a11 * a24 * a43 + a04 * (-(a13 * a21 * a40) + a13 * a20 * a41 - a11 * a20 * a43 + a10 * a21 * a43)
				- a01 * a13 * a20 * a44 + a00 * a13 * a21 * a44
				+ a03 * (a14 * a21 * a40 - a11 * a24 * a40 - a14 * a20 * a41 + a10 * a24 * a41 + a11 * a20 * a44 - a10 * a21 * a44),
				-(a01 * a13 * a24 * a30) + a00 * a13 * a24 * a31 - a01 * a14 * a20 * a33 + a00 * a14 * a21 * a33 + a01 * a10 * a24 * a33
				- a00 * a11 * a24 * a33 + a04 * (a13 * a21 * a30 - a13 * a20 * a31 + a11 * a20 * a33 - a10 * a21 * a33) + a01 * a13 * a20 * a34
				- a00 * a13 * a21 * a34
				+ a03 * (-(a14 * a21 * a30) + a11 * a24 * a30 + a14 * a20 * a31 - a10 * a24 * a31 - a11 * a20 * a34 + a10 * a21 * a34)),
		List(
				-(a11 * a24 * a32 * a40) + a11 * a22 * a34 * a40 + a10 * a24 * a32 * a41 - a10 * a22 * a34 * a41 + a11 * a24 * a30 * a42
				- a10 * a24 * a31 * a42 - a11 * a20 * a34 * a42 + a10 * a21 * a34 * a42
				+ a14 * (-(a22 * a31 * a40) + a21 * a32 * a40 + a22 * a30 * a41 - a20 * a32 * a41 - a21 * a30 * a42 + a20 * a31 * a42)
				- a11 * a22 * a30 * a44 + a10 * a22 * a31 * a44 + a11 * a20 * a32 * a44 - a10 * a21 * a32 * a44
				+ a12 * (a24 * a31 * a40 - a21 * a34 * a40 - a24 * a30 * a41 + a20 * a34 * a41 + a21 * a30 * a44 - a20 * a31 * a44),
				a01 * a24 * a32 * a40 - a01 * a22 * a34 * a40 - a00 * a24 * a32 * a41 + a00 * a22 * a34 * a41 - a01 * a24 * a30 * a42
				+ a00 * a24 * a31 * a42 + a01 * a20 * a34 * a42 - a00 * a21 * a34 * a42
				+ a04 * (a22 * a31 * a40 - a21 * a32 * a40 - a22 * a30 * a41 + a20 * a32 * a41 + a21 * a30 * a42 - a20 * a31 * a42)
				+ (a01 * a22 * a30 - a00 * a22 * a31 - a01 * a20 * a32 + a00 * a21 * a32) * a44
				+ a02 * (-(a24 * a31 * a40) + a21 * a34 * a40 + a24 * a30 * a41 - a20 * a34 * a41 - a21 * a30 * a44 + a20 * a31 * a44),
				-(a01 * a14 * a32 * a40) + a01 * a12 * a34 * a40 + a00 * a14 * a32 * a41 - a00 * a12 * a34 * a41 + a01 * a14 * a30 * a42
				- a00 * a14 * a31 * a42 - a01 * a10 * a34 * a42 + a00 * a11 * a34 * a42
				+ a04 * (-(a12 * a31 * a40) + a11 * a32 * a40 + a12 * a30 * a41 - a10 * a32 * a41 - a11 * a30 * a42 + a10 * a31 * a42)
				- a01 * a12 * a30 * a44 + a00 * a12 * a31 * a44 + a01 * a10 * a32 * a44 - a00 * a11 * a32 * a44
				+ a02 * (a14 * a31 * a40 - a11 * a34 * a40 - a14 * a30 * a41 + a10 * a34 * a41 + a11 * a30 * a44 - a10 * a31 * a44),
				a01 * a14 * a22 * a40 - a01 * a12 * a24 * a40 - a00 * a14 * a22 * a41 + a00 * a12 * a24 * a41 - a01 * a14 * a20 * a42
				+ a00 * a14 * a21 * a42 + a01 * a10 * a24 * a42 - a00 * a11 * a24 * a42
				+ a04 * (a12 * a21 * a40 - a11 * a22 * a40 - a12 * a20 * a41 + a10 * a22 * a41 + a11 * a20 * a42 - a10 * a21 * a42)
				+ (a01 * a12 * a20 - a00 * a12 * a21 - a01 * a10 * a22 + a00 * a11 * a22) * a44
				+ a02 * (-(a14 * a21 * a40) + a11 * a24 * a40 + a14 * a20 * a41 - a10 * a24 * a41 - a11 * a20 * a44 + a10 * a21 * a44),
				-(a01 * a14 * a22 * a30) + a01 * a12 * a24 * a30 + a00 * a14 * a22 * a31 - a00 * a12 * a24 * a31 + a01 * a14 * a20 * a32
				- a00 * a14 * a21 * a32 - a01 * a10 * a24 * a32 + a00 * a11 * a24 * a32
				+ a04 * (-(a12 * a21 * a30) + a11 * a22 * a30 + a12 * a20 * a31 - a10 * a22 * a31 - a11 * a20 * a32 + a10 * a21 * a32)
				- a01 * a12 * a20 * a34 + a00 * a12 * a21 * a34 + a01 * a10 * a22 * a34 - a00 * a11 * a22 * a34
				+ a02 * (a14 * a21 * a30 - a11 * a24 * a30 - a14 * a20 * a31 + a10 * a24 * a31 + a11 * a20 * a34 - a10 * a21 * a34)),
		List(
				a33 * (a12 * a21 * a40 - a11 * a22 * a40 - a12 * a20 * a41 + a10 * a22 * a41 + a11 * a20 * a42 - a10 * a21 * a42)
				+ a13 * (a22 * a31 * a40 - a21 * a32 * a40 - a22 * a30 * a41 + a20 * a32 * a41 + a21 * a30 * a42 - a20 * a31 * a42)
				+ (-(a12 * a21 * a30) + a11 * a22 * a30 + a12 * a20 * a31 - a10 * a22 * a31 - a11 * a20 * a32 + a10 * a21 * a32) * a43,
				a33 * (-(a02 * a21 * a40) + a01 * a22 * a40 + a02 * a20 * a41 - a00 * a22 * a41 - a01 * a20 * a42 + a00 * a21 * a42)
				+ a03 * (-(a22 * a31 * a40) + a21 * a32 * a40 + a22 * a30 * a41 - a20 * a32 * a41 - a21 * a30 * a42 + a20 * a31 * a42)
				+ (a02 * a21 * a30 - a01 * a22 * a30 - a02 * a20 * a31 + a00 * a22 * a31 + a01 * a20 * a32 - a00 * a21 * a32) * a43,
				a01 * a13 * a32 * a40 - a01 * a12 * a33 * a40 - a00 * a13 * a32 * a41 + a00 * a12 * a33 * a41 - a01 * a13 * a30 * a42
				+ a00 * a13 * a31 * a42 + a01 * a10 * a33 * a42 - a00 * a11 * a33 * a42
				+ a03 * (a12 * a31 * a40 - a11 * a32 * a40 - a12 * a30 * a41 + a10 * a32 * a41 + a11 * a30 * a42 - a10 * a31 * a42)
				+ (a01 * a12 * a30 - a00 * a12 * a31 - a01 * a10 * a32 + a00 * a11 * a32) * a43
				+ a02 * (-(a13 * a31 * a40) + a11 * a33 * a40 + a13 * a30 * a41 - a10 * a33 * a41 - a11 * a30 * a43 + a10 * a31 * a43),
				a13 * (a02 * a21 * a40 - a01 * a22 * a40 - a02 * a20 * a41 + a00 * a22 * a41 + a01 * a20 * a42 - a00 * a21 * a42)
				+ a03 * (-(a12 * a21 * a40) + a11 * a22 * a40 + a12 * a20 * a41 - a10 * a22 * a41 - a11 * a20 * a42 + a10 * a21 * a42)
				+ (a02 * a11 * a20 - a01 * a12 * a20 - a02 * a10 * a21 + a00 * a12 * a21 + a01 * a10 * a22 - a00 * a11 * a22) * a43,
				a13 * (-(a02 * a21 * a30) + a01 * a22 * a30 + a02 * a20 * a31 - a00 * a22 * a31 - a01 * a20 * a32 + a00 * a21 * a32)
				+ a03 * (a12 * a21 * a30 - a11 * a22 * a30 - a12 * a20 * a31 + a10 * a22 * a31 + a11 * a20 * a32 - a10 * a21 * a32)
				+ (-(a02 * a11 * a20) + a01 * a12 * a20 + a02 * a10 * a21 - a00 * a12 * a21 - a01 * a10 * a22 + a00 * a11 * a22) * a33) } };

		const auto detA = -(a02 * a13 * a24 * a31 * a40) + a01 * a13 * a24 * a32 * a40 - a02 * a14 * a21 * a33 * a40 + a01 * a14 * a22 * a33 * a40
				+ a02 * a11 * a24 * a33 * a40 - a01 * a12 * a24 * a33 * a40 + a02 * a13 * a21 * a34 * a40 - a01 * a13 * a22 * a34 * a40
				+ a02 * a13 * a24 * a30 * a41 - a00 * a13 * a24 * a32 * a41 + a02 * a14 * a20 * a33 * a41 - a00 * a14 * a22 * a33 * a41
				- a02 * a10 * a24 * a33 * a41 + a00 * a12 * a24 * a33 * a41 - a02 * a13 * a20 * a34 * a41 + a00 * a13 * a22 * a34 * a41
				- a01 * a13 * a24 * a30 * a42 + a00 * a13 * a24 * a31 * a42 - a01 * a14 * a20 * a33 * a42 + a00 * a14 * a21 * a33 * a42
				+ a01 * a10 * a24 * a33 * a42 - a00 * a11 * a24 * a33 * a42 + a01 * a13 * a20 * a34 * a42 - a00 * a13 * a21 * a34 * a42
				+ a02 * a14 * a21 * a30 * a43 - a01 * a14 * a22 * a30 * a43 - a02 * a11 * a24 * a30 * a43 + a01 * a12 * a24 * a30 * a43
				- a02 * a14 * a20 * a31 * a43 + a00 * a14 * a22 * a31 * a43 + a02 * a10 * a24 * a31 * a43 - a00 * a12 * a24 * a31 * a43
				+ a01 * a14 * a20 * a32 * a43 - a00 * a14 * a21 * a32 * a43 - a01 * a10 * a24 * a32 * a43 + a00 * a11 * a24 * a32 * a43
				+ a02 * a11 * a20 * a34 * a43 - a01 * a12 * a20 * a34 * a43 - a02 * a10 * a21 * a34 * a43 + a00 * a12 * a21 * a34 * a43
				+ a01 * a10 * a22 * a34 * a43 - a00 * a11 * a22 * a34 * a43

				+ a04
						* (a33 * (a12 * a21 * a40 - a11 * a22 * a40 - a12 * a20 * a41 + a10 * a22 * a41 + a11 * a20 * a42 - a10 * a21 * a42)
								+ a13 * (a22 * a31 * a40 - a21 * a32 * a40 - a22 * a30 * a41 + a20 * a32 * a41 + a21 * a30 * a42 - a20 * a31 * a42)
								+ (-(a12 * a21 * a30) + a11 * a22 * a30 + a12 * a20 * a31 - a10 * a22 * a31 - a11 * a20 * a32 + a10 * a21 * a32) * a43)
				+ a13 * (-(a02 * a21 * a30) + a01 * a22 * a30 + a02 * a20 * a31 - a00 * a22 * a31 - a01 * a20 * a32 + a00 * a21 * a32) * a44
				+ (-(a02 * a11 * a20) + a01 * a12 * a20 + a02 * a10 * a21 - a00 * a12 * a21 - a01 * a10 * a22 + a00 * a11 * a22) * a33 * a44
				+ a03
						* (a12 * a24 * a31 * a40 - a11 * a24 * a32 * a40 - a12 * a21 * a34 * a40 + a11 * a22 * a34 * a40 - a12 * a24 * a30 * a41
								+ a10 * a24 * a32 * a41 + a12 * a20 * a34 * a41 - a10 * a22 * a34 * a41 + a11 * a24 * a30 * a42 - a10 * a24 * a31 * a42
								- a11 * a20 * a34 * a42 + a10 * a21 * a34 * a42
								+ a14 * (-(a22 * a31 * a40) + a21 * a32 * a40 + a22 * a30 * a41 - a20 * a32 * a41 - a21 * a30 * a42 + a20 * a31 * a42)
								+ (a12 * a21 * a30 - a11 * a22 * a30 - a12 * a20 * a31 + a10 * a22 * a31 + a11 * a20 * a32 - a10 * a21 * a32) * a44);

		auto U1 = U;
		for (int n = 0; n < 5; n++) {
			for (int m = 0; m < 5; m++) {
				U1[n] -= invAdetA[n][m] * f[m] / detA;
				U1[n] = std::max(U1[n], 0.0);
			}
		}
		U = U1;
		auto &H1 = U1[nH];
		auto &Hp1 = U1[nHP];
		auto &H21 = U1[nH2];
		auto &H2p1 = U1[nH2P];
		auto &Hn1 = U1[nHN];
		auto Ht1 = std::max(H1, Hp1);
		Ht1 = std::max(Ht1, Hn1);
		Ht1 = std::max(Ht1, H21);
		Ht1 = std::max(Ht1, H2p1);
		auto Ht2 = std::max(H1, Hp);
		Ht2 = std::max(Ht1, Hn);
		Ht2 = std::max(Ht1, H2);
		Ht2 = std::max(Ht1, H2p);
		H = H1;
		Hp = Hp1;
		H2 = H21;
		err = std::abs(std::log(std::abs(Ht1 / Ht2)));
	} while (err > 1.0e-10);
	Hn = K7 * H * ne / (K8 * H + K16 * Hp + K14 * ne);
	return ne;
}

double compute(const std::array<double, NS> U0, std::array<double, NS> &U, double T, double dt) {
	auto &H = U[nH];
	auto &Hp = U[nHP];
	auto &Hn = U[nHN];
	auto &H2 = U[nH2];
	auto &H2p = U[nH2P];
	auto &He = U[nHE];
	auto &Hep = U[nHEP];
	auto &Hepp = U[nHEPP];

	auto nemax = H + Hp + 2.0 * Hn + 2.0 * H2 + 2.0 * H2p + 2.0 * He + 2.0 * Hep + 2.0 * Hepp;
	auto nemin = 0.0;
	double nemid;
	auto U1 = U;
	do {
		nemid = 0.5 * (nemax + nemin);
		auto ne = Hp - Hn + Hep + 2.0 * Hepp + H2p;
		U1 = U;
		const auto fmax = compute_next_ne(U0, U1, nemax, T, dt) - ne;
		U1 = U;
		const auto fmid = compute_next_ne(U0, U1, nemid, T, dt) - ne;
		if (fmax * fmid > 0.0) {
			nemax = nemid;
		} else {
			nemin = nemid;
		}
	} while (nemin / nemax < 0.99999);
	U = U1;
//	printf( "%e %e\n", U[nH], U[nHP]);
	return nemid;
}

void chemistry_test() {
	std::array<double, NS> U;
	const auto dt = 1.0;
	printf("%14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s\n", "time", "T", "A", "AH", "AHe", "N", "Ne", "H", "H+", "H-", "H2",
			"H2p", "He", "He+", "He++");
	double T = 1.0e4;
	for (double dt = 1.0e+3; dt <= 1e+17; dt *= 10) {
		for (int i = 0; i < NS; i++) {
			U[i] = 0.0;
		}
		double ne = 1e-1;
		U[nH] = 0.01 * ne;
		U[nHP] = 0.91 * ne;
		U[nH2P] = 0.00 * ne;
		U[nHE] = 0.00 * ne;
		U[nHEP] = 0.08 * ne;
		U[nHEPP] = 0.0 * ne;
		auto U0 = U;
		const auto eint = 1.0e-8 * ne;
		//	T /= 5.0;
		compute(U0, U, T, dt);
		ne = U[nHP] + U[nHEP] + 2.0 * U[nHEPP] + U[nH2P] - U[nHN];
		const auto nnuc = U[nH] + U[nHP] + U[nHN] + 2.0 * U[nH2] + 2.0 * U[nH2P] + 4.0 * U[nHE] + 4.0 * U[nHEP] + 4.0 * U[nHEPP];
		const auto nnucleus = U[nH] + U[nHP] + U[nHN] + U[nH2] + U[nH2P] + U[nHE] + U[nHEP] + U[nHEPP];
		const auto n = U[nH] + U[nHP] + U[nHN] + U[nH2] + U[nH2P] + U[nHE] + U[nHEP] + U[nHEPP] + ne;
		const auto nHnuc = U[nH] + U[nHP] + U[nHN] + 2.0 * U[nH2] + 2.0 * U[nH2P];
		const auto nHenuc = 4.0 * U[nHE] + 4.0 * U[nHEP] + 4.0 * U[nHEPP];
		printf("%14.5e %14.5e %14.5e %14.5e %14.5e  %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e\n", (double) dt, (double) T,
				(double) nnuc, (double) nHnuc, (double) nHenuc, (double) n, (double) ne, (double) U[nH], (double) U[nHP], (double) U[nHN], (double) U[nH2],
				(double) U[nH2P], (double) U[nHE], (double) U[nHEP], (double) U[nHEPP]);
//		U0 = U;
	}
}
