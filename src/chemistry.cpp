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

double k1(double T) {
	T /= KperEv;
	const auto c0 = -32.71396786;
	const auto c1 = 13.536556;
	const auto c2 = -5.73932875;
	const auto c3 = 1.56315498;
	const auto c4 = -0.2877056;
	const auto c5 = 3.48255977e-2;
	const auto c6 = -2.63197617e-3;
	const auto c7 = 1.11954395e-4;
	const auto c8 = -2.03914985e-6;
	const auto lnT = std::log(T);
	double k = c8;
	k = k * lnT + c7;
	k = k * lnT + c6;
	k = k * lnT + c5;
	k = k * lnT + c4;
	k = k * lnT + c3;
	k = k * lnT + c2;
	k = k * lnT + c1;
	k = k * lnT + c0;
	return std::exp(k);
}

double k2(double T) {
	T /= KperEv;
	const auto c0 = -28.6130338;
	const auto c1 = -0.72411256;
	const auto c2 = -2.02604473e-2;
	const auto c3 = -2.38086188e-3;
	const auto c4 = -3.21260521e-4;
	const auto c5 = -1.42150291e-5;
	const auto c6 = 4.98910892e-6;
	const auto c7 = 5.75561414e-7;
	const auto c8 = -1.85676704e-8;
	const auto c9 = -3.07113524e-9;
	const auto lnT = std::log(T);
	double k = c9;
	k = k * lnT + c8;
	k = k * lnT + c7;
	k = k * lnT + c6;
	k = k * lnT + c5;
	k = k * lnT + c4;
	k = k * lnT + c3;
	k = k * lnT + c2;
	k = k * lnT + c1;
	k = k * lnT + c0;
	return std::exp(k);
}

double k3(double T) {
	T /= KperEv;
	const auto c0 = -44.09864886;
	const auto c1 = 23.91596563;
	const auto c2 = -10.7532302;
	const auto c3 = 3.05803875;
	const auto c4 = -0.56851189;
	const auto c5 = 6.79539123e-2;
	const auto c6 = -5.00905610e-3;
	const auto c7 = 2.06723616e-4;
	const auto c8 = -3.64916141e-6;
	const auto lnT = std::log(T);
	double k = c8;
	k = k * lnT + c7;
	k = k * lnT + c6;
	k = k * lnT + c5;
	k = k * lnT + c4;
	k = k * lnT + c3;
	k = k * lnT + c2;
	k = k * lnT + c1;
	k = k * lnT + c0;
	return std::exp(k);

}

double k4r(double T) {

	T /= KperEv;
	return 3.925E-13 * std::pow(T, -0.6353);
}

double k4d(double T) {

	T /= KperEv;
	return 1.544E-9 * std::pow(T, -1.5) * std::exp(-48.596 / T) * (0.3 + std::exp(8.10 / T));
}

double k4(double T) {

	return k4r(T) + k4d(T);
}

double k5(double T) {

	T /= KperEv;
	const auto c0 = -68.71040990;
	const auto c1 = 43.93347633;
	const auto c2 = -18.4806699;
	const auto c3 = 4.70162649;
	const auto c4 = -0.76924663;
	const auto c5 = 8.113042E-2;
	const auto c6 = -5.32402063E-3;
	const auto c7 = 1.97570531E-4;
	const auto c8 = -3.16558106E-6;
	const auto lnT = std::log(T);
	double k = c8;
	k = k * lnT + c7;
	k = k * lnT + c6;
	k = k * lnT + c5;
	k = k * lnT + c4;
	k = k * lnT + c3;
	k = k * lnT + c2;
	k = k * lnT + c1;
	k = k * lnT + c0;
	return std::exp(k);
}

double k6(double T) {

	return 3.36e-10 * std::pow(T, -0.5) * std::pow(T / 1000, -0.2) / (1.0 + std::pow(T / 1e6, 0.7));
}

double k7(double T) {
// VERIFIED
	if (T < 6000) {
		const auto logT = std::log10(T);
		const auto c0 = 1.429e-18;
		const auto c1 = std::pow(T, 0.762);
		const auto c2 = std::pow(T, 0.1523 * logT);
		const auto c3 = std::pow(T, -3.247e-2 * logT * logT);
		return c0 * c1 * c2 * c3;
	} else {
		const auto logT = std::log10(T);
		const auto c0 = 3.802e-17;
		const auto c1 = std::pow(T, 0.1998 * logT);
		const auto c2 = std::pow(10, 4.0415e-5 * std::pow(logT, 6));
		const auto c3 = std::pow(10, -5.447e-3 * std::pow(logT, 4));
		return c0 * c1 * c2 * c3;
	}
}

double k8(double T) {

	T /= KperEv;
	if (T > 0.1) {
		const auto c0 = -20.06913897;
		const auto c1 = 0.22898;
		const auto c2 = 3.5998377E-2;
		const auto c3 = -4.55512E-3;
		const auto c4 = -3.10511544E-4;
		const auto c5 = 1.0732940E-4;
		const auto c6 = -8.36671960E-6;
		const auto c7 = 2.23830623E-7;
		const auto lnT = std::log(T);
		double k = c7;
		k = k * lnT + c6;
		k = k * lnT + c5;
		k = k * lnT + c4;
		k = k * lnT + c3;
		k = k * lnT + c2;
		k = k * lnT + c1;
		k = k * lnT + c0;
		return std::exp(k);
	} else {
		return 1.428e-9;
	}
}

double k11(double T) {

	T /= KperEv;
	const auto c0 = -24.24914687;
	const auto c1 = 3.40082444;
	const auto c2 = -3.89800396;
	const auto c3 = 2.04558782;
	const auto c4 = -0.541618285;
	const auto c5 = 8.41077503E-2;
	const auto c6 = -7.87902615E-3;
	const auto c7 = 4.13839842E-4;
	const auto c8 = -9.36345888E-6;
	const auto lnT = std::log(T);
	double k = c8;
	k = k * lnT + c7;
	k = k * lnT + c6;
	k = k * lnT + c5;
	k = k * lnT + c4;
	k = k * lnT + c3;
	k = k * lnT + c2;
	k = k * lnT + c1;
	k = k * lnT + c0;
	return std::exp(k);

}

double k12(double T) {

	auto rc = 5.6e-11 * std::sqrt(T) * std::exp(-102124 / T);
	return rc;
}

double k14(double T) {

	T /= KperEv;
	const auto c0 = -18.01849334;
	const auto c1 = 2.3608522;
	const auto c2 = -0.28274430;
	const auto c3 = 1.62331664E-2;
	const auto c4 = -3.36501203E-2;
	const auto c5 = 1.17832978E-2;
	const auto c6 = -1.65619470E-3;
	const auto c7 = 1.06827520E-4;
	const auto c8 = -2.63128581E-6;
	const auto lnT = std::log(T);
	double k = c8;
	k = k * lnT + c7;
	k = k * lnT + c6;
	k = k * lnT + c5;
	k = k * lnT + c4;
	k = k * lnT + c3;
	k = k * lnT + c2;
	k = k * lnT + c1;
	k = k * lnT + c0;
	return std::exp(k);
}

double k16(double T) {
	return 7e-8 * std::pow(T / 100, -0.5);
}

double k9(double T) {
	;

	T = 6700;
	if (T < 6700) {
		return 1.85e-23 * std::pow(T, 1.8);
	} else {
		return 5.81e-16 * std::pow(T / 56200.0, -0.6657 * std::log10(T / 56200.0));
	}
}

double k10(double T) {

	;
	return 6.4E-10;
}

double k13(double T) {
	;

	T /= KperEv;
	return 1.067e-10 * std::pow(T, 2.012) * std::exp(-(4.463 / T) * std::pow(1 + 0.2472 * T, 3.512));
}

double k15(double T) {
	;

	T /= KperEv;
	if (T > 0.1) {
		const auto c0 = -20.37260896;
		const auto c1 = 1.13944933;
		const auto c2 = -0.14210135;
		const auto c3 = 8.4644554E-3;
		const auto c4 = -1.4327641E-3;
		const auto c5 = 2.0122503E-4;
		const auto c6 = 8.6639632E-5;
		const auto c7 = -2.5850097E-5;
		const auto c8 = 2.4555012E-6;
		const auto c9 = -8.0683825E-8;
		const auto lnT = std::log(T);
		double k = c9;
		k = k * lnT + c9;
		k = k * lnT + c7;
		k = k * lnT + c6;
		k = k * lnT + c5;
		k = k * lnT + c4;
		k = k * lnT + c3;
		k = k * lnT + c2;
		k = k * lnT + c1;
		k = k * lnT + c0;
		return std::exp(k);
	} else {
		return 2.5634E-9 * std::pow(T, 1.78186);
	}

}

double k17(double T) {
	;

	T /= KperEv;
	if (T < 1.719) {
		return 2.291e-10 * std::pow(T, -0.4);
	} else {
		return 8.4258e-10 * std::pow(T, -1.4) * std::exp(-1.301 / T);
	}
}

double k18(double T) {
	;

	if (T < 617) {
		return 1e-8;
	} else {
		return 1.32e-6 * std::pow(T, -0.76);
	}
}

double k19(double T) {
	;

	return 5e-7 * std::sqrt(100 / T);
}

double compute_next_ne(std::array<double, NS> U0, std::array<double, NS> &U, double ne, double T, double dt) {
	const auto K1 = k1(T);
	const auto K2 = k2(T);
	const auto K3 = k3(T);
	const auto K4 = k4(T);
	const auto K5 = k5(T);
	const auto K6 = k6(T);
	const auto K7 = k7(T);
	const auto K8 = k8(T);
	const auto K9 = k9(T);
	const auto K10 = k10(T);
	const auto K11 = k11(T);
	const auto K12 = k12(T);
	const auto K13 = k13(T);
	const auto K14 = k14(T);
	const auto K15 = k15(T);
	const auto K16 = k16(T);
	const auto K17 = k17(T);
	const auto K18 = k18(T);
	const auto K19 = k19(T);
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
	for (double dt = 1.0; dt <= 1e+17; dt *= 1.5) {
		for (int i = 0; i < NS; i++) {
			U[i] = 0.0;
		}
		double ne = 0.1;
		U[nH] = 0.04 * ne;
		U[nHP] = 0.70 * ne;
		U[nH2P] = 0.00 * ne;
		U[nHE] = 0.25 * ne;
		U[nHEP] = 0.00 * ne;
		U[nHEPP] = 0.0 * ne;
		auto U0 = U;
		const auto eint = 1.0e-8 * ne;
		double T = 1e9;
		compute(U0, U, T, dt);
		ne = U[nHP] + U[nHEP] + 2.0 * U[nHEPP] + U[nH2P] - U[nHN];
		const auto nnuc = U[nH] + U[nHP] + U[nHN] + 2.0 * U[nH2] + 2.0 * U[nH2P] + 4.0 * U[nHE] + 4.0 * U[nHEP] + 4.0 * U[nHEPP];
		const auto nnucleus = U[nH] + U[nHP] + U[nHN] + U[nH2] + U[nH2P] + U[nHE] + U[nHEP] + U[nHEPP];
		const auto n = U[nH] + U[nHP] + U[nHN] + U[nH2] + U[nH2P] + U[nHE] + U[nHEP] + U[nHEPP] + ne;
		const auto nHnuc = U[nH] + U[nHP] + U[nHN] + 2.0 * U[nH2] + 2.0 * U[nH2P];
		const auto nHenuc = 4.0 * U[nHE] + 4.0 * U[nHEP] + 4.0 * U[nHEPP];
		printf("%14.5e %14.5e %14.5e %14.5e %14.5e  %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e\n", dt, T, nnuc, nHnuc, nHenuc, n, ne,
				U[nH], U[nHP], U[nHN], U[nH2], U[nH2P], U[nHE], U[nHEP], U[nHEPP]);
	}
}
