#include <cmath>
#include <stdio.h>
#include <vector>
#include <array>
#include <cassert>
#include <astrotiger/chemistry.hpp>

#include <functional>
#include <limits>

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

const auto tiny = 1.0e+3 * std::numeric_limits<double>::min();
const auto dhuge = std::numeric_limits<double>::max() / 1.0e+3;

void rates(double &k1, double &k2, double &k3, double &k4, double &k5, double &k6, double &k7, double &k8, double &k9, double &k10, double &k11, double &k12,
		double &k13, double &k14, double &k15, double &k16, double &k17, double &k18, double &k19, double T, bool caseb) {
	const auto tev = T / KperEv;
	const auto logtev = std::log(T);
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

const double hplanck = 6.6261e-27;
const double clight = 2.99792458e10;
const double kb = 1.3807e-16;
const double ergtoev = 6.242e+11;
const double evtoerg = 1.0 / ergtoev;

double Bp(double T) {
	return (2 * std::pow(kb * M_PI * T, 4) / (15 * clight * clight * hplanck * hplanck * hplanck));
}

double Bp_nu(double nu, double T) {
	const auto c0 = 2.0 * hplanck * std::pow(nu, 3) / (clight * clight);
	const auto x = hplanck * nu / (kb * T);
	double c1;
	if (x < 1.0e-3) {
		c1 = x;
		return c0 / c1;
	} else if (x > 1000.0) {
		return 0.0;
	} else {
		c1 = std::exp(x) - 1.0;
		return c0 / c1;
	}
}

double dB_dT(double T) {
	return 4.0 * Bp(T) / T;
}

double dBp_nu_dT(double nu, double T) {
	return hplanck * hplanck * nu * nu * nu * nu / (clight * clight * kb * T * T * (std::cosh((hplanck * nu) / (kb * T)) - 1.0));
}

double grey_cross_section(const std::function<double(double, double)> &sigma, double T) {
	const int N = 129;
	double numax = kb * T / hplanck;
	const double dnu = numax / (N - 1);
	const double lambdamax = clight * hplanck / (kb * T);
	const double dlambda = lambdamax / (N - 1);
	double sum = 0.0;
	for (int i = 1; i < N; i++) {
		const auto nu = i * dnu;
		const auto lambda = i * dlambda;
		const auto bnu = Bp_nu(nu, T) * sigma(nu, T);
		const auto blambda = Bp_nu(clight / lambda, T) * clight / (lambda * lambda) * sigma(clight / lambda, T);
		double c0;
		if (i == N - 1) {
			c0 = 1.0 / 3.0;
		} else if (i % 2 == 1) {
			c0 = 4.0 / 3.0;
		} else {
			c0 = 2.0 / 3.0;
		}
		sum += c0 * (bnu * dnu + blambda * dlambda);
	}
	return sum / Bp(T);
}

double sigma_to_rate(const std::function<double(double)> &sigma, double T) {
	const int N = 129;
	double numax = kb * T / hplanck;
	const double dnu = numax / (N - 1);
	const double lambdamax = clight * hplanck / (kb * T);
	const double dlambda = lambdamax / (N - 1);
	double sum = 0.0;
	for (int i = 1; i < N; i++) {
		const auto nu = i * dnu;
		const auto lambda = i * dlambda;
		const auto bnu = Bp_nu(nu, T) * sigma(nu) / (hplanck * nu);
		const auto blambda = Bp_nu(clight / lambda, T) * clight / (lambda * lambda) / (hplanck * nu) * sigma(clight / lambda);
		double c0;
		if (i == N - 1) {
			c0 = 1.0 / 3.0;
		} else if (i % 2 == 1) {
			c0 = 4.0 / 3.0;
		} else {
			c0 = 2.0 / 3.0;
		}
		sum += c0 * (bnu * dnu + blambda * dlambda);
	}
	return 4.0 * M_PI * sum;
}

double sigma20_22(double nu, double Z) {
	const double A0 = 6.30e-15 / (Z * Z);
	const double nuth = 13.6 * evtoerg * Z * Z / hplanck;
	if (nu > nuth) {
		const auto c0 = std::pow(nu / nuth, 4);
		const auto c1 = std::sqrt(nu / nuth - 1.0);
		const auto c2 = std::exp(4.0 - 4.0 * std::atan(c1) / c1);
		const auto c3 = 1.0 - std::exp(-2.0 * M_PI * c1);
		return A0 * c0 * c2 / c3;
	} else {
		return 0.0;
	}
}

double sigma20(double nu) {
	return sigma20_22(nu, 1.0);
}

double sigma21(double nu) {
	const double nuth = 24.6 * evtoerg / hplanck;
	if (nu > nuth) {
		return 7.42e-18 * (1.66 * std::pow(nu / nuth, -2.05) - 0.66 * std::pow(nu / nuth, -3.05));
	} else {
		return 0.0;
	}
}

double sigma22(double nu) {
	return sigma20_22(nu, 2.0);
}

double sigma_20_21_22(double nu, double H, double He, double Hep) {
	return sigma20(nu) * H + sigma21(nu) * He + sigma22(nu) * Hep;
}

double sigma23(double nu) {
	const auto nuth = 0.755 / hplanck * evtoerg;
	if (nu > nuth) {
		return 7.928e5 * std::pow(nu - nuth, 1.5) / (nu * nu * nu);
	} else {
		return 0.0;
	}
}

double sigma24(double nu) {
	const auto hnu = hplanck * nu * ergtoev;
	if (hnu < 15.42) {
		return 0.0;
	} else if (hnu < 16.50) {
		return 6.2e-18 * hnu - 9.4e-17;
	} else if (hnu < 17.7) {
		return 1.4e-18 * hnu - 1.48e-17;
	} else {
		return 2.5e-14 * std::pow(hnu, -2.71);
	}
}

double sigma25(double nu) {
	return std::pow(10,
			-1.6547717e6 + 1.8660333e5 * std::log(nu) - 7.8986431e3 * std::pow(std::log(nu), 2) + 148.73693 * std::pow(std::log(nu), 3)
					- 1.0513032 * std::pow(std::log(nu), 4));
}

double sigma26(double nu) {
	const auto hnu = hplanck * nu * ergtoev;
	if (hnu > 30.0 && hnu < 90.0) {
		return std::pow(10, -16.926 - 4.528e-2 * hnu + 2.238e-4 * hnu * hnu + 4.245e-7 * hnu * hnu * hnu);
	} else {
		return 0.0;
	}
}

double sigma28(double nu) {
	const auto hnu = nu * hplanck * ergtoev;
	double sigmaL0, sigmaW0, sigmaL1, sigmaW1;
	if (hnu > 14.675 && hnu < 16.820) {
		sigmaL0 = 1e-18 * std::pow(10, 15.1289 - 1.05139 * hnu);
	} else if (hnu >= 16.820 && hnu < 17.6) {
		sigmaL0 = 1e-18 * std::pow(10, -31.41 + 1.8042e-2 * std::pow(hnu, 3) - 4.339e-5 * std::pow(hnu, 5));
	} else {
		sigmaL0 = 0.0;
	}
	if (hnu > 14.675 && hnu < 17.7) {
		sigmaW0 = 1e-18 * std::pow(10, 13.5331 - 0.9182618 * hnu);
	} else {
		sigmaW0 = 0.0;
	}
	if (hnu > 14.159 && hnu < 15.302) {
		sigmaL1 = 1e-18 * std::pow(10, 12.0218406 - 0.819429 * hnu);
	} else if (hnu > 15.302 && hnu < 17.2) {
		sigmaL1 = 1e-18 * std::pow(10, 16.04644 - 1.082438 * hnu);
	} else {
		sigmaL1 = 0.0;
	}
	if (hnu > 14.159 && hnu < 17.2) {
		sigmaW1 = 1e-18 * std::pow(10, 12.87367 - 0.85088597 * hnu);
	} else {
		sigmaW1 = 0.0;
	}
	return 0.25 * (sigmaL0 + sigmaW0) + 0.75 * (sigmaL1 + sigmaW1);
}

double sigma_compton(double nu) {
	const auto me = 9.10938e-28;
	const auto mr = 2.81794e-13;
	const auto x = hplanck * nu / (clight * clight * me);
	double sigma;
	if (x < 0.001) {
		sigma = (8.0 * M_PI / 3.0) * (1.0 - 7.0 / 5.0 * x + 163.0 / 35.0 * x * x);
	} else {
		const auto c0 = std::atan(2.0 * std::sqrt(x));
		const auto c1 = std::log(1.0 + 4.0 * x);
		const auto c2 = M_PI * std::pow(1 - x, 2) * c0 * std::pow(x, -1.5);
		const auto c3 = M_PI * (4.0 * (1.0 + x * x + x * x * x) - (2 + x * x) * c1) / (2.0 * x);
		sigma = c2 + c3;
	}
	return mr * mr * sigma;
}

void radiation_rates(double &i20, double &i21, double &i22, double &i23, double &i24, double &i25, double &i26, double &i27, double &i28, double T) {
	i20 = sigma_to_rate(sigma20, T);
	i21 = sigma_to_rate(sigma21, T);
	i22 = sigma_to_rate(sigma22, T);
	i23 = sigma_to_rate(sigma23, T);
	i24 = sigma_to_rate(sigma24, T);
	i25 = sigma_to_rate(sigma25, T);
	i26 = sigma_to_rate(sigma26, T);
	i27 = 1.1e8 * Bp_nu(12.27 * evtoerg / hplanck, T);
	i28 = sigma_to_rate(sigma26, T);
}

double compute_next_ne(std::array<double, NS> U0, std::array<double, NS> &U, double ne, double T, double Trad, double dt) {
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
	double I20;
	double I21;
	double I22;
	double I23;
	double I24;
	double I25;
	double I26;
	double I27;
	double I28;
	rates(K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, T, false);
	radiation_rates(I20, I21, I22, I23, I24, I25, I26, I27, I28, Trad);
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
	const auto den = (1
			+ dt * (I21 + I22 + dt * I22 * K3 * ne + dt * I21 * (I22 + (K5 + K6) * ne) + ne * (K3 + K4 + K5 + K6 + dt * (K3 * K5 + (K3 + K4) * K6) * ne)));
	He = (He0 + dt * K4 * ne * (Hep0 + dt * (Hep0 + Hepp0) * K6 * ne) + dt * He0 * (I22 + ne * (K4 + K5 + K6 + dt * K4 * K6 * ne))) / den;
	Hep = (Hep0 + dt * He0 * I21 + dt * Hep0 * I21 + dt * ((He0 + Hep0) * K3 + (Hep0 + Hepp0 + dt * (He0 + Hep0 + Hepp0) * I21) * K6) * ne
			+ std::pow(dt, 2) * (He0 + Hep0 + Hepp0) * K3 * K6 * std::pow(ne, 2)) / den;
	Hepp = (dt * (Hep0 + dt * He0 * I21 + dt * Hep0 * I21 + dt * (He0 + Hep0) * K3 * ne) * (I22 + K5 * ne)
			+ Hepp0 * (1 + std::pow(dt, 2) * (I21 + K3 * ne) * (I22 + K5 * ne) + dt * (I21 + I22 + (K3 + K4 + K5) * ne))) / den;
	double err = 0.0;
#define List(a,b,c,d,e) {a,b,c,d,e}
	do {

		const std::array<std::array<double, 5>, 5> A { {
		List(1 + dt * (I20 + H2p * K10 - 2 * H2 * K13 - Hn * K15 + Hn * K8 + Hp * K9 + K1 * ne + K7 * ne),
				dt * (-(H2 * K11) - 2 * Hn * K16 + H * K9 - K2 * ne), dt * (-I23 - H * K15 - 2 * Hp * K16 - H2p * K19 + H * K8 - K14 * ne),
				dt * (-2 * I27 - 2 * I28 - Hp * K11 - 2 * H * K13 - 2 * K12 * ne), dt * (-I25 + H * K10 - Hn * K19 - 2 * K18 * ne)),
		List(dt * (-I20 - H2p * K10 + Hp * K9 - K1 * ne), 1 + dt * (H2 * K11 + Hn * K16 + Hn * K17 + H * K9 + K2 * ne), dt * (Hp * K16 + Hp * K17),
				dt * Hp * K11, dt * (-I25 - 2 * I26 - H * K10)),
		List(dt*(Hn*K15 + Hn*K8 - K7*ne),dt*(Hn*K16 + Hn*K17),1 + dt*(I23 + H*K15 + Hp*K16 + Hp*K17 + H2p*K19 + H*K8 + K14*ne),0,dt*Hn*K19),
		List(dt * (-(H2p * K10) + H2 * K13 - Hn * K8), dt * H2 * K11, dt * (-(H2p * K19) - H * K8),
				1 + dt * (I24 + I27 + I28 + Hp * K11 + H * K13 + K12 * ne), dt * (-(H * K10) - Hn * K19)),
		List(dt * (H2p * K10 - Hp * K9), dt * (-(H2 * K11) - Hn * K17 - H * K9), dt * (-(Hp * K17) + H2p * K19), dt * (-I24 - Hp * K11),
				1 + dt * (I25 + I26 + H * K10 + Hn * K19 + K18 * ne)) } };
		const std::array<double, 5> f = { H - H0
				+ dt
						* (H * I20 - Hn * I23 - H2p * I25 - 2 * H2 * I27 - 2 * H2 * I28 + H * H2p * K10 - H2 * Hp * K11 - 2 * H * H2 * K13 - H * Hn * K15
								- 2 * Hn * Hp * K16 - H2p * Hn * K19 + H * Hn * K8 + H * Hp * K9 + H * K1 * ne - 2 * H2 * K12 * ne - Hn * K14 * ne
								- 2 * H2p * K18 * ne - Hp * K2 * ne + H * K7 * ne), Hp - Hp0
				+ dt
						* (-(H * I20) - H2p * I25 - 2 * H2p * I26 - H * H2p * K10 + H2 * Hp * K11 + Hn * Hp * K16 + Hn * Hp * K17 + H * Hp * K9 - H * K1 * ne
								+ Hp * K2 * ne), Hn - Hn0
				+ dt * (Hn * I23 + H * Hn * K15 + Hn * Hp * K16 + Hn * Hp * K17 + H2p * Hn * K19 + H * Hn * K8 + Hn * K14 * ne - H * K7 * ne), H2 - H20
				+ dt * (H2 * I24 + H2 * I27 + H2 * I28 - H * H2p * K10 + H2 * Hp * K11 + H * H2 * K13 - H2p * Hn * K19 - H * Hn * K8 + H2 * K12 * ne), H2p
				- H2p0
				+ dt * (-(H2 * I24) + H2p * I25 + H2p * I26 + H * H2p * K10 - H2 * Hp * K11 - Hn * Hp * K17 + H2p * Hn * K19 - H * Hp * K9 + H2p * K18 * ne) };
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
				//			U1[n] = std::max(U1[n], 0.0);
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
		Hn = Hn1;
		H2 = H21;
		H2p = H2p1;
		err = std::abs(std::log(std::abs(Ht1 / Ht2)));
	} while (err > 1.0e-10);
	return ne;
}

double compute(const std::array<double, NS> U0, std::array<double, NS> &U, double T, double Trad, double dt) {
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
	int iter = 0;
	do {
		nemid = 0.5 * (nemax + nemin);
		auto ne = Hp - Hn + Hep + 2.0 * Hepp + H2p;
		U1 = U;
		const auto fmax = compute_next_ne(U0, U1, nemax, T, Trad, dt) - ne;
		U1 = U;
		const auto fmid = compute_next_ne(U0, U1, nemid, T, Trad, dt) - ne;
		if (fmax * fmid > 0.0) {
			nemax = nemid;
		} else {
			nemin = nemid;
		}
		iter++;

	} while (nemin / nemax < 0.999);
	U = U1;
	double ne = nemid;
	const auto evtoerg = 1.60218e-12;
	const auto Hion = -13.6 * evtoerg;
	const auto Hnion = -0.755 * evtoerg;
	const auto Heion = -(24.6 + 13.6 * 4) * evtoerg;
	const auto Hepion = -24.6 * evtoerg;
	const auto H2ion = -15.42 * evtoerg;
	const auto kb = 1.38e-16;
	double n = U[nH] + U[nHP] + U[nHN] + U[nH2] + U[nH2P] + U[nHE] + U[nHEP] + U[nHEPP] + ne;
	double eion = Hion * U[nH] + Hnion * U[nHN] + Heion * U[nHE] + Hepion * U[nHEP] + H2ion * U[nH2];
	double nmon = U[nH] + U[nHP] + U[nHN] + U[nHE] + U[nHEP] + U[nHEPP] + ne;
	double ndia = U[nH2] + U[nH2P];
	return (1.5 * nmon + 2.5 * ndia) * kb * T + eion;
}

double compute_next_chemistry(std::array<double, NS> U0, std::array<double, NS> &U, double energy, double Trad, double dt) {
	const auto evtoerg = 1.60218e-12;
	const auto Hion = -13.6 * evtoerg;
	const auto Hnion = -0.755 * evtoerg;
	const auto Heion = -(24.6 + 13.6 * 4) * evtoerg;
	const auto Hepion = -24.6 * evtoerg;
	const auto H2ion = -15.42 * evtoerg;
	const auto kb = 1.38e-16;
	double ne = U[nHP] - U[nHN] + U[nH2P] + U[nHEP] + 2.0 * U[nHEPP];
	double n = U[nH] + U[nHP] + U[nHN] + U[nH2] + U[nH2P] + U[nHE] + U[nHEP] + U[nHEPP] + ne;
	double eion = Hion * U[nH] + Hnion * U[nHN] + Heion * U[nHE] + Hepion * U[nHEP] + H2ion * U[nH2];
	double nmon = U[nH] + U[nHP] + U[nHN] + U[nHE] + U[nHEP] + U[nHEPP] + ne;
	double ndia = U[nH2] + U[nH2P];
	double T = (energy - eion) / kb / (nmon * 1.5 + ndia * 2.5);
	double T0;
	auto U1 = U;
	do {
		T0 = T;
		const auto dT = T * 0.0001;
		U1 = U;
		const auto f1 = compute(U0, U1, T, Trad, dt) - energy;
		U1 = U;
		const auto f2 = compute(U0, U1, T + dT, Trad, dt) - energy;
		const auto dfdT = (f2 - f1) / dT;
		T -= (f1 / dfdT);
	} while (std::abs(std::log(T / T0)) > 1.0e-3);
	U = U1;
	return T;
}

double ion_energy(species s) {
	const auto evtoerg = 1.60218e-12;
	const auto Hion = -13.6 * evtoerg;
	const auto Hnion = -0.755 * evtoerg;
	const auto Heion = -(24.6 + 13.6 * 4) * evtoerg;
	const auto Hepion = -24.6 * evtoerg;
	const auto H2ion = -15.42 * evtoerg;
	return s.H * Hion + s.Hn * Hnion + s.He * Heion + s.Hep * Hepion + s.H2 * H2ion;
}

thermo_props compute_thermo_properties(const species s, double energy) {
	thermo_props p;
	const double ne = s.Hp - s.Hn + s.H2p + s.Hep + 2 * s.Hepp;
	const double nA = s.H + s.Hp + s.Hn + 2.0 * s.H2 + 2.0 * s.H2p + 4.0 * s.He + 4.0 * s.Hep + 4.0 * s.Hepp;
	p.rho = nA * 1.6605e-24;
	const double nmon = ne + s.H + s.Hp + s.Hn + s.He + s.Hep + s.Hepp;
	const double ndia = s.H2 + s.H2p;
	const double n = nmon + ndia;
	const double kb = 1.38e-16;
	p.eion = ion_energy(s);
	p.gamma = (5 * nmon + 7 * ndia) / (3 * nmon + 5 * ndia);
	p.pressure = (p.gamma - 1.0) * std::max(energy - p.eion, 0.0);
	p.sound_speed = std::sqrt(p.gamma * p.pressure / p.rho);
	p.T = p.pressure / (kb * n);
	return p;
}

species compute_next_species(const species s0, double energy, double Trad, double dt) {
	std::array<double, NS> U, U0;
	species s;
	U0[nH] = s0.H;
	U0[nHP] = s0.Hp;
	U0[nHN] = s0.Hn;
	U0[nH2] = s0.H2;
	U0[nH2P] = s0.H2p;
	U0[nHE] = s0.He;
	U0[nHEP] = s0.Hep;
	U0[nHEPP] = s0.Hepp;
	compute_next_chemistry(U0, U, energy, Trad, dt);
	s.H = U[nH];
	s.Hp = U[nHP];
	s.Hn = U[nHN];
	s.H2 = U[nH2];
	s.H2p = U[nH2P];
	s.He = U[nHE];
	s.Hep = U[nHEP];
	s.Hepp = U[nHEPP];
	return s;
}

double cooling_rate(species s, double T) {
	const auto ne = s.Hp - s.Hn + s.Hep + 2.0 * s.Hepp + s.H2p;
	const auto T3 = T / 1e3;
	const auto T5 = T / 1e5;
	const auto T6 = T / 1e6;
	const auto tmp = (1.0 + std::sqrt(T5));
	const auto tev = T / KperEv;
	const auto logtev = std::log(T);
	const auto tiny = std::numeric_limits<double>::min();
	const auto k1 = exp(
			-32.71396786375 + 13.53655609057 * logtev - 5.739328757388 * std::pow(logtev, 2) + 1.563154982022 * std::pow(logtev, 3)
					- 0.2877056004391 * std::pow(logtev, 4) + 0.03482559773736999 * std::pow(logtev, 5) - 0.00263197617559 * std::pow(logtev, 6)
					+ 0.0001119543953861 * std::pow(logtev, 7) - 2.039149852002e-6 * std::pow(logtev, 8));

	const auto k3 = exp(
			-44.09864886561001 + 23.91596563469 * logtev - 10.75323019821 * std::pow(logtev, 2) + 3.058038757198 * std::pow(logtev, 3)
					- 0.5685118909884001 * std::pow(logtev, 4) + 0.06795391233790001 * std::pow(logtev, 5) - 0.005009056101857001 * std::pow(logtev, 6)
					+ 0.0002067236157507 * std::pow(logtev, 7) - 3.649161410833d - 6 * std::pow(logtev, 8));

	const auto k5 = exp(
			-68.71040990212001 + 43.93347632635 * logtev - 18.48066993568 * std::pow(logtev, 2) + 4.701626486759002 * std::pow(logtev, 3)
					- 0.7692466334492 * std::pow(logtev, 4) + 0.08113042097303 * std::pow(logtev, 5) - 0.005324020628287001 * std::pow(logtev, 6)
					+ 0.0001975705312221 * std::pow(logtev, 7) - 3.165581065665e-6 * std::pow(logtev, 8));
	const auto c0 = 7.50e-19 * tmp * std::exp(-118348 / T) * s.H * ne;
	const auto c1 = 9.10e-27 * tmp * std::pow(T, -.1687) * std::exp(-13179 / T) * ne * ne * s.He;
	const auto c2 = 5.54e-17 * tmp * std::pow(T, -.397) * std::exp(-473638 / T) * ne * s.Hep;
	const auto c3 = 2.18e-11 * k1 * ne * s.H;
	const auto c4 = 3.94e-11 * k3 * ne * s.He;
	const auto c5 = 8.72e-11 * k5 * ne * s.Hep;
	const auto c6 = 5.01e-27 * tmp * std::pow(T, -.1687) * std::exp(-55338 / T) * ne * ne * s.Hep;
	const auto c7 = 8.70e-27 * std::sqrt(T) * std::pow(T3, -0.2) / (1.0 + std::pow(T6, 0.7)) * ne * s.Hp;
	const auto c8 = 1.55e-26 * std::pow(T, 0.3647) * ne * s.Hep;
	const auto c9 = 1.24e-13 * std::pow(T, -1.5) * (1 + 0.3 * std::exp(-94000 / T)) * std::exp(-470000 / T) * ne * s.Hep;
	const auto c10 = 3.48e-26 * std::sqrt(T) * std::pow(T3, -0.2) / (1.0 + std::pow(T6, 0.7)) * ne * s.Hepp;
	const auto c11 = 1.43e-27 * std::sqrt(T) * (1.1 + 0.34 * std::exp(-std::pow(5.50 - std::log10(T), 2) / 3.0)) * ne * (s.Hp + s.Hep + s.Hepp);
	const auto z = 0.0;
	const auto c12 = 5.64e-36 * std::pow(1 + z, 4) * (T - 2.73 * (1 + z)) * ne;
	const auto xx = std::log10(T / 1e4);
	double vibha = 1.1e-18 * std::exp(-std::min(std::log(dhuge), 6744 / T));
	double dum;
	dum = 8.152e-13 * (4.2 / (kb * (T + 1190.)) + 1. / (kb * T));
	double h2k01 = 1.45e-12 * std::sqrt(T) * std::exp(-std::min(log(dhuge), dum));
	if (T > 1635.) {
		dum = 1.0e-12 * std::sqrt(T) * std::exp(-1000. / T);
	} else {
		dum = 1.4e-13 * std::exp((T / 125.) - std::pow(T / 577., 2));
	}
	double hyd01 = dum * exp(-std::min(std::log(dhuge), 8.152e-13 / (kb * T)));
	double vibla = hyd01 * s.H + h2k01 * s.H2;
	double rotla, rotha;
	if (T > 1087) {
		rotha = 3.9e-19 * std::exp(-6118 / T);
	} else {
		rotha = std::pow(10.0, (-19.24 + 0.474 * xx - 1.247 * xx * xx));
	}
	if (T > 4031) {
		rotla = 1.38e-22 * exp(-9243 / T);
	} else {
		rotla = std::pow(10, (-22.9 - 0.553 * xx - 1.148 * xx * xx));
	}
	rotla *= (std::pow(s.H2, 0.77) + 1.2 * std::pow(s.H, 0.77));
	const auto c13 = s.H2 * ((vibla * vibha) / (vibha + vibla) + (rotla * rotha) / (rotha + rotla));
	const auto c14 = sigma_to_rate([s](double nu) {
		return (nu, s.H, s.He, s.Hep);
	}, T);
	return c0 + c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10 + c11 + c12 + c13 + c14;
}

double energy_mean_opacity(species s, double T) {
	const int N = 129;
	double numax = kb * T / hplanck;
	const double dnu = numax / (N - 1);
	const double lambdamax = clight * hplanck / (kb * T);
	const double dlambda = lambdamax / (N - 1);
	double sum = 0.0;
	for (int i = 1; i < N; i++) {
		const auto nu = i * dnu;
		const auto lambda = i * dlambda;
		const auto sigma_nu = s.H * sigma20(nu) + s.He * sigma21(nu) + s.Hep * sigma22(nu) + s.Hn * sigma23(nu) + s.H2 * sigma24(nu) + s.H2p * sigma25(nu)
				+ s.H2p * sigma26(nu) + s.H2 * sigma28(nu);
		const auto nu2 = clight / lambda;
		const auto sigma_lambda = s.H * sigma20(nu2) + s.He * sigma21(nu2) + s.Hep * sigma22(nu2) + s.Hn * sigma23(nu2) + s.H2 * sigma24(nu2)
				+ s.H2p * sigma25(nu2) + s.H2p * sigma26(nu2) + s.H2 * sigma28(nu2);
		const auto bnu = Bp_nu(nu, T) * sigma_nu;
		const auto blambda = Bp_nu(nu2, T) * clight / (lambda * lambda) * sigma_lambda;
		double c0;
		if (i == N - 1) {
			c0 = 1.0 / 3.0;
		} else if (i % 2 == 1) {
			c0 = 4.0 / 3.0;
		} else {
			c0 = 2.0 / 3.0;
		}
		sum += c0 * (bnu * dnu + blambda * dlambda);
	}
	return sum / Bp(T);

}

double rossalind_mean_opacity(species s, double T) {
	const int N = 129;
	double numax = kb * T / hplanck;
	const double dnu = numax / (N - 1);
	const double lambdamax = clight * hplanck / (kb * T);
	const double dlambda = lambdamax / (N - 1);
	double sum = 0.0;
	auto ne = s.Hp - s.Hn + s.Hep + 2.0 * s.Hepp + s.H2p;
	for (int i = 1; i < N; i++) {
		const auto nu = i * dnu;
		const auto lambda = i * dlambda;
		const auto sigma_nu = s.H * sigma20(nu) + s.He * sigma21(nu) + s.Hep * sigma22(nu) + s.Hn * sigma23(nu) + s.H2 * sigma24(nu) + s.H2p * sigma25(nu)
				+ s.H2p * sigma26(nu) + s.H2 * sigma28(nu) + ne * sigma_compton(nu);
		const auto nu2 = clight / lambda;
		const auto sigma_lambda = s.H * sigma20(nu2) + s.He * sigma21(nu2) + s.Hep * sigma22(nu2) + s.Hn * sigma23(nu2) + s.H2 * sigma24(nu2)
				+ s.H2p * sigma25(nu2) + s.H2p * sigma26(nu2) + s.H2 * sigma28(nu2) + ne * sigma_compton(nu2);
		const auto bnu = dBp_nu_dT(nu, T) / sigma_nu;
		const auto blambda = dBp_nu_dT(nu2, T) * clight / (lambda * lambda) / sigma_lambda;
		double c0;
		if (i == N - 1) {
			c0 = 1.0 / 3.0;
		} else if (i % 2 == 1) {
			c0 = 4.0 / 3.0;
		} else {
			c0 = 2.0 / 3.0;
		}
		sum += c0 * (bnu * dnu + blambda * dlambda);
	}
	return dB_dT(T) / sum;

}

double flux_mean_opacity(species s, double Tg, double Tr, double dTgdx) {
	const int N = 129;
	const double arad = 7.5646e-15;
	const auto Tmean = std::sqrt(Tg * Tr);
	double numax = kb * Tmean / hplanck;
	const double dnu = numax / (N - 1);
	const double lambdamax = clight * hplanck / (kb * Tmean);
	const double dlambda = lambdamax / (N - 1);
	double sum1 = 0.0;
	double sum2 = 0.0;
	auto ne = s.Hp - s.Hn + s.Hep + 2.0 * s.Hepp + s.H2p;
	for (int i = 1; i < N; i++) {
		const auto lambda = i * dlambda;
		const auto nu1 = i * dnu;
		const auto nu2 = clight / lambda;
		const auto dBp_nu_dT1 = dBp_nu_dT(nu1, Tg);
		const auto dBp_nu_dT2 = dBp_nu_dT(nu2, Tg);
		auto sigma1 = s.H * sigma20(nu1) + s.He * sigma21(nu1) + s.Hep * sigma22(nu1) + s.Hn * sigma23(nu1) + s.H2 * sigma24(nu1) + s.H2p * sigma25(nu1)
				+ s.H2p * sigma26(nu1) + s.H2 * sigma28(nu1) + ne * sigma_compton(nu1);
		auto sigma2 = s.H * sigma20(nu2) + s.He * sigma21(nu2) + s.Hep * sigma22(nu2) + s.Hn * sigma23(nu2) + s.H2 * sigma24(nu2) + s.H2p * sigma25(nu2)
				+ s.H2p * sigma26(nu2) + s.H2 * sigma28(nu2) + ne * sigma_compton(nu2);
		const auto den1 = sigma1 * Bp_nu(nu1, Tr);
		const auto den2 = sigma2 * Bp_nu(nu2, Tr);
		const auto num1 = dBp_nu_dT1 * dTgdx;
		const auto num2 = dBp_nu_dT2 * dTgdx;
		const auto Rnu1 = std::min(num1 / den1, 1e+10);
		const auto Rnu2 = std::min(num2 / den2, 1e+10);
		const auto lnu1 = (2.0 + Rnu1) / (6.0 + 2.0 * Rnu1 + Rnu1 * Rnu1);
		const auto lnu2 = (2.0 + Rnu2) / (6.0 + 2.0 * Rnu2 + Rnu2 * Rnu2);
		const auto b1num = lnu1 * dBp_nu_dT1 * dnu;
		const auto b2num = lnu2 * dBp_nu_dT2 * clight / (lambda * lambda) * dlambda;
		const auto b1den = lnu1 * dBp_nu_dT1 / sigma1 * dnu;
		const auto b2den = lnu2 * dBp_nu_dT2 * clight / (lambda * lambda) / sigma2 * dlambda;
		double c0;
		if (i == N - 1) {
			c0 = 1.0 / 3.0;
		} else if (i % 2 == 1) {
			c0 = 4.0 / 3.0;
		} else {
			c0 = 2.0 / 3.0;
		}
		sum1 += c0 * (b1num + b2num);
		sum2 += c0 * (b1den + b2den);
	}
	return sum1 / sum2;

}

void chemistry_test() {
	species s;
	const auto n = 1.0e+10;
	s.H = 1e-10 * n;
	s.Hp = 1e-10 * n;
	s.Hn = 1e-10 * n;
	s.H2 = 0.92 * n;
	s.H2p = 1e-10 * n;
	s.He = 0.08 * n;
	s.Hep = 1e-10 * n;
	s.Hepp = 1e-10 * n;
	const auto dTdx = 1.0;
	for (double T = 1.0e1; T <= 1e10; T *= 1.5) {
		const auto kr = rossalind_mean_opacity(s, T);
		const auto kp = energy_mean_opacity(s, T);
		const auto ke = energy_mean_opacity(s, T);
		const auto kf = flux_mean_opacity(s, T, T, dTdx);
		const auto R = dTdx * dB_dT(T) / kr / Bp(T);
		const auto lambda = (2.0 + R) / (6.0 + 2.0 * R + R * R);
		printf("%e %e %e %e %e %e\n", T, kr, kf, kp, ke, lambda);
	}
}

//void chemistry_test() {
//	std::array<double, NS> U;
//	const auto dt = 1.0;
//	printf("%14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s\n", "time", "T", "A", "AH", "AHe", "N", "Ne", "H", "H+", "H-", "H2",
//			"H2p", "He", "He+", "He++");
//	double T = 1.0e4;
//	for (double dt = 1.0e+3; dt <= 1e+17; dt *= 10) {
//		for (int i = 0; i < NS; i++) {
//			U[i] = 0.0;
//		}
//		double ne = 1e-1;
//		U[nH] = 0.01 * ne;
//		U[nHP] = 0.91 * ne;
//		U[nH2P] = 0.00 * ne;
//		U[nHE] = 0.00 * ne;
//		U[nHEP] = 0.08 * ne;
//		U[nHEPP] = 0.0 * ne;
//		auto U0 = U;
//		const auto eint = 1.0e-10 * ne;
//		T = compute_thermo_properties(U0, U, eint, dt);
//		ne = U[nHP] + U[nHEP] + 2.0 * U[nHEPP] + U[nH2P] - U[nHN];
//		const auto nnuc = U[nH] + U[nHP] + U[nHN] + 2.0 * U[nH2] + 2.0 * U[nH2P] + 4.0 * U[nHE] + 4.0 * U[nHEP] + 4.0 * U[nHEPP];
//		const auto nnucleus = U[nH] + U[nHP] + U[nHN] + U[nH2] + U[nH2P] + U[nHE] + U[nHEP] + U[nHEPP];
//		const auto n = U[nH] + U[nHP] + U[nHN] + U[nH2] + U[nH2P] + U[nHE] + U[nHEP] + U[nHEPP] + ne;
//		const auto nHnuc = U[nH] + U[nHP] + U[nHN] + 2.0 * U[nH2] + 2.0 * U[nH2P];
//		const auto nHenuc = 4.0 * U[nHE] + 4.0 * U[nHEP] + 4.0 * U[nHEPP];
//		printf("%14.5e %14.5e %14.5e %14.5e %14.5e  %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e\n", (double) dt, (double) T,
//				(double) nnuc, (double) nHnuc, (double) nHenuc, (double) n, (double) ne, (double) U[nH], (double) U[nHP], (double) U[nHN], (double) U[nH2],
//				(double) U[nH2P], (double) U[nHE], (double) U[nHEP], (double) U[nHEPP]);
//		U = U0;
//	}
//}
