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

void chemical_rates(double &k1, double &k2, double &k3, double &k4, double &k5, double &k6, double &k7, double &k8, double &k9, double &k10, double &k11,
		double &k12, double &k13, double &k14, double &k15, double &k16, double &k17, double &k18, double &k19, double T) {

	const auto tev = T / KperEv;
	const auto logtev = std::log(tev);
	k1 = exp(
			-32.71396786375 + 13.53655609057 * logtev - 5.739328757388 * std::pow(logtev, 2) + 1.563154982022 * std::pow(logtev, 3)
					- 0.2877056004391 * std::pow(logtev, 4) + 0.03482559773736999 * std::pow(logtev, 5) - 0.00263197617559 * std::pow(logtev, 6)
					+ 0.0001119543953861 * std::pow(logtev, 7) - 2.039149852002e-6 * std::pow(logtev, 8));

	k3 = exp(
			-44.09864886561001 + 23.91596563469 * logtev - 10.75323019821 * std::pow(logtev, 2) + 3.058038757198 * std::pow(logtev, 3)
					- 0.5685118909884001 * std::pow(logtev, 4) + 0.06795391233790001 * std::pow(logtev, 5) - 0.005009056101857001 * std::pow(logtev, 6)
					+ 0.0002067236157507 * std::pow(logtev, 7) - 3.649161410833e-6 * std::pow(logtev, 8));

	k4 = 3.92e-13 / std::pow(tev, 0.6353);

	if (tev > 0.1) {
		k4 += 1.54e-9 * (1. + 0.3 / std::exp(8.099328789667 / tev)) / (exp(40.49664394833662 / tev) * std::pow(tev, 1.5));
	}

	k5 = exp(
			-68.71040990212001 + 43.93347632635 * logtev - 18.48066993568 * std::pow(logtev, 2) + 4.701626486759002 * std::pow(logtev, 3)
					- 0.7692466334492 * std::pow(logtev, 4) + 0.08113042097303 * std::pow(logtev, 5) - 0.005324020628287001 * std::pow(logtev, 6)
					+ 0.0001975705312221 * std::pow(logtev, 7) - 3.165581065665e-6 * std::pow(logtev, 8));

	k2 = exp(
			-28.61303380689232 - 0.7241125657826851 * logtev - 0.02026044731984691 * std::pow(logtev, 2) - 0.002380861877349834 * std::pow(logtev, 3)
					- 0.0003212605213188796 * std::pow(logtev, 4) - 0.00001421502914054107 * std::pow(logtev, 5) + 4.989108920299513e-6 * std::pow(logtev, 6)
					+ 5.755614137575758e-7 * std::pow(logtev, 7) - 1.856767039775261e-8 * std::pow(logtev, 8) - 3.071135243196595e-9 * std::pow(logtev, 9));

	k6 = 3.36e-10 / sqrt(T) / std::pow((T / 1.e3), 0.2) / (1. + std::pow((T / 1.e6), 0.7));

	if (T < 3275.2) {
		k7 = 1.429e-18 * std::pow(T, 0.7620) * std::pow(T, 0.1523 * std::log10(T)) * std::pow(T, -3.274e-2 * std::pow(std::log10(T), 2));
	} else {
		k7 = 3.802e-17 * std::pow(T, 0.1998 * std::log10(T)) * std::pow(10.0, 4.0415e-5 * std::pow(std::log10(T), 6) - 5.447e-3 * std::pow(std::log10(T), 4));
	}

	auto this_tev = std::max(0.1, tev);
	auto this_logtev = std::log(this_tev);
	k8 = exp(
			-20.06913897587003 + 0.2289800603272916 * this_logtev + 0.03599837721023835 * std::pow(this_logtev, 2)
					- 0.004555120027032095 * std::pow(this_logtev, 3) - 0.0003105115447124016 * std::pow(this_logtev, 4)
					+ 0.0001073294010367247 * std::pow(this_logtev, 5) - 8.36671960467864e-6 * std::pow(this_logtev, 6)
					+ 2.238306228891639e-7 * std::pow(this_logtev, 7));

	if (T > 7843.1) {
		k9 = 5.81e-16 * std::pow(T / 56200., (-0.6657 * std::log10(T / 56200.)));
	} else {
		k9 = 1.85e-23 * std::pow(T, 1.8);
	}

	k10 = 6.0e-10;

	this_tev = std::max(0.3, tev);
	this_logtev = std::log(this_tev);
	k13 = 1.0670825e-10 * std::pow(this_tev, 2.012) / (exp(4.463 / this_tev) * std::pow((1. + 0.2472 * this_tev), 3.512));
	k11 = exp(
			-24.24914687731536 + 3.400824447095291 * this_logtev - 3.898003964650152 * std::pow(this_logtev, 2) + 2.045587822403071 * std::pow(this_logtev, 3)
					- 0.5416182856220388 * std::pow(this_logtev, 4) + 0.0841077503763412 * std::pow(this_logtev, 5)
					- 0.007879026154483455 * std::pow(this_logtev, 6) + 0.0004138398421504563 * std::pow(this_logtev, 7)
					- 9.36345888928611e-6 * std::pow(this_logtev, 8));
	k12 = 5.6e-11 * exp(-102124. / T) * std::pow(T, 0.5);

	this_tev = std::max(0.04, tev);
	this_logtev = std::log(this_tev);
	k14 = exp(
			-18.01849334273 + 2.360852208681 * this_logtev - 0.2827443061704 * std::pow(this_logtev, 2) + 0.01623316639567 * std::pow(this_logtev, 3)
					- 0.03365012031362999 * std::pow(this_logtev, 4) + 0.01178329782711 * std::pow(this_logtev, 5)
					- 0.001656194699504 * std::pow(this_logtev, 6) + 0.0001068275202678 * std::pow(this_logtev, 7)
					- 2.631285809207e-6 * std::pow(this_logtev, 8));
	if (tev > 0.0979824) {
		k15 = exp(
				-20.37260896533324 + 1.139449335841631 * logtev - 0.1421013521554148 * std::pow(logtev, 2) + 0.00846445538663 * std::pow(logtev, 3)
						- 0.0014327641212992 * std::pow(logtev, 4) + 0.0002012250284791 * std::pow(logtev, 5) + 0.0000866396324309 * std::pow(logtev, 6)
						- 0.00002585009680264 * std::pow(logtev, 7) + 2.4555011970392e-6 * std::pow(logtev, 8) - 8.06838246118e-8 * std::pow(logtev, 9));
	} else {
		k15 = 2.56e-9 * std::pow(tev, 1.78186);
	}

	k16 = 6.5e-9 / sqrt(tev);

	if (tev < 1.74498) {
		k17 = 2.291e-10 * std::pow(tev, -0.4);
	} else {
		k17 = 8.4258e-10 * std::pow(tev, -1.4) * std::exp(-1.301 / tev);
	}

	const auto this_T = std::max(616.92, T);
	k18 = 1.32e-6 * std::pow(this_T, -0.76);

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
	} else if (x > 100.0) {
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

double radiation_heating_rate(const std::function<double(double)> &sigma, double &J, double hnuth, double T) {
	const int N = 129;
	double numax = kb * T / hplanck;
	const double dnu = numax / (N - 1);
	const double lambdamax = clight * hplanck / (kb * T);
	const double dlambda = lambdamax / (N - 1);
	for (int i = 1; i < N; i++) {
		const auto nu = i * dnu;
		const auto lambda = i * dlambda;
		const auto sigma_nu = sigma(nu);
		const auto sigma_lambda = sigma(clight / lambda);
		const auto c0_nu = 1.0 / (hplanck * nu) * std::max(1.0 - hplanck * nu / hnuth, 0.0);
		const auto c0_lambda = 1.0 / (hplanck * lambda) * std::max(1.0 - hplanck * clight / hnuth / lambda, 0.0);
		const auto bnu = Bp_nu(nu, T) * sigma_nu * c0_nu;
		const auto blambda = Bp_nu(clight / lambda, T) * c0_lambda * sigma_lambda;
		double c0;
		if (i == N - 1) {
			c0 = 1.0 / 3.0;
		} else if (i % 2 == 1) {
			c0 = 4.0 / 3.0;
		} else {
			c0 = 2.0 / 3.0;
		}
		J += c0 * (bnu * dnu + blambda * dlambda);
	}
	J *= 4.0 * M_PI / clight;

}

void sigma_to_rate(const std::function<double(double)> &sigma, double &I, double T) {
	const int N = 129;
	double numax = kb * T / hplanck;
	const double dnu = numax / (N - 1);
	const double lambdamax = clight * hplanck / (kb * T);
	const double dlambda = lambdamax / (N - 1);
	I = 0.0;
	for (int i = 1; i < N; i++) {
		const auto nu = i * dnu;
		const auto lambda = i * dlambda;
		const auto sigma_nu = sigma(nu);
		const auto sigma_lambda = sigma(clight / lambda);
		const auto c0_nu = 1.0 / (hplanck * nu);
		const auto c0_lambda = 1.0 / (hplanck * lambda);
		const auto bnu = Bp_nu(nu, T) * sigma_nu * c0_nu;
		const auto blambda = Bp_nu(clight / lambda, T) * c0_lambda * sigma_lambda;
		double c0;
		if (i == N - 1) {
			c0 = 1.0 / 3.0;
		} else if (i % 2 == 1) {
			c0 = 4.0 / 3.0;
		} else {
			c0 = 2.0 / 3.0;
		}
		I += c0 * (bnu * dnu + blambda * dlambda);
	}
	I *= 4.0 * M_PI;
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

void heating_rates(double &J20, double &J21, double &J22, double T) {
	radiation_heating_rate(sigma20, J20, 13.6 * evtoerg, T);
	radiation_heating_rate(sigma21, J21, 13.6 * 4.0 * evtoerg, T);
	radiation_heating_rate(sigma22, J22, 24.6 * evtoerg, T);
}

double power_law_spectrum(double nu) {
	const double hnu = hplanck * nu;
//	return 1e-21 * (12.86 / hnu) * std::exp(-1e22 * (sigma20(nu) + 0.08 * sigma21(nu)));
	return 3.9e-40 * std::pow(hnu / 13.6, -1.5);
}

void sigma_to_rate(const std::function<double(double)> &spec, const std::function<double(double)> &sigma, double &I) {
	const double numin = std::pow(10, 3);
	const double numax = std::pow(10, 20);
	I = 0.0;
	const int N = 1024;
	const double ymin = std::log(numin);
	const double ymax = std::log(numax);
	const auto dy = (ymax - ymin) / (N - 1);
	for (int i = 1; i <= N; i++) {
		double y = i * dy;
		double nu = std::exp(y);
		double c0;
		if (i == N) {
			c0 = 1.0 / 3.0;
		} else if (i % 2 == 1) {
			c0 = 4.0 / 3.0;
		} else {
			c0 = 2.0 / 3.0;
		}
		I += 4.0 * M_PI * c0 * spec(nu) * sigma(nu) / hplanck * dy;
	}
}

double homo_rate(double z, double alpha, double beta, double gamma, double z0, double z1) {
	auto c0 = beta * std::pow(z - z0, 2) / (1.0 + gamma * std::pow(z + z1, 2));
	c0 = std::min(c0, 100.0);
	c0 = std::max(c0, -100.0);
	return std::pow(1 + z, alpha) * std::exp(c0);
}

void radiation_rates(double &I20, double &I21, double &I22, double &I23, double &I24, double &I25, double &I26, double &I27, double &I28, double z) {
	I20 = 1.04e-12 * homo_rate(z, 0.231, -0.6818, 0.1646, 1.855, 0.3097);
	I21 = 1.84e-14 * homo_rate(z, -1.038, -1.1640, 0.1940, 1.973, -0.6561);
	I22 = 5.79e-13 * homo_rate(z, 0.278, -0.8260, 0.1730, 1.973, 0.2880);
	sigma_to_rate(power_law_spectrum, sigma23, I23);
	sigma_to_rate(power_law_spectrum, sigma24, I24);
	sigma_to_rate(power_law_spectrum, sigma25, I25);
	sigma_to_rate(power_law_spectrum, sigma26, I26);
	sigma_to_rate(power_law_spectrum, sigma28, I28);
	I27 = 1.1e8 * power_law_spectrum(12.27 * evtoerg / hplanck);
}

double ion_energy(species s) {
	const auto evtoerg = 1.60218e-12;
	const auto Hion = -13.6 * evtoerg;
	const auto Heion = -(24.6 + 13.6 * 4) * evtoerg;
	const auto Hepion = -24.6 * evtoerg;
	return s.H * Hion + s.He * Heion + s.Hep * Hepion;
}

void cooling_rate2(double &C14, double &dC14dT, double &dC14dH, double &dC14dH2, double H, double H2, double T) {

	/*********************************/
	C14 = 0.0;
	dC14dT = 0.0;
	dC14dH = 0.0;
	dC14dH2 = 0.0;
	return;
	/*********************************/

	const auto xx = std::log10(T / 1e4);
	const auto dxx_dT = 1 / (T * std::log(10));
	double vibha = 1.1e-18 * std::exp(-6744 / T);
	double dvibha_dT = 1.1e-18 * 6744 / (std::exp(6744 / T) * std::pow(T, 2));
	double dum1, dum2, ddum2_dT;
	dum1 = 8.152e-13 * (4.2 / (kb * (T + 1190.)) + 1. / (kb * T));
	double ddum1_dT = 8.152e-13 * (-1. / std::pow(T, 2) - 4.2 / std::pow(1190. + T, 2)) / kb;
	double h2k01 = 1.45e-12 * std::sqrt(T) * std::exp(dum1);
	double dh2k01_dT = 1.45e-12 * (std::exp(dum1) * (1 + 2 * T * ddum1_dT)) / (2. * std::sqrt(T));
	if (T > 1635.) {
		dum2 = 1.0e-12 * std::sqrt(T) * std::exp(-1000. / T);
		ddum2_dT = 1.0e-12 * (0.5 * (2000. + T)) / (std::exp(1000. / T) * std::pow(T, 1.5));
	} else {
		dum2 = 1.4e-13 * std::exp((T / 125.) - std::pow(T / 577., 2));
		ddum2_dT = 1.4e-13 * std::exp((0.008 - 3.003643419467814e-6 * T) * T) * (0.008 - 6.007286838935628e-6 * T);
	}
	double hyd01 = dum2 * std::exp(8.152e-13 / (kb * T));
	double dhyd01_dT = ddum2_dT * hyd01 + dum2 * (-8.152e-13 * std::exp(8.152e-13 / (kb * T))) / (kb * std::pow(T, 2));
	double vibla = hyd01 * H + h2k01 * H2;
	double dvibla_dT = dhyd01_dT * H + dh2k01_dT * H2;
	double dvibla_dH = hyd01;
	double dvibla_dH2 = h2k01;
	double rotla0, rotha, drotla0_dT, drotha_dT;
	if (T > 1087) {
		rotha = 3.9e-19 * std::exp(-6118 / T);
		drotha_dT = 3.9e-19 * 6118 / (std::exp(6118 / T) * std::pow(T, 2));
	} else {
		rotha = std::pow(10.0, (-19.24 + 0.474 * xx - 1.247 * xx * xx));
		drotha_dT = (std::pow(T, 9.45) * (8.533285780697259e-41 - 8.844652214438655e-42 * std::log(T)))
				/ std::exp(0.541565218933355 * std::pow(std::log(T), 2));
	}
	if (T > 4031) {
		rotla0 = 1.38e-22 * exp(-9243 / T);
		drotla0_dT = 9243 * rotla0 / (T * T);
	} else {
		rotla0 = std::pow(10, (-22.9 - 0.553 * xx - 1.148 * xx * xx));
		drotla0_dT = (-6.961857527221768e-24 - 1.2553250493430911e-23 * std::log(T))
				/ (std::exp(0.49857006522493297 * std::pow(std::log(T), 2)) * std::pow(T, 1.553));
	}
	const auto rotla = rotla0 * (std::pow(H2, 0.77) + 1.2 * std::pow(H, 0.77));
	const auto drotla_dT = drotla0_dT * (std::pow(H2, 0.77) + 1.2 * std::pow(H, 0.77));
	const auto drotla_dH2 = rotla0 * 0.77 * std::pow(H2, 0.77 - 1.0);
	const auto drotla_dH = rotla0 * 1.2 * 0.77 * std::pow(H, 0.77 - 1.0);

	C14 = H2 * ((vibla * vibha) / (vibha + vibla) + (rotla * rotha) / (rotha + rotla));
	dC14dT = (H2
			* (rotha * (-(rotla * drotha_dT) + rotha * drotla_dT)
					+ ((rotha + rotla)
							* (rotha * (std::pow(vibla, 2) * dvibha_dT + std::pow(vibha, 2) * dvibla_dT)
									+ rotla
											* (2 * vibha * vibla * drotha_dT + std::pow(vibla, 2) * (drotha_dT + dvibha_dT)
													+ std::pow(vibha, 2) * (drotha_dT + dvibla_dT)))) / std::pow(vibha + vibla, 2)))
			/ std::pow(rotha + rotla, 2);
	dC14dH = H2 * ((std::pow(rotha, 2) * drotla_dH) / std::pow(rotha + rotla, 2) + (std::pow(vibha, 2) * dvibla_dH) / std::pow(vibha + vibla, 2));
	dC14dH2 = C14 + H2 * ((std::pow(rotha, 2) * drotla_dH2) / std::pow(rotha + rotla, 2) + (std::pow(vibha, 2) * dvibla_dH2) / std::pow(vibha + vibla, 2));
}

void cooling_rate1(double &C1, double &C2, double &C3, double &C4, double &C5, double &C6, double &C7, double &C8, double &C9, double &C10, double &C11,
		double &C12, double &C13, double &dC1dT, double &dC2dT, double &dC3dT, double &dC4dT, double &dC5dT, double &dC6dT, double &dC7dT, double &dC8dT,
		double &dC9dT, double &dC10dT, double &dC11dT, double &dC12dT, double &dC13dT, double T, double z) {
	const auto T3 = T / 1e3;
	const auto T5 = T / 1e5;
	const auto T6 = T / 1e6;
	const auto tev = T / KperEv;
	const auto logtev = std::log(T);
	const auto tiny = std::numeric_limits<double>::min();

	const auto k1 = exp(
			-32.71396786375 + 13.53655609057 * logtev - 5.739328757388 * std::pow(logtev, 2) + 1.563154982022 * std::pow(logtev, 3)
					- 0.2877056004391 * std::pow(logtev, 4) + 0.03482559773736999 * std::pow(logtev, 5) - 0.00263197617559 * std::pow(logtev, 6)
					+ 0.0001119543953861 * std::pow(logtev, 7) - 2.039149852002e-6 * std::pow(logtev, 8));

	const auto dk1dT = (1 * 13.53655609057 - 2 * 5.739328757388 * logtev + 3 * 1.563154982022 * std::pow(logtev, 2) - 4 * 0.2877056004391 * std::pow(logtev, 3)
			+ 5 * 0.03482559773736999 * std::pow(logtev, 4) - 6 * 0.00263197617559 * std::pow(logtev, 5) + 7 * 0.0001119543953861 * std::pow(logtev, 6)
			- 8 * 2.039149852002e-6 * std::pow(logtev, 7)) * k1 / tev / KperEv;

	const auto k3 = exp(
			-44.09864886561001 + 23.91596563469 * logtev - 10.75323019821 * std::pow(logtev, 2) + 3.058038757198 * std::pow(logtev, 3)
					- 0.5685118909884001 * std::pow(logtev, 4) + 0.06795391233790001 * std::pow(logtev, 5) - 0.005009056101857001 * std::pow(logtev, 6)
					+ 0.0002067236157507 * std::pow(logtev, 7) - 3.649161410833e-6 * std::pow(logtev, 8));

	const auto dk3dT = (1 * 23.91596563469 - 2 * 10.75323019821 * logtev + 3 * 3.058038757198 * std::pow(logtev, 2)
			- 4 * 0.5685118909884001 * std::pow(logtev, 3) + 5 * 0.06795391233790001 * std::pow(logtev, 4) - 6 * 0.005009056101857001 * std::pow(logtev, 5)
			+ 7 * 0.0002067236157507 * std::pow(logtev, 6) - 8 * 3.649161410833e-6 * std::pow(logtev, 7)) * k3 / tev / KperEv;

	const auto k5 = exp(
			-68.71040990212001 + 43.93347632635 * logtev - 18.48066993568 * std::pow(logtev, 2) + 4.701626486759002 * std::pow(logtev, 3)
					- 0.7692466334492 * std::pow(logtev, 4) + 0.08113042097303 * std::pow(logtev, 5) - 0.005324020628287001 * std::pow(logtev, 6)
					+ 0.0001975705312221 * std::pow(logtev, 7) - 3.165581065665e-6 * std::pow(logtev, 8));

	const auto dk5dT = (1 * 43.93347632635 - 2 * 18.48066993568 * logtev + 3 * 4.701626486759002 * std::pow(logtev, 2)
			- 4 * 0.7692466334492 * std::pow(logtev, 3) + 5 * 0.08113042097303 * std::pow(logtev, 4) - 6 * 0.005324020628287001 * std::pow(logtev, 5)
			+ 7 * 0.0001975705312221 * std::pow(logtev, 6) - 8 * 3.165581065665e-6 * std::pow(logtev, 7)) * k5 / tev / KperEv;

	C1 = 7.50e-19 * (1.0 + std::sqrt(T / 100000)) * std::exp(-118348 / T);
	dC1dT = 7.50e-19 * ((118348. + (374.2492365256074 + 0.0015811388300841897 * T) * std::sqrt(T)) / (std::exp(118348 / T) * T * T));
	C2 = 9.10e-27 * (1.0 + std::sqrt(T / 100000)) * std::pow(T, -.1687) * std::exp(-13179 / T);
	dC2dT = 9.10e-27
			* ((13179 * std::pow(T, 1.8374) + 41.67565728335908 * std::pow(T, 2.3374) - 0.1687 * std::pow(T, 2.8374)
					+ 0.0010476625888137844 * std::pow(T, 3.3374)) / (std::exp(13179 / T) * std::pow(T, 4.0061)));
	C3 = 5.54e-17 * (1.0 + std::sqrt(T / 100000)) * std::pow(T, -.397) * std::exp(-473638 / T);
	dC3dT =
			5.54e-17
					* ((473638. * std::pow(T, 2.294) + 1497.774866406831 * std::pow(T, 2.794) - 0.397 * std::pow(T, 3.294)
							+ 0.00032571459899734327 * std::pow(T, 3.794)) / (std::exp(473638 / T) * std::pow(T, 4.691)));
	C4 = 2.18e-11 * k1;
	dC4dT = 2.18e-11 * dk1dT;
	C5 = 3.94e-11 * k3;
	dC5dT = 3.94e-11 * dk3dT;
	C6 = 8.72e-11 * k5;
	dC6dT = 8.72e-11 * dk5dT;
	C7 = 5.01e-27 * (1.0 + std::sqrt(T / 100000)) * std::pow(T, -.1687) * std::exp(-55338 / T);
	dC7dT = 5.01e-27
			* ((55338 * std::pow(T, 1.8374) + 174.9941211583978 * std::pow(T, 2.3374) - 0.1687 * std::pow(T, 2.8374)
					+ 0.0010476625888137844 * std::pow(T, 3.3374)) / (std::exp(55338 / T) * std::pow(T, 4.0061)));
	C8 = 8.70e-27 * std::sqrt(T) * std::pow(T3, -0.2) / (1.0 + std::pow(T6, 0.7));
	dC8dT = 8.70e-27 * (-((25238.293779207717 - 3e8 / std::pow(T, 0.7)) / (2.5118864315095773e8 + 31697.863849222253 * std::pow(T, 0.7) + std::pow(T, 1.4))));
	C9 = 1.55e-26 * std::pow(T, 0.3647);
	dC9dT = 0.3647 * C9 / T;
	C10 = 1.24e-13 * std::pow(T, -1.5) * (1 + 0.3 * std::exp(-94000 / T)) * std::exp(-470000 / T);
	dC10dT = 1.24e-13
			* (((169200 + 470000 * std::exp(94000 / T)) * std::pow(T, 2.5) + (-0.45 - 1.5 * std::exp(94000 / T)) * std::pow(T, 3.5))
					/ (std::exp(564000 / T) * std::pow(T, 6.)));
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
	chemical_rates(K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, T);
	radiation_rates(J20, J21, J22, J23, J24, J25, J26, J27, J28, z);
//	J20 = J21 = J22 = 0.0;
//	K3 = K4 = K5 = K6 = 0.0;
	double err;
#define List(a0,a1,a2,a3,a4,a5) {a0,a1,a2,a3,a4,a5}
#define Power std::pow
	do {
		const std::array<double, NF> f = { -H + H0 - Hp + Hp0, Hp - Hp0 + dt * Hp * K2 * ne - dt * H * (J20 + K1 * ne), -He + He0 - Hep + Hep0 - Hepp + Hepp0,
				Hep - Hep0 + dt * Hep * (J22 + (K4 + K5) * ne) - dt * (He * J21 + He * K3 * ne + Hepp * K6 * ne), Hepp - Hepp0 + dt * Hepp * K6 * ne
						- dt * Hep * (J22 + K5 * ne), -Hep - 2 * Hepp - Hp + ne };
		const std::array<std::array<double, NF>, NF> A = { { List(-1,-1,0,0,0,0), List(-(dt*(J20 + K1*ne)),1 + dt*K2*ne,0,0,0,dt*(-(H*K1) + Hp*K2)),
				List(0, 0, -1, -1, -1, 0),
		List(0,0,-(dt*(J21 + K3*ne)),1 + dt*(J22 + (K4 + K5)*ne),-(dt*K6*ne),-(dt*(He*K3 - Hep*(K4 + K5) + Hepp*K6))),
				List(0, 0, 0, -(dt * (J22 + K5 * ne)), 1 + dt * K6 * ne, dt * (-(Hep * K5) + Hepp * K6)), List(0,-1,0,-1,-2,1) } };
		const auto Ainv = matrix_inverse<NF>(A);
		for (int n = 0; n < NF; n++) {
			auto dU = 0.0;
			for (int m = 0; m < NF; m++) {
				dU += Ainv[n][m] * f[m];
			}
			U[n] = U[n] - dU;
			if( U[n] < 0.0 ) {
				return false;
			}
		}
		err = 0.0;
		for (int i = 0; i < NF; i++) {
			err += f[i] * f[i];
		}
		err = std::sqrt(err);
		err /= (H + Hp + He + Hep + Hepp + ne);
	} while (err > 1.0e-7);
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
	s0.H = 0.0 * n;
	s0.Hp = 0.92 * n;
	s0.He = 0.00 * n;
	s0.Hep = 0.00 * n;
	s0.Hepp = 0.08 * n;
	double T0 = 1e+3;
	auto s = s0;
	printf("%14s %14s %14s %14s %14s %14s %14s %14s %14s %14s\n", "t", "T0", "T", "na", "ne", "H", "Hp", "He", "Hep", "Hepp");
	for (double dt = 1e1; dt < 1e18; dt *= 10.0) {
		auto T = chemistry_update(s0, s, species_energy(s0, T0), 1.0, 1.0, dt);
		const auto ne = s.Hp + s.Hep + 2 * s.Hepp;
		const auto na = s.Hp + s.H + 4 * s.He + 4 * s.Hep + 4 * s.Hepp;
		printf("%14.4e %14.4e %14.4e %14.4e %14.4e %14.4e %14.4e %14.4e %14.4e %14.4e\n", dt / (3600 * 24 * 365), T0, T, na, ne, s.H, s.Hp, s.He, s.Hep,
				s.Hepp);
	}

}
