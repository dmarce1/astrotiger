/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once
#include "constants.hpp"
#include "opacity.hpp"
#include "srhd.hpp"
#include "conserved.hpp"
#include "primitive.hpp"
#include "eos.hpp"
#include "fpe.hpp"

#include <iostream>
#include <climits>
#include <functional>
#include <tuple>
#include "autodiff.hpp"



struct ImplicitEnergyFailure: public std::runtime_error {
	explicit ImplicitEnergyFailure(std::string const &msg, double rhoVal, double epsVal, double EVal, double fVal, double dfVal, double stepVal, int iter) :
			std::runtime_error(makeMessage(msg, rhoVal, epsVal, EVal, fVal, dfVal, stepVal, iter)) {
	}
private:
	static std::string makeMessage(std::string const &msg, double rhoVal, double epsVal, double EVal, double fVal, double dfVal, double stepVal, int iter) {
		std::ostringstream out;
		out << msg << "\n" << "  ρ       = " << rhoVal << "\n" << "  ε       = " << epsVal << "\n" << "  E_r     = " << EVal << "\n" << "  f       = " << fVal
				<< "\n" << "  df/de   = " << dfVal << "\n" << "  Δe      = " << stepVal << "\n" << "  iter    = " << iter;
		return out.str();
	}
};

template<int NDIM>
struct RadConserved {
	EnergyDensityType E;
	Vector<EnergyFluxType, NDIM> F;
	auto temperature() const;
	auto pressure() const;
	auto stressEnergy() const;
};

inline EnergyDensityType temperature2radiationEnergy(TemperatureType const &T) {
	using namespace Constants;
	return aR * sqr(sqr(T));
}

template<int NDIM>
auto RadConserved<NDIM>::stressEnergy() const {
	using namespace Constants;
	return space2spaceTime<EnergyDensityType, NDIM>(E, ic * F, pressure());
}

template<int NDIM>
auto RadConserved<NDIM>::pressure() const {
	using namespace Constants;
	using std::max;
	auto const magF = sqrt(F.dot(F));
	auto const imagF = 1 / max(magF, EnergyFluxType(tiny));
	auto const iE = 1 / E;
	auto const f = magF * ic * iE;
	auto const n = F * imagF;
	auto const f2 = sqr(f);
	auto const ξ = (3 + 4 * f2) / (5 + 2 * sqrt(4 - 3 * f2));
	auto const Ωdif = 0.5 * (1 - ξ);
	auto const Ωstr = 0.5 * (3 * ξ - 1);
	auto const P = E * (Ωdif * δ < NDIM > +Ωstr * (n * n));
	return P;
}

template<int NDIM>
auto RadConserved<NDIM>::temperature() const {
	using namespace Constants;
	return pow<Rational(1, 4)>(E / aR);
}

template<int NDIM>
void implicitEnergySolve(SpecificEnergyType &ε, EnergyDensityType &E, MassDensityType const &ρ_, SpecificAreaType const &κₐ, EquationOfState const &eos,
		TimeType const &dt) {
	using namespace Constants;
	using std::abs;
	using std::copysign;
	using std::max;
	using std::min;

	constexpr double toler = 16 * eps;
	constexpr int maxIter = 100;
	using Dimensionless = FwdAutoDiff<DimensionlessType, 2, std::tuple<DimensionlessType>>;
	using MassDensity = FwdAutoDiff<MassDensityType, 2, std::tuple<DimensionlessType>>;
	using EnergyDensity = FwdAutoDiff<EnergyDensityType, 2, std::tuple<DimensionlessType>>;
	using SpecificEnergy = FwdAutoDiff<SpecificEnergyType, 2, std::tuple<DimensionlessType>>;
	constexpr auto one = Dimensionless(DimensionlessType(1.0));
	auto const ρ = MassDensity(MassDensityType(ρ_));
	auto const Eg0 = EnergyDensity(ρ * ε);
	auto const Er0 = EnergyDensity(E);
	auto const Etot0 = EnergyDensity(Er0 + Eg0);
	auto const λ = Dimensionless(ρ * c * κₐ * dt);
	auto const Fgas = [eos, Etot0, Er0, Eg0, λ, ρ, one](SpecificEnergyType ε_) {
		auto const ε = SpecificEnergy::template independent<0>(ε_);
		auto const T = eos.energy2temperature(ρ, ε);
		auto const T2 = sqr(T);
		auto const T4 = sqr(T2);
		auto const B = aR * T4;
		auto const E = Etot0 - ρ * ε;
		auto const f = ((one + λ) * E - (Er0 + λ * B)) / (Er0 + Eg0);
		return f;
	};
	auto const Frad = [eos, Etot0, Er0, Eg0, λ, ρ, one](EnergyDensityType E_) {
		auto const E = EnergyDensity::template independent<0>(E_);
		auto const ρε = Etot0 - E;
		auto const ε = ρε / ρ;
		auto const T = eos.energy2temperature(ρ, ε);
		auto const T2 = sqr(T);
		auto const T4 = sqr(T2);
		auto const B = aR * T4;
		auto const f = (λ * B - ((one + λ) * E - Er0)) / (Er0 + Eg0);
		return f;
	};
	double error = 1.0;
	int iter = 0;
	for (; error > toler; ++iter) {
		if (iter > maxIter) {
			throw ImplicitEnergyFailure("implicitEnergySolve failed to converge", ρ_.value(), ε.value(), E.value(), std::numeric_limits<double>::quiet_NaN(),
					std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), iter);
		}
		bool gasForm = (ρ_ * ε) < E;
		Dimensionless fRoot = gasForm ? Fgas(ε) : Frad(E);
		auto const f = fRoot.template get<0>().value();
		auto const dfde = fRoot.template get<1>().value();
		auto const d2fde2 = fRoot.template get<2>().value();
		if (!std::isfinite(f) || !std::isfinite(dfde) || !std::isfinite(d2fde2)) {
			throw ImplicitEnergyFailure("NaN or Inf detected in implicitEnergySolve derivatives", ρ_.value(), ε.value(), E.value(), f, dfde, 0.0, iter);
		}
		auto const den = double(max(sqr(dfde), 2 * sqr(dfde) - f * d2fde2));
		auto const num = -2 * f * dfde;
		auto const de = num / (den + copysign(tiny, den));
		if (!std::isfinite(de)) {
			throw ImplicitEnergyFailure("Non-finite update step in implicitEnergySolve", ρ_.value(), ε.value(), E.value(), f, dfde, de, iter);
		}
		auto ε1 = ε;
		auto E1 = E;
		if (gasForm) {
			ε1 += SpecificEnergyType(de);
			E1 = (Etot0 - ρ * ε1).template get<0>();
		} else {
			E1 += EnergyDensityType(de);
			ε1 = ((Etot0 - E1) / ρ_).template get<0>();
		}
		ε = max(0.5 * ε, ε1);
		E = max(0.5 * E, E1);
		double const e = gasForm ? ε.value() : E.value();
		error = abs((de / e));
		if (!std::isfinite(error) || error > 1e10) {
			throw ImplicitEnergyFailure("Divergence detected in implicitEnergySolve", ρ_.value(), ε.value(), E.value(), f, dfde, de, iter);
		}
	}
}

template<int NDIM>
void implicitRadiationSolve(GasConserved<NDIM> &gasCon, RadConserved<NDIM> &radCon, Opacity const &opac, EquationOfState const &eos, TimeType const &dt) {
	using std::max;
	using std::min;
	using namespace Constants;
	constexpr double toler = 4 * eps;
	enableFPE();
	GasPrimitive<NDIM> gasPrim;
	EnergyDensityType &τ = gasCon.τ;
	EnergyDensityType& E = radCon.E;
	EntropyDensityType &K = gasCon.K;
	MassDensityType &D = gasCon.D;
	MassDensityType &ρ = gasPrim.ρ;
	SpecificEnergyType &ε = gasPrim.ε;
	Vector<DimensionlessType, NDIM> &β = gasPrim.β;
	Vector<EnergyFluxType, NDIM>& F = radCon.F;
	Vector<MomentumDensityType, NDIM> &S = gasCon.S;
	auto const norm = sqrt(sqr(sqrt(F.dot(F))) * ic2 + sqr(E) + c2 * sqr(sqrt(S.dot(S))) + sqr(τ)).value();
	double error = 1.0;
	auto const F0 = F;
	auto const E0 = E;
	auto const S0 = S;
	auto const τ0 = τ;
	auto R0 = radCon.stressEnergy();
	for (int j = 0; error > toler; j++) {
		auto const F1 = F;
		auto const E1 = E;
		gasPrim = gasCon.toPrimitive(eos);
		auto const γ = gasPrim.lorentzFactor();
		auto const Λ = space2spaceTime<DimensionlessType>(γ, γ * β, δ < NDIM > + (sqr(γ) / (γ + 1)) * (β * β));
		auto const dt0 = dt / γ;
		auto const iΛ = space2spaceTime<DimensionlessType>(γ, -γ * β, δ < NDIM > +(sqr(γ) / (γ + 1)) * (β * β));
		auto const λ = ρ * c * (opac.κₐ + opac.κₛ) * dt0;
		auto R = Λ * R0 * Λ;
		E = spaceTime2Tensor0(R);
		F = c * spaceTime2Tensor1(R);
		implicitEnergySolve<NDIM>(ε, E, ρ, opac.κₛ, eos, dt0);
		F = F / (1 + λ);
		R = radCon.stressEnergy();
		R = iΛ * R * iΛ;
		E = spaceTime2Tensor0(R);
		F = c * spaceTime2Tensor1(R);
		S = S0 + (F - F0) * ic2;
		auto const τo = c2 * (γ - 1.0) * D;
		auto const τmin = 0.5 * (τo + τ);
		τ = τ0 + E0 - E;
		τ = max(τmin, τ);
		auto const T = eos.energy2temperature(ρ, ε);
		K = D * eos.temperature2entropy(ρ, T);
		error = double(sqrt(sqr(sqrt((F - F1).dot(F - F1))) * ic2 + sqr(E - E1)).value() / norm);
	}
}

#undef XXX

