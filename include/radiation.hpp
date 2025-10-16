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

template<typename Type>
struct ImplicitEnergyFailure: public std::runtime_error {
	explicit ImplicitEnergyFailure(std::string const &msg, Type rhoVal, Type epsVal, Type EVal, Type fVal, Type dfVal, Type stepVal, int iter) :
			std::runtime_error(makeMessage(msg, rhoVal, epsVal, EVal, fVal, dfVal, stepVal, iter)) {
	}
private:
	static std::string makeMessage(std::string const &msg, Type rhoVal, Type epsVal, Type EVal, Type fVal, Type dfVal, Type stepVal, int iter) {
		std::ostringstream out;
		out << msg << "\n" << "  ρ       = " << rhoVal << "\n" << "  ε       = " << epsVal << "\n" << "  E_r     = " << EVal << "\n" << "  f       = " << fVal
				<< "\n" << "  df/de   = " << dfVal << "\n" << "  Δe      = " << stepVal << "\n" << "  iter    = " << iter;
		return out.str();
	}
};

template<typename Type, int dimensionCount>
struct RadConserved {
	EnergyDensityType<Type> E;
	Vector<EnergyFluxType<Type>, dimensionCount> F;
	auto temperature() const;
	auto pressure() const;
	auto stressEnergy() const;
};

template<typename Type>
inline EnergyDensityType<Type> temperature2radiationEnergy(TemperatureType<Type> const &T) {
	using namespace Constants;
	return aR * sqr(sqr(T));
}

template<typename Type, int dimensionCount>
auto RadConserved<Type, dimensionCount>::stressEnergy() const {
	using namespace Constants;
	return space2spaceTime<EnergyDensityType<Type>, dimensionCount>(E, ic * F, pressure());
}

template<typename Type, int dimensionCount>
auto RadConserved<Type, dimensionCount>::pressure() const {
	using namespace Constants;
	using std::max;
	auto const magF = sqrt(F.dot(F));
	auto const imagF = 1 / max(magF, EnergyFluxType<Type>(tiny));
	auto const iE = 1 / E;
	auto const f = magF * ic * iE;
	auto const n = F * imagF;
	auto const f2 = sqr(f);
	auto const ξ = (3 + 4 * f2) / (5 + 2 * sqrt(4 - 3 * f2));
	auto const Ωdif = 0.5 * (1 - ξ);
	auto const Ωstr = 0.5 * (3 * ξ - 1);
	auto const P = E * (Ωdif * δ<Type, dimensionCount> + Ωstr * sqr(n));
	return P;
}

template<typename Type, int dimensionCount>
auto RadConserved<Type, dimensionCount>::temperature() const {
	using namespace Constants;
	return pow<Rational(1, 4)>(E / aR);
}

template<typename Type, int dimensionCount>
void implicitEnergySolve(SpecificEnergyType<Type> &ε, EnergyDensityType<Type> &E, MassDensityType<Type> const &ρ_, SpecificAreaType<Type> const &κₐ,
		EquationOfState<Type> const &eos, TimeType<Type> const &dt) {
	using namespace Constants;
	using std::abs;
	using std::copysign;
	using std::max;
	using std::min;

	constexpr Type toler = 16 * eps;
	constexpr int maxIter = 100;
	using Dimensionless = FwdAutoDiff<DimensionlessType<Type>, 2, std::tuple<DimensionlessType<Type>>>;
	using MassDensity = FwdAutoDiff<MassDensityType<Type> , 2, std::tuple<DimensionlessType<Type>>>;
	using EnergyDensity = FwdAutoDiff<EnergyDensityType<Type> , 2, std::tuple<DimensionlessType<Type>>>;
	using SpecificEnergy = FwdAutoDiff<SpecificEnergyType<Type> , 2, std::tuple<DimensionlessType<Type>>>;
	constexpr auto one = Dimensionless(DimensionlessType(1.0));
	auto const ρ = MassDensity(MassDensityType<Type>(ρ_));
	auto const Eg0 = EnergyDensity(ρ * ε);
	auto const Er0 = EnergyDensity(E);
	auto const Etot0 = EnergyDensity(Er0 + Eg0);
	auto const λ = Dimensionless(ρ * c * κₐ * dt);
	auto const Fgas = [eos, Etot0, Er0, Eg0, λ, ρ, one](SpecificEnergyType<Type> ε_) {
		auto const ε = SpecificEnergy::template independent<0>(ε_);
		auto const T = eos.energy2temperature(ρ, ε);
		auto const T2 = sqr(T);
		auto const T4 = sqr(T2);
		auto const B = aR * T4;
		auto const E = Etot0 - ρ * ε;
		auto const f = ((one + λ) * E - (Er0 + λ * B)) / (Er0 + Eg0);
		return f;
	};
	auto const Frad = [eos, Etot0, Er0, Eg0, λ, ρ, one](EnergyDensityType<Type> E_) {
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
	Type error = 1.0;
	int iter = 0;
	for (; error > toler; ++iter) {
		if (iter > maxIter) {
			throw ImplicitEnergyFailure("implicitEnergySolve failed to converge", ρ_.value(), ε.value(), E.value(), std::numeric_limits < Type > ::quiet_NaN(),
					std::numeric_limits < Type > ::quiet_NaN(), std::numeric_limits < Type > ::quiet_NaN(), iter);
		}
		bool gasForm = (ρ_ * ε) < E;
		Dimensionless fRoot = gasForm ? Fgas(ε) : Frad(E);
		auto const f = fRoot.template get<0>().value();
		auto const dfde = fRoot.template get<1>().value();
		auto const d2fde2 = fRoot.template get<2>().value();
		if (!std::isfinite(f) || !std::isfinite(dfde) || !std::isfinite(d2fde2)) {
			throw ImplicitEnergyFailure("NaN or Inf detected in implicitEnergySolve derivatives", ρ_.value(), ε.value(), E.value(), f, dfde, 0.0, iter);
		}
		auto const den = Type(max(sqr(dfde), 2 * sqr(dfde) - f * d2fde2));
		auto const num = -2 * f * dfde;
		auto const de = num / (den + copysign(tiny, den));
		if (!std::isfinite(de)) {
			throw ImplicitEnergyFailure("Non-finite update step in implicitEnergySolve", ρ_.value(), ε.value(), E.value(), f, dfde, de, iter);
		}
		auto ε1 = ε;
		auto E1 = E;
		if (gasForm) {
			ε1 += SpecificEnergyType<Type>(de);
			E1 = (Etot0 - ρ * ε1).template get<0>();
		} else {
			E1 += EnergyDensityType<Type>(de);
			ε1 = ((Etot0 - E1) / ρ_).template get<0>();
		}
		ε = max(0.5 * ε, ε1);
		E = max(0.5 * E, E1);
		Type const e = gasForm ? ε.value() : E.value();
		error = abs((de / e));
		if (!std::isfinite(error) || error > 1e10) {
			throw ImplicitEnergyFailure("Divergence detected in implicitEnergySolve", ρ_.value(), ε.value(), E.value(), f, dfde, de, iter);
		}
	}
}

template<typename Type, int dimensionCount>
void implicitRadiationSolve(GasConserved<Type, dimensionCount> &gasCon, RadConserved<Type, dimensionCount> &radCon, Opacity<Type> const &opac,
		EquationOfState<Type> const &eos, TimeType<Type> const &dt) {
	using std::max;
	using std::min;
	using namespace Constants;
	constexpr Type toler = 4 * eps;
	enableFPE();
	GasPrimitive<Type, dimensionCount> gasPrim;
	EnergyDensityType<Type> &τ = gasCon.τ;
	EnergyDensityType<Type> &E = radCon.E;
	EntropyDensityType<Type> &K = gasCon.K;
	MassDensityType<Type> &D = gasCon.D;
	MassDensityType<Type> &ρ = gasPrim.ρ;
	SpecificEnergyType<Type> &ε = gasPrim.ε;
	Vector<DimensionlessType<Type>, dimensionCount> &β = gasPrim.β;
	Vector<EnergyFluxType<Type>, dimensionCount> &F = radCon.F;
	Vector<MomentumDensityType<Type>, dimensionCount> &S = gasCon.S;
	auto const norm = sqrt(sqr(sqrt(F.dot(F))) * ic2 + sqr(E) + c2 * sqr(sqrt(S.dot(S))) + sqr(τ)).value();
	Type error = 1.0;
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
		auto const dt0 = dt / γ;
		auto const Λ = space2spaceTime<DimensionlessType<Type>>(γ, γ * β, δ<Type, dimensionCount> + (sqr(γ) / (γ + 1)) * sqr(β));
		auto const iΛ = space2spaceTime<DimensionlessType<Type>>(γ, -γ * β, δ<Type, dimensionCount> + (sqr(γ) / (γ + 1)) * sqr(β));
		auto const λ = ρ * c * (opac.κₐ + opac.κₛ) * dt0;
		auto R = symmetric(Λ * R0 * Λ);
		E = spaceTime2Tensor0(R);
		F = c * spaceTime2Tensor1(R);
		implicitEnergySolve<Type, dimensionCount>(ε, E, ρ, opac.κₛ, eos, dt0);
		F = F / (1 + λ);
		R = radCon.stressEnergy();
		R = symmetric(iΛ * R * iΛ);
		E = spaceTime2Tensor0(R);
		F = c * spaceTime2Tensor1(R);
		S = S0 + (F - F0) * ic2;
		auto const τo = c2 * (γ - 1.0) * D;
		auto const τmin = 0.5 * (τo + τ);
		τ = τ0 + E0 - E;
		τ = max(τmin, τ);
		auto const T = eos.energy2temperature(ρ, ε);
		K = D * eos.temperature2entropy(ρ, T);
		error = Type(sqrt(sqr(sqrt((F - F1).dot(F - F1))) * ic2 + sqr(E - E1)).value() / norm);
	}
}

#undef XXX

