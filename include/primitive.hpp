/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/
#ifndef INCLUDE_SRHD_PRIMITIVE_HPP_
#define INCLUDE_SRHD_PRIMITIVE_HPP_

#include "srhd.hpp"

template<typename Type, int dimensionCount>
struct GasPrimitive {
	MassDensityType<Type> ρ;
	SpecificEnergyType<Type> ε;
	Vector<DimensionlessType<Type>, dimensionCount> β;
	auto dualEnergySwitch(EquationOfState<Type> const &eos) const {
		using namespace Constants;
		auto const h = enthalpy(eos);
		auto const γ = lorentzFactor();
		auto const num = ρ * ε;
		auto const den = ρ * γ * (γ * h - c2);
		return num / den;
	}
	SpecificEnergyType<Type> enthalpy(EquationOfState<Type> const &eos) const {
		using namespace Constants;
		return c2 + ε + eos.energy2pressure(ρ, ε) / ρ;
	}
	DimensionlessType<Type> lorentzFactor() const {
		return sqrt(1 / (1 - β.dot(β)));
	}
	GasConserved<Type, dimensionCount> toConserved(EquationOfState<Type> const &eos) const {
		using namespace Constants;
		GasConserved<Type, dimensionCount> con;
		MassDensityType<Type> &D = con.D;
		EnergyDensityType<Type> &τ = con.τ;
		EntropyDensityType<Type> &K = con.K;
		Vector<MomentumDensityType<Type>, dimensionCount> &S = con.S;
		auto const γ = lorentzFactor();
		auto const γ2 = sqr(γ);
		auto const iρ = 1 / ρ;
		auto const p = eos.energy2pressure(ρ, ε);
		auto const T = eos.energy2temperature(ρ, ε);
		auto const ϰ = eos.temperature2entropy(ρ, T);
		auto const h = c2 + ε + p * iρ;
		S = ρ * γ2 * h * β * ic;
		D = ρ * γ;
		τ = ρ * γ2 * ε + (γ2 - 1) * p + c2 * D * (γ - 1);
		K = D * ϰ;
		return con;
	}
	auto temperature(EquationOfState<Type> const &eos) const {
		return eos.energy2temperature(ρ, ε);
	}
};

#endif /* INCLUDE_SRHD_PRIMITIVE_HPP_ */
