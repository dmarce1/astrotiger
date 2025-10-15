/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/
#ifndef INCLUDE_SRHD_PRIMITIVE_HPP_
#define INCLUDE_SRHD_PRIMITIVE_HPP_

#include "srhd.hpp"

template<int NDIM>
struct GasPrimitive {
	MassDensityType ρ;
	SpecificEnergyType ε;
	Vector<DimensionlessType, NDIM> β;
	auto dualEnergySwitch(EquationOfState const &eos) const {
		using namespace Constants;
		auto const h = enthalpy(eos);
		auto const γ = lorentzFactor();
		auto const num = ρ * ε;
		auto const den = ρ * γ * (γ * h - c2);
		return num / den;
	}
	SpecificEnergyType enthalpy(EquationOfState const &eos) const {
		using namespace Constants;
		return c2 + ε + eos.energy2pressure(ρ, ε) / ρ;
	}
	DimensionlessType lorentzFactor() const {
		return sqrt(1 / (1 - β.dot(β)));
	}
	GasConserved<NDIM> toConserved(EquationOfState const &eos) const {
		using namespace Constants;
		GasConserved<NDIM> con;
		MassDensityType &D = con.D;
		EnergyDensityType &τ = con.τ;
		EntropyDensityType &K = con.K;
		Vector<MomentumDensityType, NDIM> &S = con.S;
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
	auto temperature(EquationOfState const &eos) const {
		return eos.energy2temperature(ρ, ε);
	}
};

#endif /* INCLUDE_SRHD_PRIMITIVE_HPP_ */
