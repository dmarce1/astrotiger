/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/
#ifndef INCLUDE_SRHD_PRIMITIVE_HPP_
#define INCLUDE_SRHD_PRIMITIVE_HPP_

#include "tensor.hpp"

template<typename Type, int dimensionCount>
struct GasPrimitive {
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
		using namespace Constants;
		return sqrt(1 / (1 - v.dot(v) * ic2));
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
		S = ρ * γ2 * h * v * ic2;
		D = ρ * γ;
		τ = ρ * γ2 * ε + (γ2 - 1) * p + c2 * D * (γ - 1);
		K = D * ϰ;
		return con;
	}
	auto temperature(EquationOfState<Type> const &eos) const {
		return eos.energy2temperature(ρ, ε);
	}
	void setMassDensity(MassDensityType<Type> const &ρ_) {
		ρ = ρ_;
	}
	void setSpecificEnergy(SpecificEnergyType<Type> const &ε_) {
		ε = ε_;
	}
	void setVelocity(Vector<VelocityType<Type>, dimensionCount> const &v_) {
		v = v_;
	}
	void setVelocity(int k, VelocityType<Type> const &v_) {
		v[k] = v_;
	}
	MassDensityType<Type> getMassDensity() const {
		return ρ;
	}
	SpecificEnergyType<Type> getSpecificEnergy() const {
		return ε;
	}
	Vector<VelocityType<Type>, dimensionCount> getVelocity() const {
		return v;
	}
	VelocityType<Type> getVelocity(int k) const {
		return v[k];
	}
	friend GasConserved<Type, dimensionCount> ;
	friend RadConserved<Type, dimensionCount> ;
private:
	MassDensityType<Type> ρ;
	SpecificEnergyType<Type> ε;
	Vector<VelocityType<Type>, dimensionCount> v;
};

#endif /* INCLUDE_SRHD_PRIMITIVE_HPP_ */
