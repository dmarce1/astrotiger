/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/
#ifndef INCLUDE_SRHD_PRIMITIVE_HPP_
#define INCLUDE_SRHD_PRIMITIVE_HPP_

#include "autodiff.hpp"
#include "tensor.hpp"

template<typename Type, int dimensionCount>
struct GasPrimitive {
	auto dualEnergySwitch(EquationOfState<Type> const &eos) const {
		auto const h = eos.energy2enthalpy(ρ, ε);
		auto const γ = lorentzFactor();
		auto const num = ρ * ε;
		auto const den = ρ * γ * (γ * h - sqr(pc.c));
		return num / den;
	}
	DimensionlessType<Type> lorentzFactor() const {
		return sqrt(1 / (1 - v.dot(v) / sqr(pc.c)));
	}
	GasConserved<Type, dimensionCount> toConserved(EquationOfState<Type> const &eos) const {
		GasConserved<Type, dimensionCount> con;
		MassDensityType<Type> &D = con.D;
		EnergyDensityType<Type> &τ = con.τ;
		EntropyDensityType<Type> &K = con.K;
		Vector<MomentumDensityType<Type>, dimensionCount> &S = con.S;
		auto const γ = lorentzFactor();
		auto const γ2 = sqr(γ);
		auto const p = eos.energy2pressure(ρ, ε);
		auto const T = eos.energy2temperature(ρ, ε);
		auto const ϰ = eos.temperature2entropy(ρ, T);
		auto const h = eos.energy2enthalpy(ρ, ε);
		S = ρ * γ2 * h * v / sqr(pc.c);
		D = ρ * γ;
		τ = ρ * γ2 * ε + (γ2 - 1) * p + sqr(pc.c) * D * (γ - 1);
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
//	auto eigenvalues(EquationOfState<Type> const &eos, int dim) const {
//		constexpr int fieldCount = 2 + dimensionCount;
//		constexpr auto ic2 = 1 / sqr(pc.c);
//		Vector<VelocityType<Type>, fieldCount> λ;
//		auto const u = v[dim];
//		auto const a = eos.energy2soundSpeed(ρ, ε);
//		λ.front() = (u - a) / (1 - u * a * ic2);
//		for (int k = 1; k + 1 < fieldCount; k++) {
//			λ[k] = u;
//		}
//		λ.back() = (u + a) / (1 + u * a * ic2);
//		return λ;
//	}
	auto flux(EquationOfState<Type> const &eos, int direction) const {
		constexpr auto δ = identity<DimensionlessType<Type>, dimensionCount>();
		constexpr auto c = pc.c;
		constexpr auto c2 = c * c;
		GasFlux<Type, dimensionCount> flux;
		auto const γ = lorentzFactor();
		auto const γ2 = sqr(γ);
		auto const u = v[direction];
		auto const p = eos.energy2pressure(ρ, ε);
		auto const h = c2 + ε + p / ρ;
		auto const W = ρ * γ2 * h;
		flux.D = γ * ρ * u;
		flux.S = W * v * u;
		flux.τ = (W - p - ρ * γ) * u;
		flux.S[direction] += p;
		return flux;
	}
	friend GasConserved<Type, dimensionCount> ;
	friend RadConserved<Type, dimensionCount> ;
	template<typename T, int D>
	friend auto hllc(GasConserved<T, D> const&, GasConserved<T, D> const&, EquationOfState<T> const&, int);
private:
	static constexpr PhysicalConstants<Type> pc { };
	MassDensityType<Type> ρ;
	SpecificEnergyType<Type> ε;
	Vector<VelocityType<Type>, dimensionCount> v;
};

#endif /* INCLUDE_SRHD_PRIMITIVE_HPP_ */
