/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/
#ifndef INCLUDE_SRHD_PRIMITIVE_HPP_
#define INCLUDE_SRHD_PRIMITIVE_HPP_

#include <cassert>
#include <iostream>

#include "fpe.hpp"
#include "autodiff.hpp"
#include "eos.hpp"
#include "gas_conserved.hpp"
#include "tensor.hpp"
#include "velocity.hpp"

template<typename Type, int dimensionCount>
struct GasPrimitive {
	auto dualEnergySwitch(EquationOfState<Type> const &eos) const {
		auto const h = eos.energy2enthalpy(ρ, ε);
		auto const γ = fourVelocity2LorentzFactor(u);
		return ε / (γ * (γ * h - c2));
	}
	GasConserved<Type, dimensionCount> toConserved(EquationOfState<Type> const &eos) const {
		GasConserved<Type, dimensionCount> con;
		MassDensityType<Type> &D = con.D;
		EnergyDensityType<Type> &τ = con.τ;
		EntropyDensityType<Type> &K = con.K;
		Vector<MomentumDensityType<Type>, dimensionCount> &S = con.S;
		auto const γ = fourVelocity2LorentzFactor(u);
		auto const v = fourVelocity2CoordVelocity(u);
		auto const γ2 = sqr(γ);
		auto const p = eos.energy2pressure(ρ, ε);
		auto const T = eos.energy2temperature(ρ, ε);
		auto const ϰ = eos.temperature2entropy(ρ, T);
		auto const h = c2 + ε + p / ρ;
		auto const W = ρ * γ2 * h;
		S = W * v * ic2;
		D = ρ * γ;
		τ = ρ * γ2 * ε + (γ2 - 1) * p + c2 * D * (γ - 1);
		K = D * ϰ;
		return con;
	}
	auto temperature(EquationOfState<Type> const &eos) const {
		return eos.energy2temperature(ρ, ε);
	}
	auto flux(EquationOfState<Type> const &eos, int direction) const {
		constexpr auto δ = identity<DimensionlessType<Type>, dimensionCount>();
		constexpr auto c2 = c * c;
		GasFlux<Type, dimensionCount> flux;
		auto const γ = fourVelocity2LorentzFactor(u);
		auto const v = fourVelocity2CoordVelocity(u);
		auto const γ2 = sqr(γ);
		auto const vˣ = v[direction];
		auto const T = eos.energy2temperture(ρ, ε);
		auto const κ = eos.temperture2entropy(ρ, T);
		auto const p = eos.energy2pressure(ρ, ε);
		auto const h = c2 + ε + p / ρ;
		auto const W = ρ * γ2 * h;
		flux.D = γ * ρ * vˣ;
		flux.K = γ * κ * vˣ;
		flux.S = W * v * vˣ;
		flux.τ = (W - p - c2 * ρ * γ) * vˣ;
		flux.S[direction] += p;
		return flux;
	}
	auto jacobian(EquationOfState<Type> const &eos, int ni) const {
		using std::abs;
		enableFPE();
		using Auto1 = FwdAutoDiff<Type, 1, std::array<Type, fieldCount>>;
		constexpr Type c = PhysicalConstants<Type>::c.value();
		constexpr Type c2 = sqr(c);
		constexpr Auto1 half(0.5);
		constexpr Auto1 one(1);
		constexpr int ti = dimensionCount;
		Vector<Auto1, dimensionCount> ρcu;
		Auto1 const ρc2 = Auto1::independent(Type(ρ * c2), Di);
		for (int d = 0; d < dimensionCount; d++) {
			ρcu[d] = Auto1::independent(Type(ρ * c * u[d]), 1 + d);
		}
		Auto1 const ρε = Auto1::independent(Type(ρ * ε), τi);
		auto const u0 = ρcu / ρc2;
		Vector<Auto1, dimensionCount + 1> Uᵛ;
		for (int k = 0; k < dimensionCount; k++) {
			Uᵛ[k] = u0[k];
		}
		Uᵛ[ti] = sqrt(one + u0.dot(u0));
		Tensor2<Auto1, dimensionCount + 1> const UᵛUᵛ = symmetrize(Uᵛ * Uᵛ);
		Auto1 const p = eos.energy2pressure(ρc2 / c2, ρε * c2 / ρc2);
		Auto1 const ρh = ρc2 + ρε + p;
		auto const T = ρh * UᵛUᵛ + η * p;
		SquareMatrix<Type, fieldCount> dUdV;
		SquareMatrix<Type, fieldCount> dFdV;
		for (int m = 0; m < fieldCount; m++) {
			auto const μ = Auto1::index_type::unit(m);
			dUdV(0, m) = (ρc2 * Uᵛ[ti])[μ];
			dFdV(0, m) = (ρc2 * Uᵛ[ni])[μ];
			for (int n = 0; n <= dimensionCount; n++) {
				dUdV(n + 1, m) = T(n, ti)[μ];
				dFdV(n + 1, m) = T(n, ni)[μ];
			}
			dUdV(τi, m) -= dUdV(0, m);
			dFdV(τi, m) -= dFdV(0, m);
		}
		auto dFdU = dFdV * inverse(dUdV);
		return dFdU;
	}
	auto eigenSystem(EquationOfState<Type> const &eos, int ni) const {
		enableFPE();
		using AutoMassDensity = FwdAutoDiff<MassDensityType<Type>, 1, std::tuple<MassDensityType<Type>, SpecificEnergyType<Type>>>;
		using AutoSpecificEnergy = FwdAutoDiff<SpecificEnergyType<Type>, 1, std::tuple<MassDensityType<Type>, SpecificEnergyType<Type>>>;
		using index_type = typename AutoMassDensity::index_type;
		Vector<DimensionlessType<Type>, transverseCount> βt;
		auto ti = transverseIndices[ni];
		auto const ϱ = AutoMassDensity::template independent<0>(ρ);
		auto const ϵ = AutoSpecificEnergy::template independent<1>(ε);
		auto const ᵱ = eos.energy2pressure(ϱ, ϵ);
		auto const p = ᵱ.template multiGet<index_type::zero()>();
		auto const χ = ᵱ.template multiGet<index_type::unit(0)>();
		auto const κ = ᵱ.template multiGet<index_type::unit(1)>() / ρ;
		auto const h = 1 + ε / c2 + p / (ρ * c2);
		auto const α2 = (χ + (p / ρ) * κ) / (c2 * h);
		auto const α = sqrt(α2);
		auto const v = fourVelocity2CoordVelocity(u);
		auto const β = v / c;
		auto const β2 = β.dot(β);
		auto const βx = β[ni];
		for (int k = 0; k < transverseCount; k++) {
			βt[k] = β[ti[k]];
		}
		auto const βxβx = sqr(βx);
		auto const βtβt = βt.dot(βt);
		auto const β2α2 = β2 * α2;
		auto const id = 1 / (1 - β2α2);
		auto const n1 = βx * (1 - α2);
		auto const n2 = α * sqrt((1 - β2) * (1 - βxβx - βtβt * α2));
		auto const λp = (n1 + n2) * id;
		auto const λm = (n1 - n2) * id;
		auto const λo = βx;
		auto const Ƙ = κ / (κ - α2);
		auto const Ap = (1 - βxβx) / (1 - βx * λp);
		auto const Am = (1 - βxβx) / (1 - βx * λm);
		auto const γ = u[dimensionCount] / c;
		auto const γ2 = sqr(γ);
		SquareMatrix<Type, fieldCount> R;
		Vector<Type, fieldCount> λ;
		ni++;
		for (auto &i : ti) {
			i++;
		}
		R(Di, Di) = 1;
		R(ni, Di) = Type(h * γ * Am * λm);
		R(τi, Di) = Type(h * γ * Am - 1);
		R(Di, ni) = Type(Ƙ / (h * γ));
		R(ni, ni) = Type(βx);
		R(τi, ni) = Type(1 - Ƙ / (h * γ));
		R(Di, τi) = 1;
		R(ni, τi) = Type(h * γ * Ap * λp);
		R(τi, τi) = Type(h * γ * Ap - 1);
		for (int k = 0; k < transverseCount; k++) {
			R(ti[k], Di) = Type(h * γ * βt[k]);
			R(ti[k], ni) = Type(βt[k]);
			R(ti[k], τi) = Type(h * γ * βt[k]);
			R(Di, ti[k]) = Type(γ * βt[k]);
			R(τi, ti[k]) = Type((2 * h * γ - 1) * γ * βt[k]);
			R(ni, ti[k]) = Type(2 * h * γ2 * βt[k] * βx);
			for (int j = 0; j < transverseCount; j++) {
				R(ti[j], ti[k]) = Type(h * (2 * γ2 * βt[j] * βt[k] + Type(j == k)));
			}
		}
		for (int k = 0; k < fieldCount; k++) {
			R.setColumn(k, normalize(R.getColumn(k)));
		}
		λ[Di] = Type(λm);
		λ[τi] = Type(λp);
		λ[ni] = Type(λo);
		for (int k = 0; k < transverseCount; k++) {
			λ[ti[k]] = Type(λo);
		}
		return std::tuple(λ, R);

	}
	void setMassDensity(MassDensityType<Type> const &ρ_) {
		ρ = ρ_;
	}
	void setSpecificEnergy(SpecificEnergyType<Type> const &ε_) {
		ε = ε_;
	}
	void setVelocity(Vector<VelocityType<Type>, dimensionCount> const &v) {
		u = coordVelocity2FourVelocity(v);
	}
	void setVelocity(int k, VelocityType<Type> const &vk) {
		auto v = fourVelocity2CoordVelocity(u);
		v[k] = vk;
		u = coordVelocity2FourVelocity(v);
	}
	MassDensityType<Type> getMassDensity() const {
		return ρ;
	}
	SpecificEnergyType<Type> getSpecificEnergy() const {
		return ε;
	}
	Vector<VelocityType<Type>, dimensionCount> getVelocity() const {
		return fourVelocity2CoordVelocity(u);
	}
	VelocityType<Type> getVelocity(int k) const {
		return fourVelocity2CoordVelocity(u)[k];
	}
	friend GasConserved<Type, dimensionCount> ;
	template<typename, int>
	friend struct RadConserved;
	template<typename T, int D>
	friend auto hllc(GasConserved<T, D> const&, GasConserved<T, D> const&, EquationOfState<T> const&, int);
private:
	// @formatter:off
	static constexpr auto transverseIndices =
			(dimensionCount == 1) ? std::array<std::array<int, 2>, 3> {{{-1,-1}, {-1,-1}, {-1,-1}}} :
		   ((dimensionCount == 2) ? std::array<std::array<int, 2>, 3> {{{ 1,-1}, { 0,-1}, {-1,-1}}} :
		   /*dimensionCount == 3)*/	std::array<std::array<int, 2>, 3> {{{ 1, 2}, { 0, 2}, { 0, 1}}});
																			// @formatter:on
	static constexpr int Di = 0;
	static constexpr int τi = 1 + dimensionCount;
	static constexpr int transverseCount = dimensionCount - 1;
	static constexpr int fieldCount = 2 + dimensionCount;
	static constexpr auto η = minkowski<Type, dimensionCount>();
	static constexpr auto c = PhysicalConstants<Type>::c;
	static constexpr auto ic = 1 / c;
	static constexpr auto c2 = sqr(c);
	static constexpr auto ic2 = sqr(ic);
	MassDensityType<Type> ρ;
	SpecificEnergyType<Type> ε;
	Vector<VelocityType<Type>, dimensionCount + 1> u;
};

#endif /* INCLUDE_SRHD_PRIMITIVE_HPP_ */
