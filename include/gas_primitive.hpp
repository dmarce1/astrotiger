/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/
#ifndef INCLUDE_SRHD_PRIMITIVE_HPP_
#define INCLUDE_SRHD_PRIMITIVE_HPP_

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
	auto eigenvalues(EquationOfState<Type> const &eos, int dim) const {
		constexpr int fieldCount = 2 + dimensionCount;
		Vector<VelocityType<Type>, fieldCount> λ;
		auto const a = eos.energy2soundSpeed(ρ, ε);
		auto const α = a * ic;
		auto const v = fourVelocity2CoordVelocity(u);
		auto const β = v * ic;
		auto const βˣ = β[dim];
		auto const β2 = β.dot(β);
		auto const α2 = sqr(α);
		auto const α2β2 = α2 * β2;
		auto const a_eff = c * α * (sqrt((1 - β2) * (1 - α2β2 - (1 - α2) * sqr(βˣ)))) / (1 - α2β2);
		auto const v_eff = c * βˣ / (1 - α2β2);
		λ.front() = v_eff - a_eff;
		for (int k = 1; k + 1 < fieldCount; k++) {
			λ[k] = c * βˣ;
		}
		λ.back() = v_eff + a_eff;
		return λ;
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
		enableFPE();
		using Auto1 = FwdAutoDiff<Type, 1, std::array<Type, fieldCount>>;
		constexpr Auto1 c = PhysicalConstants<Type>::c.value();
		constexpr Auto1 c2 = sqr(c);
		Vector<Auto1, dimensionCount> ρcv, β;
		Auto1 const ρc2 = Auto1::independent(Type(ρ * c2), Di);
		auto const v = fourVelocity2CoordVelocity(u);
		for (int d = 0; d < dimensionCount; d++) {
			ρcv[d] = Auto1::independent(Type(ρ * c * v[d]), 1 + d);
		}
		Auto1 const ρε = Auto1::independent(Type(ρ * ε), τi);
		β = ρcv / ρc2;
		Auto1 const β2 = β.dot(β);
		Auto1 const γ2 = Auto1(1) / (Auto1(1) - β2);
		Auto1 const γ = sqrt(γ2);
		Auto1 const p = eos.energy2pressure(ρc2 / c2, ρε / ρc2 * c2);
		Auto1 const W = ρc2 * γ2 * (1 + (ρε + p) / ρc2);
		Vector<Auto1, fieldCount> U;
		U[Di] = (ρc2 * γ);
		U[τi] = W - p - ρc2 * γ;
		for (int d = 0; d < dimensionCount; d++) {
			U[1 + d] = W * β[d];
		}
		auto F = U * c * β[ni];
		F[1 + ni] += c * p;
		F[τi] += c * β[ni] * p;
		SquareMatrix<Type, fieldCount> dUdV;
		SquareMatrix<Type, fieldCount> dFdV;
		for (int n = 0; n < fieldCount; n++) {
			for (int m = 0; m < fieldCount; m++) {
				auto const μ = Auto1::index_type::unit(m);
				dUdV(n, m) = U[n][μ];
				dFdV(n, m) = F[n][μ];
			}
		}
		auto dFdU = dFdV * inverse(dUdV);
		return dFdU;
	}
	auto eigenSystem(EquationOfState<Type> const &eos, int ni) const {
		enableFPE();
		using AutoMassDensity = FwdAutoDiff<MassDensityType<Type>, 1, std::tuple<MassDensityType<Type>, SpecificEnergyType<Type>>>;
		using AutoSpecificEnergy = FwdAutoDiff<SpecificEnergyType<Type>, 1, std::tuple<MassDensityType<Type>, SpecificEnergyType<Type>>>;
		using index_type = typename AutoMassDensity::index_type;
		std::array<DimensionlessType<Type>, transverseCount> βt;
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
		auto const βtβt = β2 - βxβx;
		auto const β2α2 = β2 * α2;
		auto const id = 1 / (1 - β2α2);
		auto const n1 = βx * (1 - α2) * id;
		auto const n2 = α * sqrt((1 - β2) * (1 - βxβx - βtβt * α2));
		auto const λp = (n1 + n2) * id;
		auto const λm = (n1 - n2) * id;
		auto const λo = βx;
		auto const Ƙ = κ / (κ - α2);
		auto const Ap = (1 - βxβx) / (1 - βx * λp);
		auto const Am = (1 - βxβx) / (1 - βx * λm);
		auto const iγ2 = 1 - β2;
		auto const γ2 = 1 / iγ2;
		auto const γ = sqrt(γ2);
		SquareMatrix<Type, fieldCount> R;
		Vector<Type, fieldCount> λ;
		ni++;
		for (auto &i : ti) {
			i++;
		}
		R(Di, Di) = 1;
		R(τi, Di) = Type(h * γ * Am - 1);
		R(ni, Di) = Type(h * γ * Am * λm);
		R(Di, τi) = 1;
		R(τi, τi) = Type(h * γ * Ap - 1);
		R(ni, τi) = Type(h * γ * Ap * λp);
		R(Di, ni) = Type(Ƙ / (h * γ));
		R(ni, ni) = Type(βx);
		R(τi, ni) = Type(1 - Ƙ / (h * γ));
		for (int k = 0; k < transverseCount; k++) {
			R(ti[k], Di) = Type(h * γ * βt[k]);
			R(ti[k], τi) = Type(h * γ * βt[k]);
			R(ti[k], ni) = Type(βt[k]);
			R(Di, ti[k]) = Type(γ * βt[k]);
			R(τi, ti[k]) = Type(2 * h * γ2 * βt[k] - γ * βt[k]);
			R(ni, ti[k]) = Type(2 * h * γ2 * βt[k] * βx);
			for (int j = 0; j < transverseCount; j++) {
				R(ti[j], ti[k]) = Type(h * (2 * γ2 * βt[j] * βt[k] + Type(j == k)));
			}
		}
		λ[Di] = Type(c * λm);
		λ[τi] = Type(c * λp);
		λ[ni] = Type(c * λo);
		for (int k = 0; k < transverseCount; k++) {
			λ[ti[k]] = Type(c * λo);
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
	static constexpr auto c = PhysicalConstants<Type>::c;
	static constexpr auto ic = 1 / c;
	static constexpr auto c2 = sqr(c);
	static constexpr auto ic2 = sqr(ic);
	MassDensityType<Type> ρ;
	SpecificEnergyType<Type> ε;
	Vector<VelocityType<Type>, dimensionCount + 1> u;
};

#endif /* INCLUDE_SRHD_PRIMITIVE_HPP_ */
