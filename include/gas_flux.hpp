/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "eos.hpp"
#include "gas_conserved.hpp"
#include "gas_primitive.hpp"
#include "left_right_state.hpp"

template<typename Type, int dimensionCount>
struct GasFlux {
	MomentumDensityType<Type> D;
	Vector<EnergyDensityType<Type>, dimensionCount> S;
	EnergyFluxType<Type> τ;
	EntropyFluxType<Type> K;
	// @formatter:off
	DEFINE_VECTOR_OPERATORS(GasFlux, Type, a.D, a.τ, a.K, a.S);
						// @formatter:on
};

template<typename Type, int dimensionCount>
auto flux(GasPrimitive<Type, dimensionCount> const &prim, EquationOfState<Type> const &eos, int n) {
	constexpr PhysicalConstants<Type> pc { };
	constexpr auto c = pc.c;
	constexpr auto c2 = sqr(c);
	GasFlux<Type, dimensionCount> flux;
	auto const &ρ = prim.ρ;
	auto const &ε = prim.ε;
	auto const &v = prim.v;
	auto const β = v / c;
	auto const β2 = β.dot(β);
	auto const γ = sqrt(1 / (1 - β2));
	auto const γ2 = sqr(γ);
	auto const vˣ = v[n];
	auto const T = eos.ε2T(ρ, ε);
	auto const κ = eos.T2s(ρ, T);
	auto const p = eos.ε2p(ρ, ε);
	auto const h = c2 + ε + p / ρ;
	auto const W = ρ * γ2 * h;
	flux.D = γ * ρ * vˣ;
	flux.K = γ * ρ * κ * vˣ;
	flux.S = W * v * vˣ / c2;
	flux.τ = (γ2 * (ρ * ε + β2 * (c2 * ρ * γ / (γ + 1) + p))) * vˣ;
	flux.S[n] += p;
	return flux;
}

template<typename Type, int dimensionCount>
auto riemannHLLC(GasConserved<Type, dimensionCount> const &UL, GasConserved<Type, dimensionCount> const &UR, EquationOfState<Type> const &eos, int n) {
	using std::max;
	using std::min;
	using PrimType = GasPrimitive<Type, dimensionCount>;
	using ConsType = GasConserved<Type, dimensionCount>;
	constexpr VelocityType<Type> zero(Type(0.0));
	constexpr PhysicalConstants<Type> pc { };
	constexpr LRS<DimensionlessType<Type>> one(Type(1.0));
	constexpr LRS<DimensionlessType<Type>> two(Type(2.0));
	constexpr LRS<DimensionlessType<Type>> three(Type(3.0));
	constexpr LRS<DimensionlessType<Type>> four(Type(4.0));
	constexpr LRS<DimensionlessType<Type>> five(Type(5.0));
	constexpr auto c = pc.c;
	constexpr auto c2 = sqr(c);
	ConsType U0;
	PrimType const primL = UL.toPrimitive(eos);
	PrimType const primR = UR.toPrimitive(eos);
	auto const U = LRS<ConsType>(UL, UR);
	auto const V = LRS<PrimType>(primL, primR);
	auto const D = LRS<MassDensityType<Type>>(U.L.D, U.R.D);
	auto const S = LRS<Vector<MomentumDensityType<Type>, dimensionCount>>(U.L.S, U.R.S);
	auto const τ = LRS<EnergyDensityType<Type>>(U.L.τ, U.R.τ);
	auto const K = LRS<EntropyDensityType<Type>>(U.L.K, U.R.K);
	auto const ρ = LRS<MassDensityType<Type>>(V.L.ρ, V.R.ρ);
	auto const v = LRS<Vector<VelocityType<Type>, dimensionCount>>(V.L.v, V.R.v);
	auto const ε = LRS<SpecificEnergyType<Type>>(V.L.ε, V.R.ε);
	auto const vx = v[n];
	auto const p = eos.ε2p(ρ, ε);
	auto const χ = eos.ε2χ(ρ, ε);
	auto const κ = eos.ε2κ(ρ, ε);
	auto const h = one + ε / c2 + p / (ρ * c2);
	auto const α2 = (χ + (p / ρ) * κ) / (c2 * h);
	auto const α = sqrt(α2);
	auto const β = v / c;
	auto const β2 = β.dot(β);
	auto const βx = β[n];
	auto const βxβx = sqr(βx);
	auto const βtβt = β2 - βxβx;
	auto const β2α2 = β2 * α2;
	auto const id = one / (one - β2α2);
	auto const n1 = βx * (one - α2);
	auto const n2 = α * sqrt((one - β2) * (one - βxβx - βtβt * α2));
	auto const λp = c * (n1 + n2) * id;
	auto const λm = c * (n1 - n2) * id;
	auto const γ2 = one / (one - β2);
	auto const γ = sqrt(γ2);
	auto const W = ρ * sqr(γ) * h;
	auto const s = LRS<VelocityType<Type>>(min(zero, min(λm.L, λm.R)), max(zero, max(λp.L, λp.R)));
	auto const s0 = (p.R - p.L + s.L * W.L * vx.L - s.R * W.R * vx.R) / (W.L * (s.L - vx.L) - W.R * (s.R - vx.R));
	auto const p0 = p.L + c2 * W.L * (s.L - vx.L) * (s0 - vx.L) / (c2 - vx.L * s.L);
	bool const b = zero > s0;
	auto const Ξ = (s(b) - vx(b)) / (s(b) - s0);
	U0.D = D(b) * Ξ;
	U0.K = K(b) * Ξ;
	U0.S = S(b) * Ξ;
	U0.τ = τ(b) * Ξ;
	U0.S[n] += (p0 - p(b)) / (s(b) - s0);
	U0.τ += (p0 * s0 - p(b) * vx(b)) / (s(b) - s0);
	auto F = flux(V(b), eos, n);
//	F += s(b) * (U0 - U(b));
//	return F;
}

