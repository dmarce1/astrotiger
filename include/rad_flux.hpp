/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "constants.hpp"
#include "left_right_state.hpp"
#include "operators.hpp"
#include "rad_conserved.hpp"
#include "vector.hpp"

template<typename Type, int dimensionCount>
struct RadFlux {
	EnergyFluxType<Type> E;
	Vector<FluxFluxType<Type>, dimensionCount> F;
	// @formatter:off
	DEFINE_VECTOR_OPERATORS(RadFlux, Type, a.E, a.F)
	template<typename Type1, int dimensionCount1>
	friend auto riemannHLLC(RadConserved<Type1, dimensionCount1> const &Ul, RadConserved<Type1, dimensionCount1> const &Ur, int n) ;
																																			// @formatter:on
};

template<typename Type, int dimensionCount>
auto flux(RadConserved<Type, dimensionCount> const &rad, int ni) {
	constexpr Type t1 = sqrt(std::numeric_limits < Type > ::min());
	constexpr PhysicalConstants<Type> pc { };
	enableFPE();
	using std::max;
	constexpr auto c2 = sqr(pc.c);
	RadFlux<Type, dimensionCount> flux;
	auto const &E = rad.E;
	auto const &F = rad.F;
	auto const un = unitVector<Type, dimensionCount>(ni);
	auto const magF = sqrt(F.dot(F));
	auto const imagF = 1 / max(magF, EnergyFluxType<Type>(t1));
	auto const iE = 1 / E;
	auto const f = magF * iE / pc.c;
	auto const N = F * imagF;
	auto const f2 = sqr(f);
	auto const ξ = (3 + 4 * f2) / (5 + 2 * sqrt(4 - 3 * f2));
	auto const Ωdif = 0.5 * (1 - ξ);
	auto const Ωstr = 0.5 * (3 * ξ - 1);
	auto const P = E * (Ωstr * N * N[ni] + Ωdif * un);
	flux.E = F[ni];
	for (int k = 0; k < dimensionCount; k++) {
		flux.F = c2 * P[ni][k];
	}
	return flux;
}

template<typename T>
using LRS = LeftRightState<T>;

template<typename Type, int dimensionCount>
auto riemannHLLC(RadConserved<Type, dimensionCount> const &Ul, RadConserved<Type, dimensionCount> const &Ur, int n) {
	using std::abs;
	using std::min;
	using std::max;
	constexpr PhysicalConstants<Type> pc { };
	constexpr VelocityType<Type> zero(Type(0.0));
	constexpr LeftRightState<DimensionlessType<Type>> one(Type(1.0));
	constexpr LeftRightState<DimensionlessType<Type>> two(Type(2.0));
	constexpr LeftRightState<DimensionlessType<Type>> three(Type(3.0));
	constexpr LeftRightState<DimensionlessType<Type>> four(Type(4.0));
	constexpr LeftRightState<DimensionlessType<Type>> five(Type(5.0));
	constexpr auto half = one / two;
	constexpr auto third = one / three;
	constexpr auto c = pc.c;
	constexpr auto c2 = sqr(c);
	RadFlux<Type, dimensionCount> flux;
	Vector<EnergyFluxType<Type>, dimensionCount> F0, F_;
	auto const E = LRS<EnergyDensityType<Type>>(Ul.E, Ur.E);
	auto const F = LRS<Vector<EnergyFluxType<Type>, dimensionCount>>(Ul.F, Ur.F);
	LRS<Vector<DimensionlessType<Type>, dimensionCount>> const f = F / (c * E);
	auto const F2 = F.dot(F);
	LRS<DimensionlessType<Type>> const f2 = f.dot(f);
	auto const f1 = sqrt(f2);
	auto const s = sqrt(four - three * f2);
	auto const ξ = (three + four * f2) / (five + two * s);
	auto const Fx = F[n];
	auto const βx = c2 * half * (3 * ξ - 1) * Fx * E / F2;
	auto const Π = half * (one - ξ) * E;
	auto const μ = Fx / sqrt(F2);
	auto const ζ = sqrt(two * third * (sqr(s) - s) + two * sqr(μ) * (two - f2 - s));
	auto const λₘ = c * (f1 * μ - ζ) / s;
	auto const λₚ = c * (f1 * μ + ζ) / s;
	auto const λ = LRS<VelocityType<Type>>(min(min(λₘ.L, λₘ.R), zero), max(max(λₚ.L, λₚ.R), zero));
	auto const A = λ * E - Fx;
	auto const B = (λ - βx) * Fx / c2 - Π;
	auto const α = (A.L - A.R) - (B.R * λ.L - B.L * λ.R);
	auto const ω = 2 * (A.L * λ.R - A.R * λ.L) / c;
	auto const δ = sqrt(2 * c * (B.R - B.L) * ω + sqr(α));
	auto const λ1 = c * (α - δ) / ω;
	auto const λ2 = c * (α + δ) / ω;
	auto const λ0 = (abs(λ1) < abs(λ2)) ? λ1 : λ2;
	bool const b = zero > λ0;
	auto const A_ = A(b);
	auto const B_ = B(b);
	auto const λ_ = λ(b);
	auto const βx_ = βx(b);
	auto const E_ = E(b);
	auto const Π_ = Π(b);
	for (int d = 0; d < dimensionCount; d++) {
		F_[d] = F[d](b);
	}
	auto const Π0 = (A_ * λ0 - c2 * B_) / (c2 - λ_ * λ0);
	auto const Ω = (λ_ - βx_) / (λ_ - λ0);
	auto E0 = E_ * Ω;
	for (int d = 0; d < dimensionCount; d++) {
		F0[d] = F[d](b) * Ω;
	}
	auto &F0x = F0[n];
	E0 = E0 - (βx_ * Π_ - Π0 * λ0) / (λ_ - λ0);
	F0x = F0x - c2 * (Π_ - Π0) / (λ_ - λ0);
	auto const f0 = F0 / (c * E0);
	auto const f02 = f0.dot(f0);
	auto const F02 = F0.dot(F0);
	auto const s0 = sqrt(4 - 3 * f02);
	auto const ξ0 = (3 + 4 * f02) / (5 + 2 * s0);
	auto const βx0 = c2 * 0.5 * (3 * ξ0 - 1) * F0x * E0 / F02;
	flux.E = F0x;
	flux.F = F0 * βx0;
	flux.F[n] += c2 * Π0;
	return flux;
}
