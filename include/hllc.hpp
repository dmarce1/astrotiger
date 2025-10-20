/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "conserved.hpp"
#include "primitive.hpp"
#include "eos.hpp"

template<typename Type, int dimensionCount>
auto hllc(GasConserved<Type, dimensionCount> const &UL, GasConserved<Type, dimensionCount> const &UR, EquationOfState<Type> const &eos, int direction) {
	using std::max;
	using std::min;
	constexpr VelocityType<Type> zero(Type(0.0));
	constexpr PhysicalConstants<Type> pc { };
	constexpr Type c = pc.c;
	constexpr Type c2 = sqr(c);
	GasFlux<Type, dimensionCount> F;
	GasConserved<Type, dimensionCount> U0;
	auto const primL = UL.toPrimitive(eos);
	auto const primR = UR.toPrimitive(eos);
	auto const ρL = primL.ρ;
	auto const ρR = primR.ρ;
	auto const uL = primL.v;
	auto const uR = primR.v;
	auto const εL = primL.ε;
	auto const εR = primR.ε;
	auto const γL = primL.lorentzFactor();
	auto const γR = primR.lorentzFactor();
	auto const λL = primL.eigenvalues(eos, direction);
	auto const λR = primR.eigenvalues(eos, direction);
	auto const hL = eos.energy2enthalpy(ρL, εL);
	auto const hR = eos.energy2enthalpy(ρR, εR);
	auto const wL = ρL * sqr(γL) * hL;
	auto const wR = ρR * sqr(γR) * hR;
	auto const pL = eos.energy2pressure(ρL, εL);
	auto const pR = eos.energy2pressure(ρR, εR);
	auto const aL = eos.energy2soundSpeed(ρL, εL);
	auto const aR = eos.energy2soundSpeed(ρR, εR);
	auto const vL = uL[direction];
	auto const vR = uR[direction];
	auto const sL = min(zero, min(λL.front(), λR.front()));
	auto const sR = max(zero, max(λL.back(), λR.back()));
	auto const s0num = pR - pL + sL * wL * uL - sR * wR * uR;
	auto const s0den = wL * (sL - uL) - wR * (sR - uR);
	auto const s0 = s0num / s0den;
	auto const p0 = pL + wL * (sL - uL) * (s0 - uL) / (1 - uL * sL);
	if (zero >= s0) {
		F = primL.flux(direction);
		if (zero >= sL) {
			U0.D = UL.D * (sL - uL) / (sL - s0);
			U0.K = UL.K * (U0.D / UL.D);
			U0.S = UL.S * (U0.D / UL.D);
			U0.τ = U0.D * hL - p0 - U0.D;
			U0.S[direction] = U0.D * hL * sqr(wL) * s0;
			F += sL * (U0 - UL);
		}
	} else {
		F = primR.flux(direction);
		if (zero >= sR) {
			U0.D = UR.D * (sL - uL) / (sL - s0);
			U0.K = UR.K * (U0.D / UR.D);
			U0.S = UR.S * (U0.D / UR.D);
			U0.τ = U0.D * hR - p0 - U0.D;
			U0.S[direction] = U0.D * hR * sqr(wR) * s0;
			F += sR * (U0 - UR);
		}
	}
	return F;
}
