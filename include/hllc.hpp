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
	GasConserved<Type, dimensionCount> U0L, U0R;
	auto const primL = UL.toPrimitive(eos);
	auto const primR = UR.toPrimitive(eos);
	auto const ρL = primL.ρ;
	auto const ρR = primR.ρ;
	auto const uL = primL.v;
	auto const uR = primR.v;
	auto const εL = primL.ε;
	auto const εR = primR.ε;
	auto const FR = primR.flux(direction);
	auto const FL = primL.flux(direction);
	auto const γL = primL.lorentzFactor();
	auto const γR = primR.lorentzFactor();
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
	auto const λmL = c2 * (vL - aL) / (c2 - vL * aL);
	auto const λmR = c2 * (vR - aR) / (c2 - vR * aR);
	auto const λpL = c2 * (vL + aL) / (c2 + vL * aL);
	auto const λpR = c2 * (vR + aR) / (c2 + vR * aR);
	auto const sL = min(zero, min(λmL, λmR));
	auto const sR = max(zero, max(λpL, λpR));
	auto const s0num = pR - pL + sL * wL * uL - sR * wR * uR;
	auto const s0den = wL * (sL - uL) - wR * (sR - uR);
	auto const s0 = s0num / s0den;
	auto const p0 = pL + wL * (sL - uL) * (s0 - uL) / (1 - uL * sL);
	U0R.D = UR.D * (sR - uR) / (sR - s0);
	U0L.D = UL.D * (sL - uL) / (sL - s0);
	U0R.S = UR.S * (U0R.D / UR.D);
	U0L.S = UL.S * (U0L.D / UL.D);
	U0R.S[direction] = U0R.D * hR * sqr(wR) * s0;
	U0L.S[direction] = U0L.D * hL * sqr(wL) * s0;
	U0R.τ = U0R.D * hR - p0 - U0R.D;
	U0L.τ = U0L.D * hL - p0 - U0L.D;
	if (zero <= sL) {
		return FL;
	} else if (s0 >= zero) {
		return FL + sL * (U0L - UL);
	} else if (sR >= zero) {
		return FR + sR * (U0R - UR);
	} else {
		return FR;
	}
}
