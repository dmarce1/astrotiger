/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "constants.hpp"
#include "quantity.hpp"
#include "vector.hpp"

template<typename T, int D4>
auto fourVelocity2CoordVelocity(Vector<VelocityType<T>, D4> const &u) {
	constexpr int D3 = D4 - 1;
	auto const γ = fourVelocity2LorentzFactor(u);
	Vector<VelocityType<T>, D3> v;
	for (int k = 0; k < D3; k++) {
		v[k] = u[k] / γ;
	}
	return v;
}

template<typename T, int D4>
auto fourVelocity2LorentzFactor(Vector<VelocityType<T>, D4> const &u) {
	constexpr int D3 = D4 - 1;
	return u[D3] / PhysicalConstants<T>::c;
}

template<typename T, int D3>
auto coordVelocity2FourVelocity(Vector<VelocityType<T>, D3> const &v) {
	constexpr int D4 = D3 + 1;
	Vector<VelocityType<T>, D4> u;
	auto const γ = coordVelocity2LorentzFactor(v);
	for (int k = 0; k < D3; k++) {
		u[k] = γ * v[k];
	}
	u[D3] = γ * PhysicalConstants<T>::c;
	return u;
}

template<typename T, int D3>
auto coordVelocity2LorentzFactor(Vector<VelocityType<T>, D3> const &v) {
	constexpr auto c2 = sqr(PhysicalConstants<T>::c);
	auto const v2 = v.dot(v);
	auto const γ2 = c2 / (c2 - v2);
	return sqrt(γ2);
}

//₀₁₂₃₄₅₆₇₈₉ₐₑₒₓₕₖₗₘₙₚₛₜ
//⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ᴬᴮᴰᴱᴳᴴᴵᴶᴷᴸᴹᴺᴼᴾᴿᵀᵁᵂᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻ

