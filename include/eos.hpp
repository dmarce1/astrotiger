/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_EOS_EOS_HPP_
#define INCLUDE_EOS_EOS_HPP_

#include "constants.hpp"
#include "cgs.hpp"
#include "rational.hpp"

template<typename Type>
struct EquationOfState {
	constexpr EquationOfState() = default;
	constexpr EquationOfState(EquationOfState const&) = default;
	constexpr EquationOfState& operator=(EquationOfState const&) = default;
	constexpr EquationOfState(MolarMassType<Type> μ_) :
			μ(μ_) {
	}
	template<typename Tρ, typename Tε>
	auto ε2p(Tρ const &ρ, Tε const &ε) const {
		return (ρ * ε) * Type(Γ - 1);
	}
	template<typename Tρ, typename Tε>
	auto ε2χ(Tρ const &ρ, Tε const &ε) const {
		return ε * Type(Γ - 1);
	}
	template<typename Tρ, typename Tε>
	auto ε2κ(Tρ const &ρ, Tε const &ε) const {
		return Type(Γ - 1);
	}
	template<typename Tρ, typename Tε>
	auto ε2T(Tρ const &ρ, Tε const &ε) const {
		return (Γ - 1) * μ * ε / (pc.nₐ * pc.kB);
	}
	template<typename Tρ, typename TT>
	auto T2s(Tρ const &ρ, TT const &T) const {
		constexpr Rational Cv = 1 / (Γ - 1);
		constexpr Rational Cp = Γ / (Γ - 1);
		return pc.nₐ * pc.kB / μ * (log(μ / (pc.nₐ * ρ) * pow<Cv>((4 * pc.π * Cv * pc.kB * T * μ) / (3 * pc.h * pc.h * pc.nₐ))) + Cp);
	}
	template<typename Tρ, typename TT>
	auto T2ε(Tρ const &ρ, TT const &T) const {
		return (pc.nₐ * pc.kB) / ((Γ - 1) * μ) * T;
	}
	template<typename Tρ, typename Tϰ>
	auto s2T(Tρ const &ρ, Tϰ const &ϰ) const {
		constexpr Rational Cv = 1 / (Γ - 1);
		constexpr Rational Cp = Γ / (Γ - 1);
		auto const term1 = ((3 * pc.h * pc.h * pc.nₐ) / (4 * pc.π * Cv * pc.kB * μ));
		auto const term2 = pow<1 / Cv>(pc.nₐ * ρ / μ * exp(μ * ϰ / (pc.nₐ * pc.kB)) / exp(Cp));
		return term1 * term2;
	}
private:
	static constexpr PhysicalConstants<Type> pc { };
	static constexpr Rational Γ = Rational(5, 3);
	MolarMassType<Type> const μ;
};

#endif /* INCLUDE_EOS_EOS_HPP_ */
