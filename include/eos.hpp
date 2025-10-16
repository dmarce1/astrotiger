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
	decltype(Tε() * Tρ()) energy2pressure(Tρ const &ρ, Tε const &ε) const {
		return (ρ * ε) * (Γ - 1);
	}
	template<typename Tρ, typename Tε>
	auto energy2temperature(Tρ const &ρ, Tε const &ε) const {
		using namespace Constants;
		return (Γ - 1) * μ * ε / (nₐ * kB);
	}
	template<typename Tρ, typename TT>
	auto temperature2entropy(Tρ const &ρ, TT const &T) const {
		using namespace Constants;
		constexpr Rational Cv = 1 / (Γ - 1);
		constexpr Rational Cp = Γ / (Γ - 1);
		return nₐ * kB / μ * (log(μ / (nₐ * ρ) * pow<Cv>((4 * π * Cv * kB * T * μ) / (3 * _h * _h * nₐ))) + Cp);
	}
	template<typename Tρ, typename Tϰ>
	auto entropy2temperature(Tρ const &ρ, Tϰ const &ϰ) const {
		using namespace Constants;
		constexpr Rational Cv = 1 / (Γ - 1);
		constexpr Rational Cp = Γ / (Γ - 1);
		auto const term1 = ((3 * _h * _h * nₐ) / (4 * π * Cv * kB * μ));
		auto const term2 = pow<1 / Cv>(nₐ * ρ / μ * exp(μ * ϰ / (nₐ * kB)) / exp(Cp));
		return term1 * term2;
	}
	template<typename Tρ, typename TT>
	auto temperature2energy(Tρ const &ρ, TT const &T) const {
		using namespace Constants;
		return (nₐ * kB) / ((Γ - 1) * μ) * T;
	}
private:
	static constexpr Rational Γ = Rational(5, 3);
	MolarMassType<Type> const μ;
};



#endif /* INCLUDE_EOS_EOS_HPP_ */
