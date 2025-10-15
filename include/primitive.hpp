/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
*******************************************************************************/
#ifndef INCLUDE_SRHD_PRIMITIVE_HPP_
#define INCLUDE_SRHD_PRIMITIVE_HPP_

#include "srhd.hpp"




template<int NDIM>
struct GasPrimitive {
#define XXX(name, type, arg, dim) std::conditional_t<dim == 1, Tensor0<type, dim>, Tensor1<type, dim>> name;
	GAS_PRIMITIVE_PROPERTIES(0, NDIM)
#undef XXX
	auto dualEnergySwitch(EquationOfState const &eos) const {
		using namespace Constants;
		auto const h = c2 + ε + eos.energy2pressure(ρ, ε) / ρ;
		auto const num = ρ * ε;
		auto const den = ρ * γ * (γ * h - c2);
		return num / den;
	}

	void updateLorentz() {
		auto const β2 = β.dot(β);
		γ = DimensionlessType(1.0) / sqrt(DimensionlessType(1.0) - β2);
	}
	GasConserved<NDIM> toConserved(EquationOfState const &eos) const {
#define XXX(name, type, arg, dim) CopyConstType<decltype(arg), std::conditional_t<dim == 1, Tensor0<type, dim>, Tensor1<type, dim>>>& name = arg.name;
		using namespace Constants;
		GasConserved<NDIM> con;
		GAS_CONSERVED_PROPERTIES(con, NDIM);
#undef XXX
		auto const β2 = double(β.dot(β));
		auto const γ2 = 1 / (1 - β2);
		auto const γ = sqrt(γ2);
		auto const iρ = 1 / ρ;
		auto const p = eos.energy2pressure(ρ, ε);
		auto const T = eos.energy2temperature(ρ, ε);
		auto const ϰ = eos.temperature2entropy(ρ, T);
		auto const h = c2 + ε + p * iρ;
		S = ρ * γ2 * h * β * ic;
		D = ρ * γ;
		τ = ρ * γ2 * ε + (γ2 - 1) * p + c2 * D * (γ - 1);
		K = D * ϰ;
		return con;
	}
	auto temperature(EquationOfState const &eos) const {
		return eos.energy2temperature(ρ, ε);
	}
};


#endif /* INCLUDE_SRHD_PRIMITIVE_HPP_ */
