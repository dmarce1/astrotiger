/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_GAS_FLUX_HPP_
#define INCLUDE_GAS_FLUX_HPP_

template<typename Type, int dimensionCount>
struct GasFlux {
	MomentumDensityType<Type> D;
	EnergyDensityType<Type> S;
	EnergyFluxType<Type> τ;
	EntropyFluxType<Type> K;
	// @formatter:on
	DEFINE_VECTOR_OPERATORS(GasFlux, Type, a.D, a.τ, a.K, a.S);
	// @formatter:off
};





#endif /* INCLUDE_GAS_FLUX_HPP_ */
