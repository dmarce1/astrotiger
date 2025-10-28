/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_CONSTANTS_HPP_
#define INCLUDE_CONSTANTS_HPP_

#include "cgs.hpp"
#include "math.hpp"
#include "matrix.hpp"

template<typename Type>
struct PhysicalConstants {
	static constexpr Type π = acos(-Type(1));
	static constexpr PerMoleType<Type> nₐ = 6.02214076e23;
	static constexpr VelocityType<Type> c = 2.99792458e10;
	static constexpr EntropyType<Type> kB = 1.380658e-16;
	static constexpr TorqueType<Type> h = 6.62607015e-27;
	static constexpr GravityConstantType<Type> G = 6.67430e-8;
	static constexpr auto ℏ = h / (Type(2) * π);
	static constexpr auto σ = (Type(2) * pow<5>(π) * pow<4>(kB)) / (Type(15) * pow<3>(h) * sqr(c));
	static constexpr auto aR = Type(4) * σ / c;
};

#endif /* INCLUDE_CONSTANTS_HPP_ */
