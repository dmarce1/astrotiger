/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_CONSTANTS_HPP_
#define INCLUDE_CONSTANTS_HPP_

#include "cgs.hpp"
#include "matrix.hpp"

namespace Constants {

static constexpr double π = 4 * atan(1);
static constexpr auto amu = MassType(1.0 / 6.02214076e23);
static constexpr auto nₐ = PerMoleType(6.02214076e23);
static constexpr auto c = VelocityType(2.99792458e10);
static constexpr auto kB = EntropyType(1.380658e-16);
static constexpr auto aR = EntropyDensityPerKelvinCubedType(7.5657e-15);
static constexpr auto _h = TorqueType(6.62607015e-27);
static constexpr auto ℏ = _h / (2 * π);
static constexpr auto c2 = c * c;
static constexpr auto ic = 1.0 / c;
static constexpr auto ic2 = ic * ic;
static constexpr double tiny = sqrt(std::numeric_limits<double>::min());
static constexpr double huge = sqrt(std::numeric_limits<double>::max());
static constexpr double eps = std::numeric_limits<double>::epsilon();
template<int D>
static constexpr auto δ = identity<DimensionlessType, D>();

}
#endif /* INCLUDE_CONSTANTS_HPP_ */
