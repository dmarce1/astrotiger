#pragma once

#include "Util.hpp"

#include <numbers>

using Real = REAL_DATA_TYPE;


namespace Global {
constexpr auto simdWidth = SIMD_WIDTH;
constexpr auto dimensionCount = DIMENSION_COUNT;
constexpr auto modeCount = MODE_COUNT;
constexpr auto nodeCount = MODE_COUNT;
constexpr auto modeCount3d = binomialCoefficient(modeCount + dimensionCount - 1, dimensionCount);
constexpr auto nodeCount3d = power(nodeCount, dimensionCount);
};

namespace Constants {
constexpr Real zero { Real(0) };
constexpr Real one { Real(1) };
constexpr Real two { Real(2) };
constexpr Real half = one / two;
constexpr Real pi = std::numbers::pi_v<Real>;
}

