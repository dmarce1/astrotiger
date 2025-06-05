#pragma once

#include "Util.hpp"

static constexpr auto maximumSimdWidth = SIMD_WIDTH;
static constexpr auto dimensionCount = DIMENSION_COUNT;
static constexpr auto modeCount = MODE_COUNT;
static constexpr auto nodeCount = MODE_COUNT + 1;
static constexpr auto nodeCount3d = ipow(nodeCount, dimensionCount);
static constexpr auto modeCount3d = binco(modeCount + dimensionCount - 1, dimensionCount);

