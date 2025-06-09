#pragma once

#include "Definitions.hpp"

#include <vector>

struct QuadraturePoint {
	Real x;
	Real w;
};

std::vector<QuadraturePoint> getLobattoQuadraturePoints(int);
std::vector<QuadraturePoint> getLegendreQuadraturePoints(int);
