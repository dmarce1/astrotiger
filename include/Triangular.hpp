#pragma once

#include <array>

template<int dimensionCount>
int triangularToFlat(std::array<int, dimensionCount> const &triangleIndex) {
	int flatIndex = 0;
	for (int dimension = dimensionCount - 1; dimension >= 0; dimension--) {
		auto thisIndex = triangleIndex[dimension];
		thisIndex -= triangleIndex[dimension + 1];
		flatIndex += binomialCoefficient(thisIndex + dimensionCount - 1, dimensionCount);
	}
	return flatIndex;
}

template<int modeCount, int dimensionCount>
bool triangularIndexIncrement(std::array<int, dimensionCount> &triIndex) {
	for (int d = dimensionCount - 2; d >= 0; d--) {
		triIndex[d] += triIndex[d + 1];
	}
	int d = dimensionCount - 1;
	while (++triIndex[d] >= (d ? std::min(modeCount, triIndex[d - 1] + 1) : modeCount)) {
		if (d <= 0) {
			return false;
		}
		triIndex[d--] = 0;
	}
	for (int d = 0; d < dimensionCount - 1; d++) {
		triIndex[d] -= triIndex[d + 1];
	}
	return true;
}

template<int dimensionCount>
constexpr int triangularSize(int maximumOrder) {
	return binomialCoefficient(maximumOrder + dimensionCount - 1, dimensionCount);
}
