#pragma once

#include "Definitions.hpp"
#include "Quadrature.hpp"

#include <numbers>

std::vector<Real> getLagrangePolynomial(int j, int N);

template<typename Type, int modeCount, int nodeCount>
constexpr auto legendreTransform() {
	constexpr auto quad = getLegendreQuadraturePoints<Type, nodeCount>();
	constexpr auto qPoint = quad.first;
	constexpr auto qWeight = quad.second;
	std::array<std::array<Type, nodeCount>, modeCount> transform;
	for (int mi = 0; mi < modeCount; mi++) {
		for (int ni = 0; ni < nodeCount; ni++) {
			transform[mi][ni] = qWeight[ni] * legendre<Type>(mi, qPoint[ni]) * Type(2 * mi + 1) / Type(2);
		}
	}
	return transform;
}

template<typename Type, int dimCount, int modeCount, int nodeCount>
auto nodalToModal(std::array<Type, power(nodeCount, dimCount)> x) {
	constexpr auto lT = legendreTransform<Type, modeCount, nodeCount>();
	std::array<Type, power(nodeCount, dimCount)> y;
	std::array<Type, binomialCoefficient(modeCount + dimCount - 1, dimCount)> z;
	int blockCount = power(nodeCount, dimCount - 1);
	int loSize = 1; //binomialCoefficient(modeCount + cdim - 1, cdim);
	for (int dim = dimCount - 1; dim >= 0; dim--) {
		int oi = 0;
		int const cdim = dimCount - 1 - dim;
		for (int hi = 0; hi < blockCount; hi++) {
				for (int mi = 0; mi < modeCount; mi++) {
				int loMax = binomialCoefficient(modeCount + cdim - mi - 1, cdim);
				for (int lo = 0; lo < loMax; lo++) {
					y[oi] = Type(0);
					for (int ni = 0; ni < nodeCount; ni++) {
						int const xIndex = (hi * nodeCount + ni) * loSize + lo;
						y[oi] += lT[mi][ni] * x[xIndex];
					}
					oi++;
				}
			}
		}
		loSize *= modeCount + cdim;
		loSize /= cdim + 1;
		blockCount /= nodeCount;
		std::swap(x, y);
	}
	std::copy(x.begin(), x.begin() + z.size(), z.begin());
	return z;
}
