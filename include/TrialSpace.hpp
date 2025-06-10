#pragma once

#include "Definitions.hpp"
#include "Quadrature.hpp"

#include <numbers>

template<typename Type, int modeCount, int nodeCount>
constexpr auto legendreAnalysis() {
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

template<typename Type, int modeCount, int nodeCount>
constexpr auto legendreSynthesis() {
	constexpr auto quad = getLegendreQuadraturePoints<Type, nodeCount>();
	constexpr auto qPoint = quad.first;
	std::array<std::array<Type, modeCount>, nodeCount> transform;
	for (int mi = 0; mi < modeCount; mi++) {
		for (int ni = 0; ni < nodeCount; ni++) {
			transform[ni][mi] = legendre<Type>(mi, qPoint[ni]);
		}
	}
	return transform;
}

template<typename Type, int dimCount, int modeCount, int nodeCount>
auto nodalToModal(std::array<Type, power(nodeCount, dimCount)> x) {
	constexpr auto lT = legendreAnalysis<Type, modeCount, nodeCount>();
	std::array<Type, power(nodeCount, dimCount)> y;
	std::array<Type, binomialCoefficient(modeCount + dimCount - 1, dimCount)> z;
	int blockCount = power(nodeCount, dimCount - 1);
	int loSize = 1;
	for (int dim = dimCount - 1; dim >= 0; dim--) {
		int oi = 0;
		int const cdim = dimCount - 1 - dim;
		for (int hi = 0; hi < blockCount; hi++) {
			int loMax = loSize;
			for (int mi = 0; mi < modeCount; mi++) {
				for (int lo = 0; lo < loMax; lo++) {
					y[oi] = Type(0);
					for (int ni = 0; ni < nodeCount; ni++) {
						int const xIndex = (hi * nodeCount + ni) * loSize + lo;
						y[oi] += lT[mi][ni] * x[xIndex];
					}
					oi++;
				}
				if (mi + 1 < modeCount) {
					loMax *= modeCount - mi - 1;
					loMax /= modeCount + cdim - mi - 1;
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

template<typename Type, int dimCount, int modeCount, int nodeCount>
auto modalToNodal(std::array<Type, binomialCoefficient(modeCount + dimCount - 1, dimCount)> z) {
	auto lT = legendreSynthesis<Type, modeCount, nodeCount>();
	std::array<Type, power(nodeCount, dimCount)> y;
	std::array<Type, power(nodeCount, dimCount)> x;
	std::copy(z.begin(), z.end(), y.begin());
	int blockCount = 1;
	int loSize = binomialCoefficient(modeCount + dimCount - 2, dimCount - 1);
	for (int dim = 0; dim < dimCount; dim++) {
		std::fill(x.begin(), x.end(), Type(0));
		int oi = 0;
		int const cdim = dimCount - 1 - dim;
		for (int hi = 0; hi < blockCount; hi++) {
			int loMax = loSize;
			for (int mi = 0; mi < modeCount; mi++) {
				printf("%i %i %i\n", dim, loMax, loSize);
				for (int lo = 0; lo < loMax; lo++) {
					for (int ni = 0; ni < nodeCount; ni++) {
						int const xIndex = (hi * nodeCount + ni) * loSize + lo;
						x[xIndex] += lT[ni][mi] * y[oi];
					}
					oi++;
				}
				if (mi + 1 < modeCount) {
					loMax *= modeCount - mi - 1;
					loMax /= modeCount + cdim - mi - 1;
				}
			}
		}
		/** (modeCount + dimCount - 2, dimCount - 1)*/
		if (cdim) {
			loSize *= cdim;
			loSize /= modeCount + cdim - 1;
		}
		blockCount *= nodeCount;
		std::swap(x, y);
	}
	return y;
}
