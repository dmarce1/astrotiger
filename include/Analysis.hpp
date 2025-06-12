/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_ANALYSIS_HPP_
#define INCLUDE_ANALYSIS_HPP_

#include "Combinatorics.hpp"
#include "Util.hpp"

#include <array>

namespace detail {

template<typename T>
constexpr T legendre(int N, T x) {
	T Pnp1, Pnm1;
	T Pn = T(1);
	if (N >= 1) {
		Pnp1 = x;
		Pnm1 = Pn;
		Pn = Pnp1;
		for (int n = 1; n < N; n++) {
			Pnp1 = (T(2 * n + 1) * x * Pn - T(n) * Pnm1) / T(n + 1);
			Pnm1 = Pn;
			Pn = Pnp1;
		}
	}
	return Pn;
}
}

template<typename T, int N>
struct QuadraturePoints {
	std::array<T, N> x;
	std::array<T, N> w;
};

template<typename T, int N>
constexpr auto getLegendreQuadraturePoints() {
	using std::cos;
	std::array<T, N> x;
	std::array<T, N> w;
	T const Np1 = T(N + 1);
	T const N_ = T(N);
	for (int j = 0; j < N; j++) {
		x[j] = cos(std::numbers::pi_v < T > *(T(1) - (T(j) + (T(1) / T(2))) / N_));
		T xPrevious, Pnp1;
		do {
			xPrevious = x[j];
			Pnp1 = detail::legendre(N + 1, x[j]);
			T const Pn = detail::legendre(N, x[j]);
			T const x2 = sqr(x[j]);
			T const num = (T(1) - x2) * Pn;
			T const den = Np1 * (Pnp1 - x[j] * Pn);
			x[j] += num / den;
		} while (T(x[j]) != T(xPrevious));
		T const dpdx = Np1 * (detail::legendre(N + 1, x[j]) - x[j] * detail::legendre(N, x[j])) / (sqr(x[j]) - T(1));
		w[j] = T(2) / ((T(1) - sqr(x[j])) * sqr(dpdx));
	}
	return QuadraturePoints<T, N>( { x, w });
}

namespace detail {

template<typename Type, int modeCount, int nodeCount>
constexpr auto legendreAnalysis() {
	constexpr auto quad = getLegendreQuadraturePoints<Type, nodeCount>();
	constexpr auto qPoint = quad.x;
	constexpr auto qWeight = quad.w;
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
	constexpr auto qPoint = quad.x;
	std::array<std::array<Type, modeCount>, nodeCount> transform;
	for (int mi = 0; mi < modeCount; mi++) {
		for (int ni = 0; ni < nodeCount; ni++) {
			transform[ni][mi] = legendre<Type>(mi, qPoint[ni]);
		}
	}
	return transform;
}

template<typename Type, int Cols, int J = 0>
inline constexpr Type computeRowSum(std::array<Type, Cols> const &row, std::array<Type, Cols> const &x) {
	if constexpr (J >= Cols) {
		return Type(0);
	} else {
		Type const term = (row[J] == Type(0)) ? Type(0) : row[J] * x[J];
		return term + computeRowSum<Type, Cols, J + 1>(row, x);
	}
}

template<typename Type, int Rows, int Cols>
inline constexpr auto multiplySparseMatrixVector(std::array<std::array<Type, Cols>, Rows> const &matrix, std::array<Type, Cols> const &x) {
	std::array<Type, Rows> result;
	for (int i = 0; i < Rows; ++i) {
		result[i] = computeRowSum<Type, Cols>(matrix[i], x);
	}
	return result;
}
}

template<typename Type, int dimCount, int modeCount, int nodeCount = modeCount, int dimension = dimCount - 1>
auto modalAnalysis(std::array<Type, power(nodeCount, dimension + 1) * binco(modeCount + dimCount - dimension - 2, dimCount - dimension - 1)> const &x) {
	constexpr int loSize = binco(modeCount + dimCount - dimension - 2, dimCount - dimension - 1);
	constexpr int inSize = loSize * power(nodeCount, dimension + 1);
	constexpr int outSize = loSize * power(nodeCount, dimension) * (modeCount + dimCount - 1 - dimension) / (1 + dimCount - 1 - dimension);
	constexpr auto T = []() {
		constexpr auto lT = detail::legendreAnalysis<Type, modeCount, nodeCount>();
		constexpr int blockCount = power(nodeCount, dimension);
		std::array<std::array<Type, inSize>, outSize> T;
		for (auto &row : T) {
			std::fill(row.begin(), row.end(), Type(0));
		}
		int oi = 0;
		for (int hi = 0; hi < blockCount; hi++) {
			int loMax = loSize;
			for (int mi = 0; mi < modeCount; mi++) {
				for (int lo = 0; lo < loMax; lo++) {
					for (int ni = 0; ni < nodeCount; ni++) {
						int const xIndex = (hi * nodeCount + ni) * loSize + lo;
						T[oi][xIndex] = lT[mi][ni];
					}
					oi++;
				}
				if (mi + 1 < modeCount) {
					loMax *= modeCount - mi - 1;
					loMax /= modeCount + dimCount - 1 - dimension - mi - 1;
				}
			}
		}
		return T;
	}();
	auto y = detail::multiplySparseMatrixVector<Type, outSize, inSize>(T, x);
	if constexpr (dimension > 0) {
		return modalAnalysis<Type, dimCount, modeCount, nodeCount, dimension - 1>(y);
	} else {
		return y;
	}
}

template<typename Type, int dimCount, int modeCount, int nodeCount = modeCount, int dimension = 0>
auto modalSynthesis(std::array<Type, binco(modeCount + dimCount - dimension - 1, dimCount - dimension) * power(nodeCount, dimension)> const &x) {
	constexpr int loSize = binco(modeCount + dimCount - dimension - 2, dimCount - dimension - 1);
	constexpr int outSize = loSize * power(nodeCount, dimension + 1);
	constexpr int inSize = (modeCount + dimCount - dimension - 1) * loSize * power(nodeCount, dimension) / (dimCount - dimension);
	constexpr auto T = []() {
		auto lT = detail::legendreSynthesis<Type, modeCount, nodeCount>();
		constexpr int blockCount = power(nodeCount, dimension);
		std::array<std::array<Type, inSize>, outSize> T { };
		for (auto &row : T) {
			std::fill(row.begin(), row.end(), Type(0));
		}
		int oi = 0;
		for (int hi = 0; hi < blockCount; hi++) {
			int loMax = loSize;
			for (int mi = 0; mi < modeCount; mi++) {
				for (int lo = 0; lo < loMax; lo++) {
					for (int ni = 0; ni < nodeCount; ni++) {
						int const xIndex = (hi * nodeCount + ni) * loSize + lo;
						T[xIndex][oi] = lT[ni][mi];
					}
					oi++;
				}
				if (mi + 1 < modeCount) {
					loMax *= modeCount - mi - 1;
					loMax /= modeCount + dimCount - dimension - mi - 2;
				}
			}
		}
		return T;
	}();
	auto y = detail::multiplySparseMatrixVector<Type, outSize, inSize>(T, x);
	if constexpr (dimension != dimCount - 1) {
		return modalSynthesis<Type, dimCount, modeCount, nodeCount, dimension + 1>(y);
	} else {
		return y;
	}
}

#endif /* INCLUDE_ANALYSIS_HPP_ */
