/*
 * IntegerPartition.hpp
 *
 *  Created on: Dec 28, 2025
 *      Author: dmarce1
 */

#pragma once

#include "Math.hpp"

#include <array>
#include <utility>

template <int... λ>
	requires(nonIncreasing<λ...>())
class IntegerPartition {
	static constexpr auto conjParts() {
		constexpr int conjSize = values[0];
		std::array<int, conjSize> conjPart;
		int n = size();
		for (int col = 0; col < conjSize; col++) {
			while (col + 1 > values[n - 1]) {
				n--;
			}
			conjPart[col] = n;
		}
		return conjPart;
	}
	static constexpr std::array<int, sizeof...(λ)> values = {λ...};

public:
	static constexpr int size() {
		return sizeof...(λ);
	}
	static constexpr auto conj() {
		constexpr int conjSize = values[0];
		auto const create = []<int... δ>(std::integer_sequence<int, δ...>) {
			constexpr auto conjPart = conjParts();
			return IntegerPartition<conjPart[δ]...>{};
		};
		constexpr auto seq = std::make_integer_sequence<int, conjSize>{};
		return create(seq);
	}
	constexpr int operator[](int i) const {
		return values[i];
	}
};

template <int... λ>
constexpr auto conj(IntegerPartition<λ...> const &part);
