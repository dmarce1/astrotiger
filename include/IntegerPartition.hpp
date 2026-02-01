/*
 * IntegerPartition.hpp
 *
 *  Created on: Dec 28, 2025
 *      Author: dmarce1
 */

#pragma once

#include "Math.hpp"

#include <array>
#include <concepts>
#include <tuple>
#include <utility>

template <int... λ>
requires (nonIncreasing<λ...>())
class IntegerPartition {
	static constexpr auto conjParts() {
		constexpr int conjSize = values[0];
		std::array<int, conjSize> conjPart;
		int n = count();
		for (int col = 0; col < conjSize; col++) {
			while (col + 1 > values[n - 1]) {
				n--;
			}
			conjPart[col] = n;
		}
		return conjPart;
	}
	static constexpr std::array<int, sizeof...(λ)> values = { λ... };

public:
	static constexpr int count() {
		return sizeof...(λ);
	}
	static constexpr int size() {
		return (0 + ... + λ);
	}
	static constexpr auto conj() {
		constexpr int conjSize = values[0];
		auto const create = []<int... δ>(std::integer_sequence<int, δ...>) {
			constexpr auto conjPart = conjParts();
			return IntegerPartition<conjPart[δ]...> {};
		};
		constexpr auto seq = std::make_integer_sequence<int, conjSize> { };
		return create(seq);
	}
	constexpr int operator[](int i) const {
		return values[i];
	}
	constexpr auto begin() const {
		return values.cbegin();
	}
	constexpr auto end() const {
		return values.cend();
	}
};

template<typename T>
struct IsIntegerPartition {
	static constexpr bool value = false;
};

template<int ... λ>
struct IsIntegerPartition<IntegerPartition<λ...>> {
	static constexpr bool value = true;
};

namespace detail {
template<typename ... Partitions>
struct PartitionList {
	static constexpr int size = sizeof...(Partitions);
	template<int I>
	using at = std::tuple_element_t<I, std::tuple<Partitions...>>;

	template<typename F>
	static constexpr void forEachType(F &&f) {
		(static_cast<void>(f.template operator()<Partitions>()), ...);
	}
};
}

template<int ... λ>
constexpr auto conj(IntegerPartition<λ...> const &part);

template<typename ... Partitions>
using IntegerPartitionList = std::enable_if_t<(IsIntegerPartition<Partitions>::value && ...), detail::PartitionList<Partitions...>>;

template<typename T>
struct IsIntegerPartitionList {
	static constexpr bool value = false;
};

template <typename... Partitions>
struct IsIntegerPartitionList<detail::PartitionList<Partitions...>> {
  static constexpr bool value = true;
};


template <typename T>
concept IntegerPartitionType =
  IsIntegerPartition<std::remove_cvref_t<T>>::value;

template <typename T>
concept IntegerPartitionListType =
  IsIntegerPartitionList<std::remove_cvref_t<T>>::value;
