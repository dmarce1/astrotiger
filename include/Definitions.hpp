/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <ranges>
#include <type_traits>

using Dimension = uint_fast8_t;
using Rank = uint_fast8_t;
using Index = uint_fast64_t;
using Sign = int_fast8_t;

template<typename, Dimension>
struct Box;

template<Dimension>
struct Child;

template<Dimension>
struct Face;

template<int, typename = double>
struct FixedPrecision;

template<Dimension>
struct Leaf;

template<Dimension>
struct MortonKey;

template<Dimension>
struct Octree;

template<typename, Dimension>
struct Point;

template<Dimension dimCount>
constexpr size_t faceCount = (dimCount << 1);

template<Dimension dimCount>
constexpr size_t childCount = (1 << dimCount);

template<class T>
concept Container = std::ranges::range < T > &&std::default_initializable < T > &&std::move_constructible < T > &&std::assignable_from<T&, T>;

template <typename T, typename I>
concept Indexed = std::integral<I> && requires(T t, I i) {t[i];};

template <std::unsigned_integral auto I>
using int_least = std::conditional_t<
I <= std::numeric_limits<uint_least8_t>::max(), uint_least8_t,
std::conditional_t<I <= std::numeric_limits<uint_least16_t>::max(), uint_least16_t,
std::conditional_t<I <= std::numeric_limits<uint_least32_t>::max(), uint_least32_t,
std::conditional_t<I <= std::numeric_limits<uint_least64_t>::max(), uint_least64_t, void>>>>;

template <std::unsigned_integral auto I>
using int_fast = std::conditional_t<
I <= std::numeric_limits<uint_fast8_t>::max(), uint_fast8_t,
std::conditional_t<I <= std::numeric_limits<uint_fast16_t>::max(), uint_fast16_t,
std::conditional_t<I <= std::numeric_limits<uint_fast32_t>::max(), uint_fast32_t,
std::conditional_t<I <= std::numeric_limits<uint_fast64_t>::max(), uint_fast64_t, void>>>>;

template<typename T>
struct ConvertsToLiteral {
	static constexpr bool value = false;
};

constexpr void constexprAssert(bool cond) {
	if (std::is_constant_evaluated()) {
		if (!cond) {
			throw 0;
		}
	} else {
		assert(cond);
	}
}

template<typename T>
concept Number = std::integral<T> || std::floating_point<T>;

template<Number T>
constexpr int largest = std::numeric_limits<T>::max();

template<Number T>
constexpr int smallest = std::is_integral<T>::value ? std::numeric_limits<T>::min() : -largest<T>;

template<std::floating_point T>
constexpr int tiny = std::numeric_limits<T>::min();

#define ASSERT(b)                                 \
	if(!(b)) {                                    \
		std::string estr = "Exception throw in "; \
		estr += __FILE__;                         \
		estr += " on line ";                      \
		estr += std::to_string(__LINE__);         \
		estr += ".\n";                            \
		throw std::runtime_error(estr);           \
	}

#define NUMERICAL_CONSTANTS(T)          \
		constexpr static auto zero = static_cast<T>(0); \
		constexpr static auto one = static_cast<T>(1);

