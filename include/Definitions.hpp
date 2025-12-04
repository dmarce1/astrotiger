/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include <cstddef>
#include <cstdint>

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
