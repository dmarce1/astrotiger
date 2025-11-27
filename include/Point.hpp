/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "Definitions.hpp"
#include "FixedPrecision.hpp"

#include <array>

template<typename Type, size_t dimCount>
struct Point {
	static constexpr auto zero = Type(0.0);
	template<typename Arc>
	void serialize(Arc &&arc, unsigned) {
		arc & coords_;
	}
	constexpr Point() = default;
	constexpr Point(Point const&) = default;
	constexpr Point(Point&&) = default;
	constexpr Point& operator=(Point const&) = default;
	constexpr Point& operator=(Point&&) = default;
	constexpr Type operator[](size_t i) const {
		return coords_[i];
	}
	constexpr Type& operator[](size_t i) {
		return coords_[i];
	}
	constexpr Point operator-(Point const &other) const {
		Point dif;
		for(Dimension dim = 0; dim < dimCount; dim++){
			dif.coords_[dim] = coords_[dim] - other.coords_[dim];
		}
		return dif;
	}
	constexpr bool operator==(Point const &other) const {
		return coords_ == other.coords_;
	}
	constexpr bool operator!=(Point const &other) const {
		return coords_ != other.coords_;
	}
private:
	constexpr Point(std::array<Type, dimCount> const &coords) :
			coords_(coords) {
	}
	std::array<Type, dimCount> coords_;
};

