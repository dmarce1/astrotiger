/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "Definitions.hpp"
#include "Point.hpp"

#include <climits>
#include <cmath>
#include <numeric>

template<typename Type, size_t dimCount>
struct Box {
	static constexpr auto zero = Type(0.0);
	template<typename Arc>
	void serialize(Arc &&arc, unsigned) {
		arc & a_;
		arc & b_;
	}
	constexpr Box() {
		*this = null();
	}
	constexpr Box(Box const&) = default;
	constexpr Box(Box&&) = default;
	constexpr Box& operator=(Box const&) = default;
	constexpr Box& operator=(Box&&) = default;
	constexpr bool operator==(Box const &other) const {
		return (a_ == other.a_) && (b_ == other.b_);
	}
	constexpr bool operator!=(Box const &other) const {
		return (a_ != other.a_) || (b_ != other.b_);
	}
	constexpr std::pair<Box, Box> split(Dimension dim) const {
		using std::midpoint;
		return split(dim, midpoint(a_[dim], b_[dim]));
	}
	constexpr std::pair<Box, Box> split(Dimension dim, Type pos) const {
		using std::numeric_limits;
		std::pair<Box, Box> children(*this, *this);
		children.first.b_[dim] = pos;
		children.second.a_[dim] = nexttoward(pos, numeric_limits<Type>::max());
		return children;
	}
	constexpr bool contains(Point<Type, dimCount> const &point) const {
		for (Dimension k = 0; k < dimCount; k++) {
			if (point[k] < a_[k]) {
				return false;
			}
			if (point[k] > b_[k]) {
				return false;
			}
		}
		return true;
	}
	constexpr bool contains(Box const &box) const {
		return contains(box.a_) && contains(box.b_);
	}
	constexpr bool intersects(Box const &other) const {
		for (Dimension k = 0; k < dimCount; k++) {
			auto const a = max(a_[k], other.a_[k]);
			auto const b = min(b_[k], other.b_[k]);
			if (a >= b) {
				return false;
			}
		}
		return true;
	}
	constexpr Box intersection(Box const &other) const {
		Box iBox;
		for (Dimension k = 0; k < dimCount; k++) {
			iBox.a_[k] = max(a_[k], other.a_[k]);
			iBox.b_[k] = min(b_[k], other.b_[k]);
			if (iBox.a_[k] >= iBox.b_[k]) {
				return null();
			}
		}
		return iBox;
	}
	constexpr Type span(Dimension dim) const {
		return Real(a_[dim] - b_[dim]);
	}
	constexpr Type volume() const {
		Type vol = Type(1);
		for (Dimension k = 0; k < dimCount; k++) {
			if (a_[k] >= b_[k]) {
				return Type(0);
			}
			vol *= span(k);
		}
		return vol;
	}
	friend constexpr Box bounding(Box const &boxA, Box const &boxB) {
		using std::min;
		using std::max;
		Point<Type, dimCount> a, b;
		for (Dimension k = 0; k < dimCount; k++) {
			a[k] = min(boxA.a_[k], boxB.a_[k]);
			b[k] = max(boxA.b_[k], boxB.b_[k]);
		}
		return Box(a, b);
	}
	static constexpr Box null() {
		using std::min;
		Point<Type, dimCount> a, b;
		for (Dimension k = 0; k < dimCount; k++) {
			b[k] = a[k] = Type::min();
		}
		return Box(a, b);
	}
	static constexpr Box cube(Type minValue, Type maxValue) {
		Point<Type, dimCount> a, b;
		for (Dimension k = 0; k < dimCount; k++) {
			a[k] = minValue;
			b[k] = maxValue;
		}
		return Box(a, b);
	}
private:
	constexpr Box(Point<Type, dimCount> const &a, Point<Type, dimCount> const &b) :
			a_(a), b_(b) {
	}
	Point<Type, dimCount> a_;
	Point<Type, dimCount> b_;
};

