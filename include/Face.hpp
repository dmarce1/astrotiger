/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "Definitions.hpp"

template<Dimension dimCount>
struct Face {
	using Type = signed char;

	constexpr Face(Face const&);
	constexpr Face(Face&&);
	constexpr Face& operator=(Face const&);
	constexpr Face& operator=(Face&&);
	template<typename Arc>
	void serialize(Arc &&arc, unsigned) {
		arc & bits_;
	}
	constexpr Face() :
			bits_(dimCount << 1) {
	}
	constexpr Face(Dimension dim, Sign dir) {
		bits_ = (dim << 1) + ((dir + 1) >> 1);
	}
	constexpr void setDimension(Dimension dim) {
		bits_ = (dim << 1) | (bits_ & 1);
	}
	constexpr void setDirection(Sign dir) {
		bits_ = ((dir + 1) >> 1) | (bits_ & ~1);
	}
	constexpr Dimension getDimension() const {
		return bits_ >> 1;
	}
	constexpr Sign getDirection() const {
		return 2 * (bits_ & 1) - 1;
	}
	constexpr bool valid() const {
		return bool((0 <= bits_) && (bits_ < (dimCount << 1)));
	}
	constexpr explicit operator Index() const {
		return bits_;
	}
	constexpr Face& operator++() {
		bits_++;
		return *this;
	}
	constexpr Face operator++(int) const {
		Face original = *this;
		operator++();
		return original;
	}
	constexpr Face flip() const {
		return Face(bits_ ^ 1);
	}
	static constexpr Face begin() {
		return Face(0);
	}
	static constexpr Face end() {
		return Face(dimCount << 1);
	}
	static constexpr size_t count() {
		return dimCount << 1;
	}
private:
	Type bits_;
	constexpr Face(Type b) :
			bits_(b) {
	}
};

