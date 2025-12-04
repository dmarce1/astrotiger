/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "Definitions.hpp"

#include <array>

#include "Face.hpp"

template<Dimension dimCount>
struct Child {
	using Type = unsigned char;
	using Face_t = Face<dimCount>;
	constexpr Child(Child const&);
	constexpr Child(Child&&);
	constexpr Child& operator=(Child const&);
	constexpr Child& operator=(Child&&);
	template<typename Arc>
	void serialize(Arc &&arc, unsigned) {
		arc & bits_;
	}
	constexpr Child() :
			bits_(1 << dimCount) {
	}
	constexpr Child(std::array<Sign, dimCount> const &signs) :
			bits_(0) {
		set(signs);
	}
	constexpr void set(std::array<Sign, dimCount> const &signs) {
		for (Index dim = 0; dim < dimCount; dim++) {
			set(dim, signs[dim]);
		}
	}
	constexpr void set(Dimension dim, Sign dir) {
		bits_ &= ~(1 << (maxDim - dim));
		bits_ |= (dir + 1) >> 1;
	}
	constexpr Dimension get(Dimension dim) const {
		return bits_ >> (maxDim - dim);
	}
	constexpr std::array<Sign, dimCount> get() const {
		std::array<Sign, dimCount> signs;
		for (Index dim = 0; dim < dimCount; dim++) {
			signs[dim] = get(dim);
		}
		return signs;
	}
	constexpr bool valid() const {
		return bool((0 <= bits_) && (bits_ < (1 << dimCount)));
	}
	constexpr bool shares(Face_t const &Face) const {
		Type mask = Type(1) << Face.getDimension();
		if (Face.getDirection() < 0) {
			return bool(bits_ ^ mask);
		} else {
			return bool(bits_ & mask);
		}
	}
	constexpr explicit operator Index() const {
		return bits_;
	}
	constexpr Child& operator++() {
		bits_++;
		return *this;
	}
	constexpr Child operator++(int) const {
		Child original = *this;
		operator++();
		return original;
	}
	constexpr Child flip() const {
		return Child(bits_ ^ bitMask);
	}
	static constexpr Child begin() {
		return Child(0);
	}
	static constexpr Child end() {
		return Child(1 << dimCount);
	}
	static constexpr size_t count() {
		return 1 << dimCount;
	}
private:
	static constexpr Dimension maxDim = dimCount - 1;
	static constexpr Type bitMask = (1 << dimCount) - 1;
	Type bits_;
	constexpr Child(Type b) :
			bits_(b) {
	}
};

