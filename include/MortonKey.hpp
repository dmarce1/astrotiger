/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "Definitions.hpp"

#include <bit>

#include "Child.hpp"
#include "Face.hpp"

template<Dimension dimCount>
struct MortonKey {
	using Integer = std::uint64_t;
	static constexpr Index maxLevel = sizeof(Integer) * CHAR_BIT - 1;
	constexpr MortonKey(MortonKey const&);
	constexpr MortonKey(MortonKey&&);
	constexpr MortonKey& operator=(MortonKey const&);
	constexpr MortonKey& operator=(MortonKey&&);
	constexpr MortonKey() :
			bits_(1) {
	}
	constexpr MortonKey genChild(Child<dimCount> ci) const {
		return MortonKey((bits_ << dimCount) | Index(ci));
	}
	constexpr MortonKey genFace(Face<dimCount> si) const {
		auto const idim = maxDim - si.getDimension();
		auto const mask = getBitMask() << idim;
		auto const one = 1 << idim;
		if (si.getDirection() > 0) {
			return MortonKey((one + (bits_ | ~mask)) | (bits_ & ~mask));
		} else {
			return MortonKey((-one + (bits_ | mask)) | (bits_ & ~mask));
		}
	}
	constexpr MortonKey genParent() const {
		return MortonKey(bits_ >> dimCount);
	}
	constexpr bool level() const {
		return std::bit_width(bits_);
	}
	constexpr bool valid() const {
		return bool(bits_ != 0);
	}
	template<typename Arc>
	void serialize(Arc &&arc, unsigned) {
		arc & bits_;
	}
private:
	static constexpr Dimension maxDim = dimCount - 1;
	Integer bits_;
	constexpr auto getBitMask(Integer mask = Integer(1)) {
		if (mask) {
			return 1 | getBitMask(mask << dimCount);
		} else {
			return 0;
		}
	}
	constexpr MortonKey(Integer b) :
			bits_(b) {
	}
};
