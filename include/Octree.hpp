/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "Definitions.hpp"

#include <memory>

#include "Box.hpp"
#include "FixedPrecision.hpp"
#include "Leaf.hpp"
#include "MortonKey.hpp"

template<Dimension dimCount>
using OctreePointer = std::shared_ptr<Octree<dimCount>>;

template<Dimension dimCount>
struct Octree {
	using box_t = Box<FixedPrecision<32, double>, dimCount>;
	using child_t = Child<dimCount>;
	using fixed_t = FixedPrecision<32, double>;
	using key_t = MortonKey<dimCount>;
	using pointer_t = OctreePointer<dimCount>;
	static constexpr size_t nchildren = childCount<dimCount>;
	Octree& operator=(Octree const&) = default;
	Octree& operator=(Octree&&) = default;
	void refine() {
		auto const cdomains = domain_.split();
		for (auto ci = child_t::begin(); ci != child_t::end(); ci++) {
			children_[ci] = Octree::create(cdomains[ci], key_.getChild(ci));
		}
	}
	static pointer_t create() {
		box_t box = box_t::cube(fixed_t::min(), fixed_t::max());
		key_t key { };
		return create(box, key);
	}
	static pointer_t create(box_t const &domain, key_t key) {
		pointer_t ptr = std::make_shared<Octree>();
		ptr->domain_ = domain;
		ptr->key_ = key;
	}
private:
	Octree() = default;
	std::array<pointer_t, nchildren> children_;
	box_t domain_;
	key_t key_;
};

