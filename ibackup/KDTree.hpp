/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "Definitions.hpp"
#include "Box.hpp"
#include "FixedPrecision.hpp"

#include <memory>

template<size_t>
struct KDTree;

template<size_t dimCount>
struct KDTree {
	using FixedType = FixedPrecision<32>;
	using BoxType = Box<FixedType, dimCount>;
	using Pointer = std::shared_ptr<KDTree>;
private:
	static constexpr size_t childCount = 2;
	struct Private {
		explicit Private() = default;
	};
	BoxType box_;
	Pointer left_;
	Pointer right_;
	uint64_t key_;
	template<typename Arc>
	void serialize(Arc &&arc, unsigned) {
		arc & box_;
		arc & left_;
		arc & right_;
		arc & key_;
	}
public:
	KDTree(Private) :
			box_(BoxType::cube(FixedType::min(), FixedType::max())), left_(nullptr), right_(nullptr), key_(1) {
	}
	static Pointer create() {
		return std::make_shared<KDTree>(Private());
	}
};

