/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
*******************************************************************************/

#pragma once

#include "Definitions.hpp"

#include <memory>

template<Dimension dimCount>
using LeafPointer = std::shared_ptr<Leaf<dimCount>>;


template<Dimension dimCount>
struct Leaf {
	using pointer_t = LeafPointer<dimCount>;
	Leaf& operator=(Leaf const&) = default;
	Leaf& operator=(Leaf&&) = default;
	static pointer_t create() {
		return Leaf{};
	}
private:
	Leaf() = default;
};

