/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "Definitions.hpp"

#include <climits>
#include <type_traits>

template<int width>
using SignedInteger = std::conditional_t<
	(width <= sizeof(signed char) * CHAR_BIT),
	char,
	std::conditional_t< (width <= sizeof(signed short) * CHAR_BIT),
		short,
		std::conditional_t< (width <= sizeof(signed long) * CHAR_BIT),
			long,
			std::conditional_t< (width <= sizeof(signed long long) * CHAR_BIT),
				long long,
				void
			>
		>
	>
>;

template<int width>
using UnsignedInteger = std::make_signed_t<SignedInteger<width>>;

