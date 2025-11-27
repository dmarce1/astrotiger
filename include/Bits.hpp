/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include <climits>

template<typename Int>
constexpr int msb(Int bits) {
	if (bits) {
		int bit = 0;
		for (int test = (sizeof(Int) * CHAR_BIT) >> 1; test; test >>= 1) {
			int testLev = bit | test;
			if (bits & ~Int((Int(1) << testLev) - 1)) {
				bit = testLev;
			}
		}
		return bit;
	} else {
		return -1;
	}
}

