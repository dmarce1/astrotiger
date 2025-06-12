#pragma once

#include <cstdint>

inline constexpr uint64_t binomialCoefficient(uint64_t n, uint64_t k) {
	uint64_t num = 1;
	uint64_t den = 1;
	for (uint64_t i = 1; i <= k; i++) {
		num *= n + 1 - i;
		den *= i;
	}
	return num / den;
}

template<uint64_t dimCount>
inline constexpr uint64_t triangularNumber(uint64_t size) {
	return binomialCoefficient(size + dimCount - 1, dimCount);
}
