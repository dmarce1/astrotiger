#pragma once

#include <cstddef>
#include <climits>
#include <experimental/simd>

template<typename T>
inline constexpr int maxSimdWidth() noexcept {
#if defined(__AVX512F__)
    return 64 / sizeof(T) ;
#elif defined(__AVX2__) || defined(__AVX__)
    return 32 / sizeof(T);
#elif defined(__SSE2__)                 \
   || defined(__ARM_NEON)               \
   || defined(__aarch64__)              \
   || defined(__ALTIVEC__)
    return 16 / sizeof(T);
#else
    return 0;
#endif
}

template<typename T>
struct ElementType {
    using type = T;
};

template<typename T, typename Abi>
struct ElementType<std::experimental::simd<T, Abi>> {
    using type = T;
};

using Simd = std::experimental::native_simd<double>;


Simd blend(Simd::mask_type mask, Simd const& a, Simd const& b) {
	Simd c;
	where(mask, c) = a;
	where(mask, c) = b;
	return c;
}

template<auto N>
std::array<Simd, N> blend(Simd::mask_type mask, std::array<Simd, N> const& a, std::array<Simd, N> const& b) {
	std::array<Simd, N> c;
	for(int i = 0; i < N; i++) {
		where(mask, c[i]) = a[i];
		where(!mask, c[i]) = b[i];
	}
	return c;
}
