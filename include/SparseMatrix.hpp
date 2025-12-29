
#pragma once

#include <concepts>
#include <tuple>
#include <type_traits>

template <int N, int M, auto X>
struct MatrixElement {
	using value_type = decltype(X);
	static constexpr int row = N;
	static constexpr int column = M;
	static constexpr value_type value = X;
};

template <int N, int M, int J, int K>
struct MatrixElementCompare {
	static constexpr bool eq = (N == J) && (M == K);
	static constexpr bool lt = (N == J) ? (M < K) : (N < J);
	static constexpr bool gt = !(lt || eq);
	static constexpr bool neq = !eq;
	static constexpr bool lte = !gt;
	static constexpr bool gte = !lt;
};

template <typename>
struct IsMatrixElement {
	static constexpr bool value = false;
};

template <int N, int M, auto X>
struct IsMatrixElement<MatrixElement<N, M, X>> {
	static constexpr bool value = true;
};

template <typename>
struct IsSparseMatrix {
	static constexpr bool value = false;
};

template <typename... Args>
struct IsSparseMatrix<std::tuple<Args...>> {
	static constexpr bool value = (true && ... && IsMatrixElement<Args>::value);
};

template <typename T>
concept SparseMatrix = IsSparseMatrix<std::remove_cvref_t<T>>::value;

//
// lo = 0
// hi = n - 1
//
// while lo <= hi:
//    mid = lo + (hi - lo) // 2   # avoids overflow
//
//    if a[mid] == x:
//        return mid
//    else if a[mid] < x:
//        lo = mid + 1
//    else:
//        hi = mid - 1
//
// return -1

template <SparseMatrix Tuple, int N, int M, int LO = 0, int HI = std::tuple_size<Tuple>::value - 1>
class FindMatrixElement {
	static constexpr int I = (HI + LO) / 2;
	static constexpr int J = std::tuple_element_t<I, Tuple>::row;
	static constexpr int K = std::tuple_element_t<I, Tuple>::column;
	using Compare = MatrixElementCompare<J, K, N, M>;

public:
	static constexpr int value =
		((((HI <= LO) && Compare::neq) ? -1
									   : (Compare::eq ? I
													  : (Compare::lt ? FindMatrixElement<Tuple, N, M, I + 1, HI>::value
																	 : FindMatrixElement<Tuple, N, M, LO, I - 1>::value))));
};
