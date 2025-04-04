#pragma once

#include <array>
#include <bitset>
#include <optional>
#include <type_traits>
#include <utility>

#include "IndexTuple.hpp"
#include "Numbers.hpp"
#include "Permutation.hpp"
#include "SymmetricGroup.hpp"

namespace Tensors {

template<char C>
struct Index {
	static constexpr char value = C;
	constexpr operator char() const {
		return C;
	}
};

template<typename T>
struct IsIndex {
	static constexpr bool value = false;
};

template<char C>
struct IsIndex<Index<C>> {
	static constexpr bool value = true;
	operator char() const {
		return value;
	}
};

struct map_index_t;

template<typename T, typename F, size_t D, size_t R, Symmetries<R> S, char ...I>
struct TensorExpr;

template<typename T, size_t D, size_t R, Symmetries<R> S>
struct Tensor;

struct map_index_t {
	signed char sign;
	unsigned index;
};

template<typename U, size_t ... I>
auto getTupleTail(const U &tup, std::index_sequence<I...>) {
	return std::make_tuple(std::get<I + 1>(tup)...);
}

template<typename ... Args>
auto tupleTail(const std::tuple<Args...> &tup) {
	return getTupleTail(tup, std::make_index_sequence<sizeof...(Args) - 1> { });
}

template<typename ... Args>
std::tuple<> createIndexTuple(std::tuple<Args...> tup) {
	return std::tuple<Args...>();
}

template<char C, typename ...Args>
auto createIndexTuple(std::tuple<Index<C>, Args...> tup, size_t i, auto ...J) {
	std::tuple<size_t> const thisTuple(i);
	auto const nextTuple = createIndexTuple(tuple_tail(tup), J...);
	return std::tuple_cat(thisTuple, nextTuple);
}

template<typename ...Args>
auto createIndexTuple(std::tuple<int, Args...> tup, auto ...is) {
	std::tuple<size_t> const thisTuple(std::get<0>(tup));
	auto const nextTuple = createIndexTuple(tuple_tail(tup), is...);
	return std::tuple_cat(thisTuple, nextTuple);
}

template<typename ...Args>
std::tuple<> extractChars(Args...) {
	return std::tuple<>();
}

template<typename ...Args>
auto extractChars(int i, Args ...args) {
	return extractChars(args...);
}

template<char C, typename ...Args>
auto extractChars(Index<C>, Args ...args) {
	return std::tuple_cat(std::tuple<Index<C>>(), extractChars(args...));
}

template<typename T, size_t R, Permutation<R> P>
constexpr auto permuteList(T const &list) {
	auto permutedList = list;
	for (size_t i = 0; i < list.size(); i++) {
		for (size_t j = 0; j < R; j++) {
			permutedList[i][j] = list[i][P[j] - 1];
		}
	}
	return permutedList;
}

template<size_t N, std::array<char, N> A, std::array<char, N> B>
constexpr Permutation<N> computePermutation() {
	Permutation<N> P;
	for (size_t i = 0; i < N; i++) {
		size_t j = 0;
		while (A[i] != B[j]) {
			j++;
		}
		static_assert(j != N);
		P[i] = j + 1;
	}
}

template<size_t D, size_t R, Symmetries<R> S>
constexpr auto computeUniqueIndices() {
	using i_type = IndexTuple<D, R>;
	static constexpr auto symmetrySpecification = genSymmetryPermutations(S);
	static constexpr size_t Size = Math::integerPower(D, R);
	std::array<i_type, Size> uniqeIndices;
	std::bitset<Size> haveVisited;
	size_t tensorIndex = 0;
	for (i_type tensorIndices = i_type::begin(); tensorIndices != i_type::end(); tensorIndices++) {
		if (!haveVisited[tensorIndex]) {
			haveVisited[tensorIndex] = true;
			uniqeIndices[tensorIndex++] = tensorIndices;
			for (size_t n = 0; n < symmetrySpecification.size(); n++) {
				auto const permutation = symmetrySpecification[n].second;
				auto permutedIndices = permutation.apply(tensorIndices);
				do {
					haveVisited[permutedIndices.flatIndex()] = true;
					permutedIndices = permutation.apply(permutedIndices);
				} while (permutedIndices != tensorIndices);
			}
		}
	}
	return uniqeIndices;

}

template<size_t D, size_t R, size_t N, Symmetries<R> S>
constexpr size_t computeUniqueIndicesSize() {
	using i_type = IndexTuple<D, R>;
	static constexpr size_t Size = Math::integerPower(D, R);
	std::bitset<Size> haveVisited;
	size_t count = 0;
	for (i_type tensorIndices = i_type::begin(); tensorIndices != i_type::end(); tensorIndices++) {
		if (!haveVisited[count]) {
			haveVisited[count] = true;
			for (size_t n = 0; n < N; n++) {
				auto const permutation = S[n].second;
				auto permutedIndices = permutation.apply(tensorIndices);
				do {
					haveVisited[permutedIndices.flatIndex()] = true;
					permutedIndices = permutation.apply(permutedIndices);
				} while (permutedIndices != tensorIndices);
			}
		}
	}
	return count;
}

template<size_t D, size_t R, Symmetries<R> S>
constexpr auto createIndexMap() {
	using i_type = IndexTuple<D, R>;
	static constexpr size_t Size = Math::integerPower(D, R);
	static constexpr auto permutations = genSymmetryPermutations<R, S>();
	static constexpr size_t N = permutations.size();
	std::array<map_index_t, Size> indexMap;
	std::bitset<Size> haveVisited;
	size_t elementCount = 0;
	for (i_type tensorIndices = i_type::begin(); tensorIndices != i_type::end(); tensorIndices++) {
		size_t const tensorIndex = tensorIndices.flatIndex();
		if (!haveVisited[tensorIndex]) {
			haveVisited[tensorIndex] = true;
			unsigned const commonMapIndex = elementCount++;
			indexMap[tensorIndex] = { +1, commonMapIndex };
			for (size_t n = 0; n < N; n++) {
				auto const &symmetrySign = permutations[n].first;
				auto const &permutation = permutations[n].second;
				auto permutedIndices = permutation.apply(tensorIndices);
				do {
					size_t const mapIndex = permutedIndices.flatIndex();
					haveVisited[mapIndex] = true;
					signed char const sign = ((symmetrySign < 0) && (permutation.parity() < 0)) ? ((permutedIndices == tensorIndices) ? 0 : -1) : +1;
					indexMap[mapIndex] = { sign, commonMapIndex };
					permutedIndices = permutation.apply(permutedIndices);
				} while (permutedIndices != tensorIndices);
			}
		}
	}
	return indexMap;
}

template<typename T, typename F, size_t D, size_t R, Symmetries<R> S, char ...I>
struct TensorExpr {
	using index_type = IndexTuple<D, R>;
	using value_type = T;
	TensorExpr(F const &h) :
			handle(h) {
	}
	auto operator()(auto ...is) const {
		return handle(is...);
	}
	auto& operator()(auto ...is) {
		return handle(is...);
	}
	template<typename F1, Symmetries<R> S1, char ...I1>
	auto& operator=(TensorExpr<T, F1, D, R, S1, I1...> const &other) {
		std::array<char, R> A = { I... };
		std::array<char, R> B = { I1... };
		static constexpr auto P = computePermutation<R, A, B>();
		static constexpr auto indicesList = computeUniqueIndices<D, R, S1>();
		static constexpr auto permutedIndicesList = permuteList<decltype(indicesList), R, P>();
		for (size_t i = 0; i < indicesList.size(); i++) {
			handle[i] = other.handle[permutedIndicesList[i]];
		}
	}
private:
	F handle;
};

template<typename T, size_t D, size_t R, Symmetries<R> S>
struct Tensor {
	using value_type = T;
	using index_type = IndexTuple<D, R>;

	template<typename ...Args, std::enable_if<(std::is_integral<Args>::value && ...), int>::type = 0>
	T& operator()(Args ...is) {
		static thread_local T zero;
		std::array<size_t, sizeof...(is)> const indices_ = { is... };
		index_type indices(indices_);
		auto const index = IndexMap[indices.flatIndex()];
		switch (index.sign) {
		case +1:
			return +data[IndexMap[index.index]];
		case -1:
			return -data[IndexMap[index.index]];
		case 0:
			zero = T(0);
			return zero;
		}
	}

	template<typename ...Args, std::enable_if<(std::is_integral<Args>::value && ...), int>::type = 0>
	T const& operator()(Args ...is) const {
		index_type const indices = { is... };
		auto const index = IndexMap[indices.flatIndex()];
		switch (index.sign) {
		case +1:
			return +data[IndexMap[index.index]];
		case -1:
			return -data[IndexMap[index.index]];
		case 0:
			return T(0);
		}
	}

	template<char ...I>
	auto operator()(Index<I> ...) const {
		auto const handle = [this](auto ...is) {
			return this->operator()(is...);
		};
		return TensorExpr<T, decltype(handle), D, R, S, I...>(handle);
	}

	template<typename ...Args>
	auto operator()(Args ...args) {
		std::tuple<Args...> const argTuple(args...);
		auto const handle = [this, argTuple](auto ...is) {
			auto const intTuple = createIndexTuple(argTuple, is...);
			return std::apply(intTuple, [this](auto ...is) {
				return this->operator()(is...);
			});
		};
		auto const chars = extractChars(args...);
		return createTensorExpression<decltype(handle)>(handle, chars);
	}

private:
	template<typename F, char ...C>
	auto createTensorExpression(F const &handle, std::tuple<Index<C> ...>) {
		return TensorExpr<T, F, D, R, S, C...>(handle);
	}
	static constexpr size_t Size = Math::integerPower(D, R);
	static constexpr auto IndexMap = createIndexMap<D, R, S>();
	std::array<T, Size> data;
};

}

