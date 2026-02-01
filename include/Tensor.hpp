/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/
#pragma once

#include "Definitions.hpp"
#include "Math.hpp"
#include "Matrix.hpp"
#include "Permutation.hpp"
#include "Rational.hpp"
#include "SparseMatrix.hpp"
#include "TensorSymmetry.hpp"
#include "Young.hpp"

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <iostream>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>

#define TENSOR_EXPRESSION_COMPOUND_ASSIGNMENT(op)                                                                                          \
	template <typename F1, std::array<char, R> C1>                                                                                         \
	Expression &operator op(Expression<F1, T, R, D, C1> const &other) {                                                                    \
		auto in = createFreeIndexMap<R, C1>();                                                                                             \
		std::array<int, R> out;                                                                                                            \
		auto const lambda1 = [this, &out, &in, other = &other]<int I>(auto const &self) {                                                  \
			if constexpr (I < R) {                                                                                                         \
				for (int i = 0; i < D; i++) {                                                                                              \
					out[I] = i;                                                                                                            \
					std::get<std::tuple_element<I, FreeIndexMapType<R, C>>>(in).second = i;                                                \
					self.template operator()<I + 1>(self);                                                                                 \
				}                                                                                                                          \
			} else {                                                                                                                       \
				std::apply(handle_, out) = std::apply(other.handle_, out);                                                                 \
			}                                                                                                                              \
		};                                                                                                                                 \
		lambda1.template operator()<0>(lambda1);                                                                                           \
		return *this;                                                                                                                      \
	}

#define TENSOR_EXPRESSION_BINARY_OPERATOR(op)                                                                                              \
	template <typename F1, std::array<char, R> C1>                                                                                         \
	auto operator op(Expression<F1, T, R, D, C1> const &other) const {                                                                     \
		auto const f = [this, other = &other](std::integral auto... args) {                                                                \
			auto iB = createFreeIndexMap<R, C1>();                                                                                         \
			std::array<int, R> const iA = {args...};                                                                                       \
			auto const lambda1 = [this, iA, &iB, other = &other]<int I>(auto const &self) {                                                \
				if constexpr (I < R) {                                                                                                     \
					std::get<std::tuple_element<I, FreeIndexMapType<R, C>>>(iB).second = iA[I];                                            \
					return self.template operator()<I + 1>(self);                                                                          \
				} else {                                                                                                                   \
					return std::apply(handle_, iA) + std::apply(other.handle_, iB);                                                        \
				}                                                                                                                          \
			};                                                                                                                             \
			lambda1.template operator()<0>(lambda1);                                                                                       \
		};                                                                                                                                 \
		return Expression<decltype(f), T, R, D, C>(f);                                                                                     \
	}

namespace Tensors {

template <char C>
struct FreeIndex {
	static constexpr char value = C;
	constexpr operator char() const {
		return C;
	}
};

template <char C>
struct FreeIndexMapElement {
	constexpr FreeIndexMapElement &operator=(int i) {
		index_ = i;
		return *this;
	}
	consteval operator int() const {
		return index_;
	}

private:
	FreeIndex<C> free_;
	int index_;
};

template <class T>
struct IsFreeIndex : std::false_type {};

template <char C>
struct IsFreeIndex<FreeIndex<C>> : std::true_type {};

template <class T>
struct IsFixedIndex : std::is_integral<T> {};

template <class T>
concept Index = IsFixedIndex<std::remove_cvref_t<T>>::value || IsFreeIndex<std::remove_cvref_t<T>>::value;

template <class T1, class... Ts>
struct IsDummyIndex {
	static constexpr int value = ((int(std::is_same_v<T1, Ts>) + ... + 0) > 1);
};

template <typename T, size_t... I>
consteval auto concatenate(std::array<T, I> const &...arrays) {
	constexpr int N = (0 + ... + I);
	std::array<T, N> out{};
	int index = 0;
	auto append = [&out, &index](auto const &a) {
		for (auto const &x : a) {
			out[index++] = x;
		}
	};
	(append(arrays), ...);
	return out;
}

template <typename F, typename T, std::integral auto R, std::integral auto D, std::array<char, R> C>
class Expression {
	F handle_;

	template <std::integral auto R1, std::array<char, R1> C1>
	static consteval auto commonIndexCount() {
		int count = 0;
		for (int i = 0; i < R; i++) {
			for (int j = 0; j < R1; j++) {
				if (C[i] == C1[j]) {
					count++;
					break;
				}
			}
		}
		return count;
	}
	template <std::integral auto R1, std::array<char, R1> C1>
	static consteval auto createFreeIndexMap() {
		auto const lambda = []<int I>(auto const &self) {
			if constexpr (I < R1) {
				return std::tuple_cat(std::tuple(FreeIndexMapElement<C1[I]>{}), self.template operator()<I + 1>(self));
			} else {
				return std::tuple{};
			}
		};
		return lambda.template operator()<0>(lambda);
	};
	template <std::integral auto R1, std::array<char, R1> C1>
	using FreeIndexMapType = decltype(createFreeIndexMap<R1, C1>());
	template <std::integral auto R1, std::array<char, R1> C1>
	static constexpr void insertFreeIndices(FreeIndexMapType<R, C> &mapA, FreeIndexMapType<R1, C1> &mapB, std::integral auto... args) {
		constexpr auto Cs = splitIndices<R, R1, C, C1>();
		constexpr int Ra = std::get<0>(Cs).size();
		constexpr int Rb = std::get<1>(Cs).size();
		constexpr int Rt = Ra + Rb;
		auto Ca = std::get<0>(Cs);
		auto Cb = std::get<1>(Cs);
		std::array<int, R + R1> indices;
		std::array<int, sizeof...(args)> free = {args...};
		auto const fill = [free, &mapA, &mapB]<int I>(auto const &self) {
			if constexpr (I < Ra) {
				std::get<std::tuple_element<I, decltype(Ca)>>(mapA) = free[I];
			} else if constexpr (I < Rt) {
				std::get<std::tuple_element<I - Ra, decltype(Cb)>>(mapB) = free[I];
			} else {
				return;
			}
			self.template operator()<I + 1>(self);
		};
		fill.template operator()<0>(fill);
	}
	template <std::integral auto R1, std::array<char, R1> C1>
	static consteval auto splitIndices() {
		constexpr auto Rc = commonIndexCount<R, R1, C, C1>();
		constexpr auto Ra = R - Rc;
		constexpr auto Rb = R1 - Rc;
		std::array<char, Ra> Ca;
		std::array<char, Rb> Cb;
		std::array<char, Rc> Cc;
		std::bitset<R> common1{};
		std::bitset<R1> common2{};
		int k = 0;
		for (int i = 0; i < R; i++) {
			for (int j = 0; j < R1; j++) {
				if (C[i] == C1[j]) {
					common1[i] = common2[j] = true;
					Cc[k++] = C[i];
					break;
				}
			}
		}
		k = 0;
		for (int i = 0; i < R; i++) {
			if (!common1[i]) {
				Ca[k++] = C1[i];
			}
		}
		k = 0;
		for (int j = 0; j < R1; j++) {
			if (!common1[j]) {
				Cb[k++] = C1[j];
			}
		}
		return std::tuple(Ca, Cb, Cc);
	}

public:
	Expression(F const &f) :
		handle_(f) {
	}
	Expression &operator*=(T value) {
		std::array<int, R> out;
		auto const lambda1 = [this, &out, value]<int I>(auto const &self) {
			if constexpr (I < R) {
				for (int i = 0; i < D; i++) {
					self.template operator()<I + 1>(self);
				}
			} else {
				std::apply(handle_, out) *= value;
			}
		};
		lambda1.template operator()<0>(lambda1);
		return *this;
	}
	Expression &operator/=(T value) {
		operator*=(T(1) / value);
		return *this;
	}
	template <typename F1, std::integral auto R1, std::array<char, R1> C1>
	auto operator*(Expression<F1, T, R1, D, C1> const &vB) const {
		constexpr auto Cs = splitIndices<R, R1, C, C1>();
		constexpr auto Ca = std::get<0>(Cs);
		constexpr auto Cb = std::get<1>(Cs);
		constexpr auto Cc = std::get<2>(Cs);
		constexpr int Ra = std::get<0>(Cs).size();
		constexpr int Rb = std::get<1>(Cs).size();
		constexpr int R2 = Ra + Rb;
		auto const f = [vB, this](std::integral auto... args) {
			auto mapA = createFreeIndexMap<R, C>();
			auto mapB = createFreeIndexMap<R1, C1>();
			auto const lambda = [&mapA, &mapB, vB, this, args...]<int I>(auto const &self) {
				using type = std::tuple_element<I, decltype(Cc)>;
				T sum = T(0);
				if constexpr (I < Cc.size()) {
					for (int i = 0; i < D; i++) {
						std::get<type>(mapA) = std::get<type>(mapB) = i;
						sum += self.template operator()<I + 1>(self);
					}
				} else {
					insertFreeIndices(mapA, mapB, args...);
					return std::apply(handle_, mapA) * std::apply(vB, mapB);
				}
			};
			lambda.template operator()<0>(lambda);
		};
		return Expression<decltype(f), T, R, D, concatenate(Ca, Cb)>(f);
	}
	auto operator+() const {
		return *this;
	}
	auto operator-() const {
		auto const f = [this](std::integral auto... args) {
			return -handle_(args...);
		};
		return Expression<decltype(f), T, R, D, C>(f);
	}
	TENSOR_EXPRESSION_BINARY_OPERATOR(+);
	TENSOR_EXPRESSION_BINARY_OPERATOR(-);
	TENSOR_EXPRESSION_COMPOUND_ASSIGNMENT(=);
	TENSOR_EXPRESSION_COMPOUND_ASSIGNMENT(+=);
	TENSOR_EXPRESSION_COMPOUND_ASSIGNMENT(-=);
};


// template <SymmetryType auto... Ss>
// constexpr auto splitSymmetries() {
//	auto tup = std::tuple();
//	((tup = std::tuple_cat(tup, std::tuple(Ss))), ...);
//	using Tuple = std::remove_cvref_t<decltype(tup)>;
//	auto const makeTuple = [tup]<int I>(auto const &self) {
//		if constexpr (I == std::tuple_size_v<Tuple>) {
//			return std::tuple();
//		} else {
//			auto const ele = std::get<I>(tup);
//			auto const nextPair = self.template operator()<I + 1>(self);
//			auto const thisTuple = std::tuple(ele);
//			if (ele.sign == +1) {
//				auto const nextTuple = std::tuple_cat(thisTuple, nextPair.first);
//				return std::make_pair(nextTuple, nextPair.second);
//			} else if (ele.sign == -1) {
//				auto const nextTuple = std::tuple_cat(thisTuple, nextPair.second);
//				return std::make_pair(nextPair.first, nextTuple);
//			}
//		}
//	};
//	return makeTuple.template operator()<0>(makeTuple);
// }

//	constexpr auto ySym = youngSymmetrizer<Λ>();
//	SparseMatrix<Rational> A(N);
//	for (unsigned i = 0; i < ySym.size(); i++) {
//		for (auto ns = IndexType::ibegin(); ns != IndexType::iend(); ns++) {
//			auto const c = ySym[i].first;
//			auto const p = ySym[i].second;
//			auto const ms = apply<R>(p, ns);
//			A(Index(ns), Index(ms)) += c;
//		}
//	}
//	rankReduce(A);

template <typename T, std::integral auto R, std::integral auto D, SymmetryType auto... Ss>
class Tensor {
	std::array<T, ipow(D, R)> data_;

	template <typename... Args>
	static auto computeStart(Args... args) {
		constexpr auto R1 = sizeof...(Args);
		constexpr auto strides = createStrides<R1>();
		int i = 0;
		int start = 0;
		((IsFixedIndex<Args>::value ? (start += args * strides[i], i++) : i++), ...);
		return start;
	};
	template <typename... Args>
	static consteval auto createChars() {
		return concatenate<char>([]() {
			if constexpr (IsFreeIndex<Args>::value && !IsDummyIndex<Args, Args...>::value) {
				return std::array<char, 1>{Args::value};
			} else {
				return std::array<char, 0>{};
			}
		}()...);
	}
	template <typename... Args>
	static consteval auto createDummyStrides() {
		constexpr int R1 = sizeof...(Args);
		constexpr auto strides = createStrides<R1>();
		constexpr auto N = (0 + ... + int(multiplicity<Args, Args...>() > 1)) / 2;
		std::array<int, N> dummyStrides;
		int i = 0;
		int j = 0;
		(
			[strides, &i, &j, &dummyStrides]() {
				constexpr int k = find<Args, Args...>();
				if ((k != j) && multiplicity<Args, Args...>()) {
					dummyStrides[i++] = strides[j] + strides[k];
				}
				j++;
			}(),
			...);
		return dummyStrides;
	}
	template <typename... Args>
	static consteval auto createFreeStrides() {
		constexpr int R1 = sizeof...(Args);
		constexpr int N = (0 + ... + int(IsFreeIndex<Args>::value && !IsDummyIndex<Args, Args...>::value));
		constexpr std::array<bool, R1> flags = {IsFreeIndex<Args>::value && !IsDummyIndex<Args, Args...>::value...};
		constexpr auto strides = createStrides<R1>();
		std::array<int, N> freeStrides;
		for (int i = 0, j = 0; i < R1; i++) {
			if (flags[i]) {
				freeStrides[j++] = strides[i];
			}
		}
		return freeStrides;
	};
	template <std::integral auto R1>
	static consteval auto createStrides() {
		std::array<int, R1> strides;
		int stride = 1;
		for (int i = 0; i < R1; i++) {
			strides[R1 - 1 - i] = stride;
			stride *= D;
		}
		return strides;
	};
	template <class T1, class... Ts>
	static consteval int find() {
		int i = 0;
		int pos = -1;
		((std::is_same<T1, Ts>::value ? (pos = i, i++) : i++), ...);
		return pos;
	}
	template <class T1, class... Ts>
	static consteval int multiplicity() {
		return (0 + ... + int(std::is_same_v<T1, Ts>));
	}

public:
	template <typename... Args>
	auto operator()(Args... args) const {
		constexpr int R1 = (0 + ... + int(IsFreeIndex<Args>::value && !IsDummyIndex<Args, Args...>::value));
		auto const lambda = [this, args...](std::integral auto... i) {
			constexpr auto freeStrides = createFreeStrides<Args...>();
			constexpr auto dummyStrides = createDummyStrides<Args...>();
			std::array<int, sizeof...(i)> indices = {i...};
			auto const sum = [this, indices, freeStrides]<int I>(auto const &self, int start) {
				if constexpr (I == 0) {
					return data_[std::inner_product(indices.begin(), indices.end(), freeStrides.begin(), start)];
				} else {
					T result = T(0);
					for (int d = 0; d < D; d++) {
						result += self(self.template operator()<dummyStrides.size() - 1>, start + dummyStrides[d]);
					}
					return result;
				}
			};
			sum(sum.template operator()<dummyStrides.size()>, computeStart(args...));
		};
		return Expression<decltype(lambda), T, R1, D, createChars<Args...>()>(lambda);
	}
	auto &operator()(std::integral auto... i) {
		int j = 0;
		((j = D * j + i), ...);
		return data_[j];
	}
	template <char... C>
	auto operator()(FreeIndex<C>...) const {
		return Expression<Tensor, T, R, D, C...>(*this);
	}
};

} // namespace Tensors
// namespace Tensors
