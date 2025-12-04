/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/
#pragma once

#include "Definitions.hpp"

#include <array>
#include <bit>
#include <bitset>
#include <concepts>
#include <cstdint>
#include <ranges>
#include <type_traits>

#include "Math.hpp"

namespace Tensors {

template<typename, int, char...>
struct Expression;

template<char>
struct FreeIndex;

template<auto, auto...>
struct IndexOf;

template<auto, typename, auto...>
struct IndexOfImpl;

template<int>
struct MultiFor;

template<auto...>
struct Multiplicity;

template<typename, auto...>
struct Permute;

template<typename, int, int>
struct Tensor;

template<char ...C>
constexpr int contractionRank = (int(Multiplicity<C, C...>::value == 1) + ...);

template<int D>
struct IndexMap {
	constexpr void set(char c, uint8_t v) {
		if (c & 1) {
			map_[c >> 1].u = v;
		} else {
			map_[c >> 1].l = v;
		}
	}
	constexpr uint8_t get(char c) const {
		if (c & 1) {
			return map_[c >> 1].u;
		} else {
			return map_[c >> 1].l;
		}
	}
	constexpr IndexMap() :
			map_ { } {

	}
private:
	union NibblePair {
		uint8_t l :4;
		uint8_t u :4;
	};
	static constexpr int size = 128;
	std::array<NibblePair, size / 2> map_;
};

template<typename F1, int D, char ...C1>
struct Expression {
	static constexpr int R1 = sizeof...(C1);
	constexpr Expression(F1 handle) :
			handle_(handle) {
	}
	Expression() = delete;
	Expression(Expression const&) = default;
	auto operator()(auto ...i) const {
		return handle_(i...);
	}
	auto& operator()(auto ...i) {
		return handle_(i...);
	}
	template<typename F2, char ...C2>
	constexpr auto& operator=(Expression<F2, D, C2...> &&other_) {
		auto const other = std::move(other_);
		using Seq = typename Permute<int, C1...>::type<C2...>;
		MultiFor<D>::template call<R1>([this, other](auto ...i) {
			constexpr Seq seq { };
			(*this)(i...) = access(other, seq, i...);
		});
		return *this;
	}
private:
	template<int ...I>
	static constexpr auto access(auto const &Expr, std::integer_sequence<int, I...>, auto ...args) {
		auto const tup = std::tuple(args...);
		return Expr(std::get<I>(tup)...);
	}
	F1 const handle_;
};


template<typename F1, int D>
struct Expression<F1, D> {
	constexpr Expression(F1 handle) :
			handle_(handle) {
	}
	Expression() = delete;
	Expression(Expression const&) = default;
	operator auto() const {
		return handle_();
	}
	template<typename F2>
	constexpr auto& operator=(Expression &&other) {
		handle_ = other.handle_;
		return *this;
	}
private:
	template<int ...I>
	static constexpr auto access(auto const &Expr, std::integer_sequence<int, I...>, auto ...args) {
		auto const tup = std::tuple(args...);
		return Expr(std::get<I>(tup)...);
	}
	F1 const handle_;
};

template<char C>
struct FreeIndex {
	static constexpr char value = C;
	constexpr operator char() const {
		return C;
	}
};

template<auto C1, auto ...Cs>
struct IndexOf {
	static constexpr int value = IndexOfImpl<C1, std::integral_constant<int, 0>, Cs...>::value;
};

template<auto C1, typename I, auto C2, auto ...Cs>
struct IndexOfImpl<C1, I, C2, Cs...> {
	static constexpr int value = (C1 == C2) ? int(I()) : IndexOfImpl<C1, std::integral_constant<int, int(I()) + 1>, Cs...>::value;
};

template<auto C1, typename I>
struct IndexOfImpl<C1, I> {
	static constexpr int value = int(I());
};

template<int, char...>
struct FixedBits;

template<int I, char C0, char ...C1>
struct FixedBits<I, C0, C1...> {
	template<char ...C2>
	static constexpr int value = ((Multiplicity<C0, C2...>::value > 1) ? 0 : (1 << I)) | FixedBits<I - 1, C1...>::template value<C2...>;
};

template<int I>
struct FixedBits<I> {
	template<char ...C2>
	static constexpr int value = 0;
};

template<char ...C>
constexpr std::bitset<sizeof...(C)> fixedBits = FixedBits<sizeof...(C) - 1, C...>::template value<C...>;

template<int D>
struct MultiFor {
	template<int R, typename F, int I = 0, typename ... Args>
	static constexpr auto call(F const &f) noexcept {
		if constexpr (I < D) {
			if constexpr (R == 0) {
				f(Args::value...);
			} else {
				call<R - 1, F, 0, Args..., std::integral_constant<int, I> >(f);
				call<R, F, I + 1, Args...>(f);
			}
		}
	}
};

template<auto A, auto B, auto ...C>
struct Multiplicity<A, B, C...> {
	static constexpr int value = int(A == B) + Multiplicity<A, C...>::value;
};

template<auto A>
struct Multiplicity<A> {
	static constexpr int value = 0;
};

template<typename T, auto ...C1>
struct Permute {
	template<auto ...C3>
	using type = std::integer_sequence<T, IndexOf<C3, C1...>::value...>;
};

template<typename T, auto ...Cs>
struct FixedSequence {
	static constexpr auto value = std::tuple_cat((std::conditional_t<Multiplicity<Cs, Cs...>::value == 1,
			std::tuple<std::integral_constant<T, IndexOf<Cs, Cs...>::value>>,
			std::tuple<> >())...);
};

template<typename T, auto ...C2>
struct FreeSequence;

template<typename T, auto C1, auto ...C2>
struct FreeSequence<T, C1, C2...> {
	template<auto ...Cs>
	struct type {
		static constexpr auto value = std::tuple_cat(
				std::conditional_t<Multiplicity<C1, C2...>::value >= 1, std::tuple<std::integral_constant<T, IndexOf<C1, Cs...>::value>>, std::tuple<> >(),
				FreeSequence<T, C2...>::template type<Cs...>::value);
	};
};

template<typename T, auto C1>
struct FreeSequence<T, C1> {
	template<auto ...Cs>
	struct type {
		static constexpr auto value = std::tuple<>();
	};
};

template<typename T, auto ...Cs>
class TraceSequence {
	static constexpr auto freeIdx = FreeSequence<T, Cs...>::template type<Cs...>::value;
public:
	static constexpr auto value = std::tuple_cat(FixedSequence<T, Cs...>::value, freeIdx, freeIdx);
};

template<typename T, int R, int D>
struct Tensor {
	T operator()(auto ...i) const {
		return data_[flatten(i...)];
	}
	T& operator()(auto ...i) {
		return data_[flatten(i...)];
	}
	template<char ...C>
	auto operator()(FreeIndex<C> ...) {
		if constexpr (contractionRank<C...> == R) {
			return Expression<Tensor&, D, C...>(*this);
		} else {
			Expression<Tensor const&, D, C...> expr(*this);
			constexpr auto seq = FixedSequence<int, C...>::value;
			constexpr auto R2 = R - std::tuple_size<decltype(seq)>::value;
			return makeExpression<std::array { C... }>([expr](auto ...i) {
				T sum = T(0);
				MultiFor<D>::template call<R2 / 2>([&sum, expr, i...](auto ...j) {
					sum += access(expr, TraceSequence<int, C...>::value, i..., j..., j...);
				});
				return sum;
			}, seq);
		}
	}
private:
	template<auto Tup, typename F, int ...I>
	static auto makeExpression(F f, std::tuple<std::integral_constant<int, I>...> seq) {
		return Expression<F, D, std::get<I>(Tup)...>(f);
	}
	template<int ...I>
	static constexpr auto access(auto const &Expr, std::tuple<std::integral_constant<int, I>...>, auto ...args) {
		auto const tup = std::tuple(args...);
		return Expr(std::get<I>(tup)...);
	}
	int flatten(auto ...i) const {
		int j = 0;
		((j = D * j + i),...);
		return j;
	}
	static constexpr int size() {
		return iPow(D, R);
	}
	std::array<T, size()> data_;
};

template<typename T, int D>
struct Tensor<T, 0, D> {
	operator T() const {
		return data_;
	}
	Tensor& operator=(T const& value) const {
		data_ = value;
		return *this;
	}
private:
	static constexpr int size() {
		return 1;
	}
	T data_;
};

template<char...>
struct CommonIndices;

template<char C0, char...C1>
struct CommonIndices<C0, C1...> {
	template<char...C2>
	struct type {
		static constexpr int value = int(((C0 == C2) || ...)) + CommonIndices<C1...>::template type<C2...>::value;
	};
};

template<>
struct CommonIndices<> {
	template<char...C2>
	struct type {
		static constexpr int value = 0;
	};
};

template<typename F1, typename F2, int D, char...C1, char...C2>
auto operator*(Expression<F1, D, C1...> const&, Expression<F2, D, C2...> const&) {

}

}
