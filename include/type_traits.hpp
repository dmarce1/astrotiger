#pragma once

#include <tuple>
#include <type_traits>

template<typename T>
struct Unwrap;

template<template<typename, typename > class Q, typename U, typename T>
struct Unwrap<Q<U, T>> {
	using type = U;
};

template<typename Tuple>
struct UnwrapTuple;

template<template<typename ...> class Tuple, typename... Ts>
struct UnwrapTuple<Tuple<Ts...>> {
	using type = Tuple<typename Unwrap<Ts>::type...>;
};

template<typename From, typename To>
using CopyConstType = std::conditional_t<std::is_const_v<std::remove_reference_t<From>>, std::add_const_t<To>, To>;

template<typename Tuple>
struct AllSame;

template<typename T, typename ... Ts>
struct AllSame<std::tuple<T, Ts...>> : std::conjunction<std::is_same<T, Ts>...> {
};

#include <tuple>
#include <utility>
#include <type_traits>

