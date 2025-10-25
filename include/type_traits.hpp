#pragma once


template<typename Tuple>
struct AllSame;

template<typename T>
struct Unwrap;

template<typename Tuple>
struct UnwrapTuple;

namespace detail {

// ---------- Repeat a type T, N times, as tuple<T, T, ..., T> (type-only) ----------

template<typename T, std::size_t N, typename Seq = std::make_index_sequence<N>>
struct RepeatTypeAsTuple;

template<typename T, std::size_t N, std::size_t... I>
struct RepeatTypeAsTuple<T, N, std::index_sequence<I...>> {
    // Use the conditional_t<true, T, ...> trick to "expand" T N times:
    using type = std::tuple<
        std::conditional_t<true, std::remove_cv_t<T>, std::integral_constant<std::size_t, I>>...
    >;
};

// Public alias (still under detail)
template<typename T, std::size_t N>
using ArrayToTupleType = typename RepeatTypeAsTuple<T, N>::type;


// ---------- Tuple concatenation that tolerates non-tuple second arg ----------

// Normalize any type to a tuple type.
template<typename T>
struct ToTuple { using type = std::tuple<T>; };

template<typename... Ts>
struct ToTuple<std::tuple<Ts...>> { using type = std::tuple<Ts...>; };

// Primary helper declaration (must be in the same namespace as the specialization)
template<typename T1, typename T2>
struct ConcatTupleTypesHelper;

// Specialization for two tuples.
template<typename... Ts, typename... Us>
struct ConcatTupleTypesHelper<std::tuple<Ts...>, std::tuple<Us...>> {
    using type = std::tuple<Ts..., Us...>;
};

} // namespace detail

// Public aliases (outside detail)

// Convert array element type T and extent N to tuple<T, T, ..., T> (N times).
template<typename T, std::size_t N>
using ArrayToTupleType = typename detail::ArrayToTupleType<T, N>;

// Concatenate tuple-like types; if either side is not a tuple, it's wrapped in a 1-elem tuple.
template<typename T1, typename T2>
using ConcatTupleTypes =
    typename detail::ConcatTupleTypesHelper<
        typename detail::ToTuple<T1>::type,
        typename detail::ToTuple<T2>::type
    >::type;


// ---------- Your existing utilities (fine to keep as you had) ----------

template<typename From, typename To>
using CopyConstType = std::conditional_t<std::is_const_v<std::remove_reference_t<From>>, std::add_const_t<To>, To>;

template<typename T, typename... Ts>
struct AllSame<std::tuple<T, Ts...>> : std::conjunction<std::is_same<T, Ts>...> { };

template<>
struct AllSame<std::tuple<>> : std::true_type { };

template<typename T>
struct AllSame<std::tuple<T>> : std::true_type { };

template<typename T>
struct Unwrap { using type = T; };

template<template<typename, typename> class Q, typename U, typename T>
struct Unwrap<Q<U, T>> { using type = U; };

template<template<typename...> class Tuple, typename... Ts>
struct UnwrapTuple<Tuple<Ts...>> { using type = Tuple<typename Unwrap<Ts>::type...>; };
