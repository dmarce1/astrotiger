/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
*******************************************************************************/

#ifndef INCLUDE_UNITS1_UNITS_HPP_
#define INCLUDE_UNITS1_UNITS_HPP_

#include "numbers.hpp"


template<typename>
struct IsUnit;

template<Rational, Rational, Rational, Rational, Rational>
struct Unit;

template<typename UnitA, std::enable_if_t<IsUnit<UnitA>::value, int> = 0>
struct UnitInverse;

template<typename UnitA, Rational, std::enable_if_t<IsUnit<UnitA>::value, int> = 0>
struct UnitPower;

template<typename UnitA, typename UnitB, std::enable_if_t<IsUnit<UnitA>::value && IsUnit<UnitB>::value, int> = 0>
struct UnitProduct;

template<typename UnitA, typename UnitB, std::enable_if_t<IsUnit<UnitA>::value && IsUnit<UnitB>::value, int> = 0>
struct UnitQuotient;



template<typename Type>
struct IsUnit {
	static constexpr bool value = false;
};

template<Rational L, Rational M, Rational T, Rational K, Rational N>
struct IsUnit<Unit<L, M, T, K, N>> {
	static constexpr bool value = true;
};

template<Rational L, Rational M, Rational T, Rational K, Rational N>
struct Unit {
	static constexpr Rational L_ = L;
	static constexpr Rational M_ = M;
	static constexpr Rational T_ = T;
	static constexpr Rational K_ = K;
	static constexpr Rational N_ = N;
};

template<typename UnitA, std::enable_if_t<IsUnit<UnitA>::value, int>>
struct UnitInverse {
	using type = typename UnitPower<UnitA, -1>::type;
};

template<typename UnitA, Rational power, std::enable_if_t<IsUnit<UnitA>::value, int>>
struct UnitPower {
	using type = Unit<power * UnitA::L_, power * UnitA::M_, power * UnitA::T_, power * UnitA::K_, power * UnitA::N_>;
};

template<typename UnitA, typename UnitB, std::enable_if_t<IsUnit<UnitA>::value && IsUnit<UnitB>::value, int>>
struct UnitProduct {
	using type = Unit<UnitA::L_ + UnitB::L_, UnitA::M_ + UnitB::M_, UnitA::T_ + UnitB::T_, UnitA::K_ + UnitB::K_, UnitA::N_ + UnitB::N_>;
};

template<typename UnitA, typename UnitB, std::enable_if_t<IsUnit<UnitA>::value && IsUnit<UnitB>::value, int>>
struct UnitQuotient {
	using type = typename UnitProduct<UnitA, typename UnitInverse<UnitB>::type>::type;
};


using NullUnitType = Unit<0, 0, 0, 0, 0>;

#endif /* INCLUDE_UNITS_UNITS_HPP_ */
