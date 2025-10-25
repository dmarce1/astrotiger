/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
*******************************************************************************/

#ifndef INCLUDE_UNITS_CGS_HPP_
#define INCLUDE_UNITS_CGS_HPP_

#include "quantity.hpp"


namespace cgs {
using null = Unit<0, 0, 0, 0, 0>;
using cm = Unit<1, 0, 0, 0, 0>;
using s = Unit<0, 1, 0, 0, 0>;
using g = Unit<0, 0, 1, 0, 0>;
using K = Unit<0, 0, 0, 1, 0>;
using mol = Unit<0, 0, 0, 0, 1>;
using erg_per_cm3 = Unit<-1, 1, -2, 0, 0>;
using g_per_cm3 = Unit<-3, 1, 0, 0, 0>;
using erg_per_g = Unit<2, 0, -2, 0, 0>;
}

template<typename T>
using DimensionlessType = Quantity<Unit<0, 0, 0, 0, 0>, T>;

template<typename T>
using LengthType = Quantity<Unit<1, 0, 0, 0, 0>, T>;

template<typename T>
using NumberDensityType = Quantity<Unit<-3, 0, 0, 0, 0>, T>;

template<typename T>
using MassType = Quantity<Unit<0, 1, 0, 0, 0>, T>;

template<typename T>
using MolarMassType = Quantity<Unit<0, 1, 0, 0, -1>, T>;

template<typename T>
using TimeType = Quantity<Unit<0, 0, 1, 0, 0>, T>;

template<typename T>
using TemperatureType = Quantity<Unit<0, 0, 0, 1, 0>, T>;

template<typename T>
using VelocityType = Quantity<Unit<1, 0, -1, 0, 0>, T>;

template<typename T>
using MomentumType = Quantity<Unit<1, 1, -1, 0, 0>, T>;

template<typename T>
using TorqueType = Quantity<Unit<2, 1, -1, 0, 0>, T>;

template<typename T>
using EnergyType = Quantity<Unit<2, 1, -2, 0, 0>, T>;

template<typename T>
using MomentumDensityType = Quantity<Unit<-2, 1, -1, 0, 0>, T>;

template<typename T>
using EnergyDensityType = Quantity<Unit<-1, 1, -2, 0, 0>, T>;

template<typename T>
using EnergyFluxType = Quantity<Unit<0, 1, -3, 0, 0>, T>;

template<typename T>
using EnergyFluxTypeSquared = Quantity<Unit<0, 2, -6, 0, 0>, T>;

template<typename T>
using FluxFluxType = Quantity<Unit<1, 1, -4, 0, 0>, T>;

template<typename T>
using SpecificEnergyType = Quantity<Unit<2, 0, -2, 0, 0>, T>;

template<typename T>
using PowerType = Quantity<Unit<2, 1, -3, 0, 0>, T>;

template<typename T>
using EntropyType = Quantity<Unit<2, 1, -2, -1, 0>, T>;

template<typename T>
using EntropyDensityType = Quantity<Unit<-1, 1, -2, -1, 0>, T>;

template<typename T>
using EntropyFluxType = Quantity<Unit<0, 1, -3, -1, 0>, T>;

template<typename T>
using MassDensityType = Quantity<Unit<-3, 1, 0, 0, 0>, T>;

template<typename T>
using SpecificEntropyType = Quantity<Unit<2, 0, -2, -1, 0>, T>;

template<typename T>
using EntropyDensityPerKelvinCubedType = Quantity<Unit<-1, 1, -2, -4, 0>, T>;

template<typename T>
using AreaType = Quantity<Unit<2, 0, 0, 0, 0>, T>;

template<typename T>
using SpecificAreaType = Quantity<Unit<2, -1, 0, 0, 0>, T>;

template<typename T>
using VolumeType = Quantity<Unit<3, 0, 0, 0, 0>, T>;

template<typename T>
using SpecificVolumeType = Quantity<Unit<3, 0, -1, 0, 0>, T>;

template<typename T>
using GravityConstantType = Quantity<Unit<3, -2, -1, 0, 0>, T>;

template<typename T>
using MoleType = Quantity<Unit<0, 0, 0, 0, 1>, T>;

template<typename T>
using PerMoleType = Quantity<Unit<0, 0, 0, 0, -1>, T>;

#endif /* INCLUDE_UNITS_CGS_HPP_ */
