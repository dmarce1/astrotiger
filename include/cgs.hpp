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

using DimensionlessType = Quantity<Unit<0, 0, 0, 0, 0>, double>;
using LengthType = Quantity<Unit<1, 0, 0, 0, 0>, double>;
using NumberDensityType = Quantity<Unit<-3, 0, 0, 0, 0>, double>;
using MassType = Quantity<Unit<0, 1, 0, 0, 0>, double>;
using MolarMassType = Quantity<Unit<0, 1, 0, 0, -1>, double>;
using TimeType = Quantity<Unit<0, 0, 1, 0, 0>, double>;
using TemperatureType = Quantity<Unit<0, 0, 0, 1, 0>, double>;
using VelocityType = Quantity<Unit<1, 0, -1, 0, 0>, double>;
using MomentumType = Quantity<Unit<1, 1, -1, 0, 0>, double>;
using MomentumDensityType = Quantity<Unit<-2, 1, -1, 0, 0>, double>;
using TorqueType = Quantity<Unit<2, 1, -1, 0, 0>, double>;
using EnergyType = Quantity<Unit<2, 1, -2, 0, 0>, double>;
using EnergyDensityType = Quantity<Unit<-1, 1, -2, 0, 0>, double>;
using SpecificEnergyType = Quantity<Unit<2, 0, -2, 0, 0>, double>;
using PowerType = Quantity<Unit<2, 1, -3, 0, 0>, double>;
using EnergyFluxType = Quantity<Unit<0, 1, -3, 0, 0>, double>;
using EntropyType = Quantity<Unit<2, 1, -2, -1, 0>, double>;
using EntropyDensityType = Quantity<Unit<-1, 1, -2, -1, 0>, double>;
using MassDensityType = Quantity<Unit<-3, 1, 0, 0, 0>, double>;
using SpecificEntropyType = Quantity<Unit<2, 0, -2, -1, 0>, double>;
using EntropyDensityPerKelvinCubedType = Quantity<Unit<-1, 1, -2, -4, 0>, double>;
using AreaType = Quantity<Unit<2, 0, 0, 0, 0>, double>;
using SpecificAreaType = Quantity<Unit<2, -1, 0, 0, 0>, double>;
using VolumeType = Quantity<Unit<3, 0, 0, 0, 0>, double>;
using SpecificVolumeType = Quantity<Unit<3, 0, -1, 0, 0>, double>;
using MoleType = Quantity<Unit<0, 0, 0, 0, 1>, double>;
using PerMoleType = Quantity<Unit<0, 0, 0, 0, -1>, double>;

#endif /* INCLUDE_UNITS_CGS_HPP_ */
