/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_RADIATION_HPP_
#define INCLUDE_RADIATION_HPP_

#include "FwdAutoDiff.hpp"
#include "Matrix.hpp"
#include "Units.hpp"
#include <iostream>
#include <climits>
#include <functional>
#include <tuple>

// Α,α,Β,β,Γ,γ,Δ,δ,Ε,ε,Ζ,ζ,Η,η,Θ,θ,ϑ,Ι,ι,Κ,κ,ϰ,Λ,λ,Μ,μ,Ν,ν,Ξ,ξ,Ο,ο,Π,π,ϖ,Ρ,ρ,ϱ,Σ,σ,ς,Τ,τ,Υ,υ,Φ,φ,ϕ,Χ,χ,Ψ,ψ,Ω,ω
// ₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₒₓₕₖₗₘₙₚₛₜⱼ

namespace Radiation {

static constexpr int ndim = 3;
using DimensionlessType = Quantity<Unit<0, 0, 0, 0, 0>, double>;
using LengthType = Quantity<Unit<1, 0, 0, 0, 0>, double>;
using MassType = Quantity<Unit<0, 1, 0, 0, 0>, double>;
using MassDensityType = Quantity<Unit<-3, 1, 0, 0, 0>, double>;
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
using EntropyType = Quantity<Unit<2, 1, -2, -1, 0>, double>;
using EntropyDensityType = Quantity<Unit<-1, 1, -2, -1, 0>, double>;
using SpecificEntropyType = Quantity<Unit<2, 0, -2, -1, 0>, double>;
using EntropyDensityPerKelvinCubedType = Quantity<Unit<-1, 1, -2, -4, 0>, double>;
using AreaType = Quantity<Unit<2, 0, 0, 0, 0>, double>;
using SpecificAreaType = Quantity<Unit<2, -1, 0, 0, 0>, double>;
using VolumeType = Quantity<Unit<3, 0, 0, 0, 0>, double>;
using SpecificVolumeType = Quantity<Unit<3, 0, -1, 0, 0>, double>;
using MoleType = cgs::Mole;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

#define GAS_CONSERVED_PROPERTIES(arg) \
		XXX(D, EnergyDensityType, arg, 1) \
		XXX(τ, EnergyDensityType, arg, 1) \
		XXX(K, EnergyDensityType, arg, 1) \
		XXX(S, EnergyDensityType, arg, 3)

#define GAS_PRIMITIVE_PROPERTIES(arg) \
		XXX(ρ, MassDensityType, arg, 1) \
		XXX(ε, SpecificEnergyType, arg, 1)  \
		XXX(γ, DimensionlessType, arg, 1)  \
		XXX(β, DimensionlessType, arg, 3)

#define RADIATION_CONSERVED_PROPERTIES(arg) \
		XXX(Er, EnergyDensityType, arg, 1) \
		XXX(F, EnergyDensityType, arg, 3)

#define OPACITY_PROPERTIES(arg) \
		XXX(κₐ, SpecificAreaType, arg, 1) \
		XXX(κₛ, SpecificAreaType, arg, 1)

#define EOS_PROPERTIES(arg) \
		XXX(μ, MolarMassType, arg, 1) \
		XXX(Γ, DimensionlessType, arg, 1)

#pragma GCC diagnostic pop

#define XXX(name, type, arg, dim) std::conditional_t<dim == 1, type, ColumnVector<type, dim>> name;

struct RadConserved {
	RADIATION_CONSERVED_PROPERTIES(0)
};

struct Opacity {
	OPACITY_PROPERTIES(0)
};

template<typename T>
struct Underlying {
	using type = T;
};

template<template<typename > class W, typename U>
struct Underlying<W<U>> {
	using type = U;
};

template<typename T>
using UnderlyingType = typename Underlying<T>::type;

struct EquationOfState {
	EOS_PROPERTIES(0)
	constexpr EquationOfState() = default;
	constexpr EquationOfState(EquationOfState const&) = default;
	constexpr EquationOfState& operator=(EquationOfState const&) = default;

	constexpr EquationOfState(double Γ_, double μ_) :
			μ(μ_), Γ(Γ_) {
	}

	template<typename Tρ, typename Tε>
	auto pressure(Tρ const&, Tε const&) const;

	template<typename Tρ, typename Tε>
	auto entropy(Tρ const&, Tε const&) const;

	template<typename Tρ, typename Tϰ>
	auto energy(Tρ const&, Tϰ const&) const;

	template<typename Tε>
	auto temperature(Tε const&) const;

};

struct GasPrimitive;

struct GasConserved {
	GAS_CONSERVED_PROPERTIES(0)
	GasPrimitive toPrimitive(EquationOfState const&, DimensionlessType ξo = DimensionlessType(1e-3)) const;
};

struct GasPrimitive {
	static constexpr double dualEnergySwitch = 1e-3;GAS_PRIMITIVE_PROPERTIES(0)
	GasConserved toConserved(EquationOfState const&) const;
};

#undef XXX

template<typename Type>
using TensorType = SquareMatrix<Type, ndim + 1>;

template<typename Type>
using SpaceTensorType = SquareMatrix<Type, ndim>;

template<typename Type>
using VectorType = ColumnVector<Type, ndim + 1>;

template<typename Type>
using SpaceVectorType = ColumnVector<Type, ndim>;

template<typename From, typename To>
using copy_const_t = std::conditional_t<std::is_const_v<std::remove_reference_t<From>>, std::add_const_t<To>, To>;

static constexpr int mantissaWidth = std::numeric_limits<double>::digits;
static constexpr double π = 4 * atan(1);
static constexpr auto amu = cgs::g / 6.02214076e23;
static constexpr auto nₐ = DimensionlessType(1.0) / MoleType(cgs::mol / 6.02214076e23);
static constexpr VelocityType c = 2.99792458e10 * cgs::cm / cgs::s;
static constexpr EntropyType kB = 1.380658e-16 * cgs::g * cgs::cm2 / (cgs::s2 * cgs::K);
static constexpr EntropyDensityPerKelvinCubedType aR = 7.5657e-15 * cgs::erg / (cgs::cm3 * cgs::K4);
static constexpr TorqueType ℎ = 6.62607015e-27 * cgs::g * cgs::cm2 / cgs::s;
static constexpr TorqueType ℏ = ℎ / (2 * π);
static constexpr auto c2 = sqr(c);
static constexpr auto ic = 1.0 / c;
static constexpr auto ic2 = sqr(ic);
static constexpr double nDigits = std::numeric_limits<double>::digits;
static constexpr double tiny = sqrt(std::numeric_limits<double>::min());
static constexpr double huge = sqrt(std::numeric_limits<double>::max());
static constexpr double eps = std::numeric_limits<double>::epsilon();
static constexpr auto δ = SquareMatrix<DimensionlessType, ndim>::identity();
static constexpr auto η = []() {
	SquareMatrix<double, ndim + 1> η;
	η(ndim, ndim) = -1;
	for (int j = 0; j < ndim; j++) {
		η(ndim, j) = η(ndim, j) = 0;
		for (int k = 0; k < j; k++) {
			η(j, k) = η(k, j) = 1;
		}
	}
	return η;
}(); /**/
// Α,α,Β,β,Γ,γ,Δ,δ,Ε,ε,Ζ,ζ,Η,η,Θ,θ,ϑ,Ι,ι,Κ,κ,ϰ,Λ,λ,Μ,μ,Ν,ν,Ξ,ξ,Ο,ο,Π,π,ϖ,Ρ,ρ,ϱ,Σ,σ,ς,Τ,τ,Υ,υ,Φ,φ,ϕ,Χ,χ,Ψ,ψ,Ω,ω
// ₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₒₓₕₖₗₘₙₚₛₜⱼ

template<typename QuantityType>
static TensorType<QuantityType> space2SpaceTime(QuantityType s, SpaceVectorType<QuantityType> const &v, SpaceTensorType<QuantityType> const &t) {
	TensorType<QuantityType> T;
	for (int j = 0; j < ndim; j++) {
		for (int k = 0; k < ndim; k++) {
			T(j, k) = t(k, j);
		}
		T(j, ndim) = T(ndim, j) = v[j];
	}
	T(ndim, ndim) = s;
	return T;
}

template<typename QuantityType>
static VectorType<QuantityType> space2SpaceTime(QuantityType s, SpaceVectorType<QuantityType> const &v) {
	VectorType<QuantityType> V;
	for (int j = 0; j < ndim; j++) {
		V[j] = v[j];
	}
	V[ndim] = s;
	return V;
}

template<typename QuantityType>
static QuantityType spaceTime2Scalar(TensorType<QuantityType> const &T) {
	return T(ndim, ndim);
}

template<typename QuantityType>
static SpaceVectorType<QuantityType> spaceTime2Vector(TensorType<QuantityType> const &T) {
	SpaceVectorType<QuantityType> v;
	for (int j = 0; j < ndim; j++) {
		v[j] = T(j, ndim);
	}
	return v;
}

template<typename QuantityType>
static SpaceTensorType<QuantityType> spaceTime2Tensor(TensorType<QuantityType> const &T) {
	SpaceTensorType<QuantityType> t;
	for (int j = 0; j < ndim; j++) {
		for (int k = 0; k < ndim; k++) {
			t(j, k) = T(j, k);
		}
	}
	return t;
}

void implicitRadiationSolve(GasConserved&, RadConserved&, Opacity const&, EquationOfState const&, cgs::Seconds const&, DimensionlessType ξo = DimensionlessType(1e-3));

#define XXX(name, type, arg, dim) copy_const_t<decltype(arg), std::conditional_t<dim == 1, type, ColumnVector<type, dim>>>& name = arg.name;

inline EnergyDensityType temperature2radiationEnergy(TemperatureType const &T) {
	return aR * sqr(sqr(T));
}

template<typename Type>
auto gasEnergy2temperature(EquationOfState const &mat, Type const &ε) {
	EOS_PROPERTIES(mat);
	auto const C = ((μ * (Γ - 1)) / (kB * nₐ));
	return C * ε;
}

template<typename Type>
auto gasEnergy2specificEntropyPotential(EquationOfState const &mat, Type const &ε) {
	using QuantityType = Quantity<Unit<-2, 0, 2, 1, 0>, double>;
	EOS_PROPERTIES(mat);
	using AutoType = std::conditional_t<std::is_same_v<Type, FwdAutoDiff<SpecificEnergyType>>, FwdAutoDiff<QuantityType>, QuantityType>;
	return AutoType(μ * (Γ - 1) / (kB * nₐ)) * ε;
}

template<typename Type>
inline auto temperature2gasEnergy(EquationOfState const &mat, Type const &T) {
	using QuantityType = SpecificEntropyType;
	EOS_PROPERTIES(mat);
	using AutoType = std::conditional_t<std::is_same_v<Type, FwdAutoDiff<TemperatureType>>, FwdAutoDiff<QuantityType>, QuantityType>;
	return AutoType(kB * nₐ / (μ * (Γ - 1))) * T;
}

template<typename Type>
auto radiationEnergy2temperature(Type const &Er) {
	return sqrt(sqrt(Er / aR));
}

#undef XXX

template<typename Tε>
auto EquationOfState::temperature(Tε const &ε) const {
	return ε * ((μ * (Γ - 1)) / (kB * nₐ));
}

template<typename Tρ, typename Tε>
auto EquationOfState::pressure(Tρ const &ρ, Tε const &ε) const {
	return ρ * ε * (Γ - 1);
}

template<typename Tρ, typename Tε>
auto EquationOfState::entropy(Tρ const &ρ, Tε const &ε) const {
	auto const ϰo = pow((nₐ * ℎ * ℎ / (2 * π * μ * kB)).value(), -Γ.value());
	auto const ϰ = pressure(ρ.value(), ε.value()) * pow(ρ.value(), -Γ.value());
	return DimensionlessType(ϰ / ϰo);
}

template<typename Tρ, typename Tϰ>
auto EquationOfState::energy(Tρ const &ρ, Tϰ const &ϰ) const {
	auto const ϰo = pow((nₐ * ℎ * ℎ / (2 * π * μ * kB)).value(), -Γ.value());
	return SpecificEnergyType(ϰo * ϰ.value() * pow(ρ.value(), Γ.value() - 1) / (Γ - 1).value());
}

}

#endif /* INCLUDE_RADIATION_HPP_ */
