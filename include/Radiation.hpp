/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_RADIATION_HPP_
#define INCLUDE_RADIATION_HPP_

#include "constants.hpp"
#include "units.hpp"
#include "numbers.hpp"

//#include "FwdAutoDiff.hpp"
#include "Matrix.hpp"
#include <iostream>
#include <climits>
#include <functional>
#include <tuple>

// Α,α,Β,β,Γ,γ,Δ,δ,Ε,ε,Ζ,ζ,Η,η,Θ,θ,ϑ,Ι,ι,Κ,κ,ϰ,Λ,λ,Μ,μ,Ν,ν,Ξ,ξ,Ο,ο,Π,π,ϖ,Ρ,ρ,ϱ,Σ,σ,ς,Τ,τ,Υ,υ,Φ,φ,ϕ,Χ,χ,Ψ,ψ,Ω,ω
// ₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₒₓₕₖₗₘₙₚₛₜⱼ

namespace Radiation {

static constexpr int ndim = 3;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

#define GAS_CONSERVED_PROPERTIES(arg) \
		XXX(D, MassDensityType, arg, 1) \
		XXX(τ, EnergyDensityType, arg, 1) \
		XXX(K, EntropyDensityType, arg, 1) \
		XXX(S, MomentumDensityType, arg, 3)

#define GAS_PRIMITIVE_PROPERTIES(arg) \
		XXX(ρ, MassDensityType, arg, 1) \
		XXX(ε, SpecificEnergyType, arg, 1)  \
		XXX(γ, DimensionlessType, arg, 1)  \
		XXX(β, DimensionlessType, arg, 3)

#define RADIATION_CONSERVED_PROPERTIES(arg) \
		XXX(Er, EnergyDensityType, arg, 1) \
		XXX(F,  EnergyFluxType, arg, 3)

#define OPACITY_PROPERTIES(arg) \
		XXX(κₐ, SpecificAreaType, arg, 1) \
		XXX(κₛ, SpecificAreaType, arg, 1)

#pragma GCC diagnostic pop

#define XXX(name, type, arg, dim) std::conditional_t<dim == 1, type, ColumnVector<type, dim>> name;

struct RadConserved {
	RADIATION_CONSERVED_PROPERTIES(0)
	auto temperature() const;
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
	constexpr EquationOfState() = default;
	constexpr EquationOfState(EquationOfState const&) = default;
	constexpr EquationOfState& operator=(EquationOfState const&) = default;

	constexpr EquationOfState(MolarMassType μ_) :
			μ(μ_) {
	}

	template<typename Tρ, typename Tε>
	auto energy2pressure(Tρ const&, Tε const&) const;

	template<typename Tρ, typename TT>
	auto temperature2entropy(Tρ const&, TT const&) const;

	template<typename Tρ, typename Tϰ>
	auto entropy2temperature(Tρ const&, Tϰ const&) const;

	template<typename Tρ, typename TT>
	auto temperature2energy(Tρ const&, TT const&) const;

	template<typename Tρ, typename Tε>
	auto energy2temperature(Tρ const&, Tε const&) const;

private:
	static constexpr Rational Γ = Rational(5, 3);
	MolarMassType μ;
};

struct GasPrimitive;

struct GasConserved {
	GAS_CONSERVED_PROPERTIES(0)
	GasPrimitive toPrimitive(EquationOfState const&) const;
};

struct GasPrimitive {
	GAS_PRIMITIVE_PROPERTIES(0)
	GasConserved toConserved(EquationOfState const&) const;
	auto dualEnergySwitch(EquationOfState const&) const;
	void updateLorentz();
	auto temperature(EquationOfState const&) const;
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

void implicitRadiationSolve(GasConserved&, RadConserved&, Opacity const&, EquationOfState const&, TimeType const&,
		DimensionlessType ξo = DimensionlessType(1e-3));

#define XXX(name, type, arg, dim) copy_const_t<decltype(arg), std::conditional_t<dim == 1, type, ColumnVector<type, dim>>>& name = arg.name;

inline EnergyDensityType temperature2radiationEnergy(TemperatureType const &T) {
	return aR * sqr(sqr(T));
}

#undef XXX

template<typename Tρ, typename Tε>
auto EquationOfState::energy2pressure(Tρ const &ρ, Tε const &ε) const {
	return ρ * ε * (Γ - 1);
}

template<typename Tρ, typename Tε>
auto EquationOfState::energy2temperature(Tρ const &ρ, Tε const &ε) const {
	return (Γ - 1) * μ * ε / (nₐ * kB);
}

template<typename Tρ, typename TT>
auto EquationOfState::temperature2entropy(Tρ const &ρ, TT const &T) const {
	constexpr Rational Cv = 1 / (Γ - 1);
	constexpr Rational Cp = Γ / (Γ - 1);
	return nₐ * kB / μ * (log(μ / (nₐ * ρ) * pow<Cv>((4 * π * Cv * kB * T * μ) / (3 * _h * _h * nₐ))) + Cp);
}

template<typename Tρ, typename Tϰ>
auto EquationOfState::entropy2temperature(Tρ const &ρ, Tϰ const &ϰ) const {
	constexpr Rational Cv = 1 / (Γ - 1);
	constexpr Rational Cp = Γ / (Γ - 1);
	auto const term1 = ((3 * _h * _h * nₐ) / (4 * π * Cv * kB * μ));
	auto const term2 = pow<1 / Cv>(nₐ * ρ / μ * exp(μ * ϰ / (nₐ * kB)) / exp(Cp));
	return term1 * term2;
}

template<typename Tρ, typename TT>
auto EquationOfState::temperature2energy(Tρ const &ρ, TT const &T) const {
	return (nₐ * kB) / ((Γ - 1) * μ) * T;
}

}

#endif /* INCLUDE_RADIATION_HPP_ */
