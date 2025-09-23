/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_RADIATION_HPP_
#define INCLUDE_RADIATION_HPP_

#include "FwdAutoDiff.hpp"
#include "Units.hpp"
#include <iostream>
#include <climits>
#include <functional>
#include <tuple>

// Α,α,Β,β,Γ,γ,Δ,δ,Ε,ε,Ζ,ζ,Η,η,Θ,θ,ϑ,Ι,ι,Κ,κ,ϰ,Λ,λ,Μ,μ,Ν,ν,Ξ,ξ,Ο,ο,Π,π,ϖ,Ρ,ρ,ϱ,Σ,σ,ς,Τ,τ,Υ,υ,Φ,φ,ϕ,Χ,χ,Ψ,ψ,Ω,ω
// ₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₒₓₕₖₗₘₙₚₛₜⱼ



namespace Radiation {

static constexpr int ndim = 3;
using DimensionlessType = Quantity<Unit<0, 0, 0, 0>>;
using LengthType = Quantity<Unit<1, 0, 0, 0>>;
using MassType = Quantity<Unit<0, 1, 0, 0>>;
using MassDensityType = Quantity<Unit<-3, 1, 0, 0>>;
using TimeType = Quantity<Unit<0, 0, 1, 0>>;
using TemperatureType = Quantity<Unit<0, 0, 0, 1>>;
using VelocityType = Quantity<Unit<1, 0, -1, 0>>;
using MomentumType = Quantity<Unit<1, 1, -1, 0>>;
using MomentumDensityType = Quantity<Unit<-2, 1, -1, 0>>;
using EnergyType = Quantity<Unit<2, 1, -2, 0>>;
using EnergyDensityType = Quantity<Unit<-1, 1, -2, 0>>;
using SpecificEnergyType = Quantity<Unit<2, 0, -2, 0>>;
using EntropyType = Quantity<Unit<2, 1, -2, -1>>;
using EntropyDensityType = Quantity<Unit<-1, 1, -2, -1>>;
using SpecificEntropyType = Quantity<Unit<2, 0, -2, -1>>;
using EntropyDensityPerKelvinCubed = Quantity<Unit<-1, 1, -2, -4>>;
using AreaType = Quantity<Unit<2, 0, 0, 0>>;
using SpecificAreaType = Quantity<Unit<2, -1, 0, 0>>;
using VolumeType = Quantity<Unit<3, 0, 0, 0>>;
using SpecificVolumeType = Quantity<Unit<3, 0, -1, 0>>;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

#define GAS_CONSERVED_PROPERTIES(arg) \
		XXX(D, MassDensityType, arg, 1) \
		XXX(τ, EnergyDensityType, arg, 1) \
		XXX(S, MomentumDensityType, arg, 3)

#define GAS_PRIMITIVE_PROPERTIES(arg) \
		XXX(ρ, MassDensityType, arg, 1) \
		XXX(ε, SpecificEnergyType, arg, 1)  \
		XXX(v, VelocityType, arg, 3)

#define RADIATION_CONSERVED_PROPERTIES(arg) \
		XXX(Er, EnergyDensityType, arg, 1) \
		XXX(F, EnergyDensityType, arg, 3)

#define OPACITY_PROPERTIES(arg) \
		XXX(κₐ, SpecificAreaType, arg, 1) \
		XXX(κₛ, SpecificAreaType, arg, 1)

#define MATERIAL_PROPERTIES(arg) \
		XXX(μ, MassType, arg, 1) \
		XXX(Γ, DimensionlessType, arg, 1)

#pragma GCC diagnostic pop

#define XXX(name, type, arg, dim) std::conditional_t<dim == 1, type, ColumnVector<type, dim>> name;
struct GasConserved {
	GAS_CONSERVED_PROPERTIES(0)
};
struct GasPrimitive {
	GAS_PRIMITIVE_PROPERTIES(0)
};
struct RadConserved {
	RADIATION_CONSERVED_PROPERTIES(0)
};
struct Opacity {
	OPACITY_PROPERTIES(0)
};
struct Material {
	MATERIAL_PROPERTIES(0)
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

struct EquationOfState {
#define XXX(name, type, arg, dim) name = arg.name;
	constexpr EquationOfState(Material const &mat) {
		MATERIAL_PROPERTIES(mat)
	}
#undef XXX
	constexpr EnergyDensityType operator()(MassDensityType const &ρ, SpecificEnergyType const &ε) const {
		return (Γ - 1) * ρ * ε;
	}
	constexpr SpecificEnergyType d_dρ(MassDensityType const &ρ, SpecificEnergyType const &ε) const {
		return (Γ - 1) * ε;
	}
	constexpr auto d2_dρ2(MassDensityType const &ρ, SpecificEnergyType const &ε) const {
		return SpecificEnergyType(0.0) / MassDensityType(1.0);
	}
	constexpr MassDensityType d_dε(MassDensityType const &ρ, SpecificEnergyType const &ε) const {
		return (Γ - 1) * ρ;
	}
	constexpr auto d2_dε2(MassDensityType const &ρ, SpecificEnergyType const &ε) const {
		return MassDensityType(0.0) / SpecificEnergyType(1.0);
	}
	constexpr auto d2_dρdε(MassDensityType const &ρ, SpecificEnergyType const &ε) const {
		return (Γ - 1);
	}
private:
#define XXX(name, type, arg, dim) type name;
	MATERIAL_PROPERTIES(0) /**/
#undef XXX
};

static constexpr int mantissaWidth = std::numeric_limits<double>::digits;
static constexpr double π = 4 * atan(1);
static constexpr double nₐ = 6.02214076e23;
static constexpr VelocityType c = 2.99792458e10 * cgs::cm / cgs::s;
static constexpr EntropyType kB = 1.380658e-16 * cgs::g * cgs::cm2 / (cgs::s2 * cgs::K);
static constexpr EntropyDensityPerKelvinCubed aR = 7.5657e-15 * cgs::erg / (cgs::cm3 * cgs::K4);
static constexpr auto c2 = sqr(c);
static constexpr auto ic = 1.0 / c;
static constexpr auto ic2 = sqr(ic);
static constexpr double nDigits = std::numeric_limits<double>::digits;
static constexpr double tiny = sqrt(std::numeric_limits<double>::min());
static constexpr double huge = sqrt(std::numeric_limits<double>::max());
static constexpr double eps = std::numeric_limits<double>::epsilon();
static constexpr auto δ = SquareMatrix<double, ndim>::identity();
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

GasPrimitive gasCon2Prim(GasConserved const&, Material const&);
GasConserved gasPrim2Con(GasPrimitive const&, Material const&);
void implicitEnergySolve(GasPrimitive&, RadConserved&, Opacity const&, Material const&, cgs::Seconds const&);
void implicitRadiationSolve(GasConserved&, RadConserved&, Opacity const&, Material const&, cgs::Seconds const&);


#define XXX(name, type, arg, dim) copy_const_t<decltype(arg), std::conditional_t<dim == 1, type, ColumnVector<type, dim>>>& name = arg.name;

inline SpecificEnergyType temperature2gasEnergy(Material const &mat, TemperatureType const &T) {
	MATERIAL_PROPERTIES(mat);
	return kB * nₐ / (μ * (Γ - 1)) * T;
}


inline EnergyDensityType temperature2radiationEnergy(TemperatureType const &T) {
	return aR * sqr(sqr(T));
}
#undef XXX

}
#endif /* INCLUDE_RADIATION_HPP_ */
