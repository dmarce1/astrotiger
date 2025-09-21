/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_RADIATION_HPP_
#define INCLUDE_RADIATION_HPP_

#include "AutoDiff.hpp"
#include "Matrix.hpp"
#include "Units.hpp"
#include <iostream>
#include <climits>
#include <functional>
#include <tuple>

// Α,α,Β,β,Γ,γ,Δ,δ,Ε,ε,Ζ,ζ,Η,η,Θ,θ,ϑ,Ι,ι,Κ,κ,ϰ,Λ,λ,Μ,μ,Ν,ν,Ξ,ξ,Ο,ο,Π,π,ϖ,Ρ,ρ,ϱ,Σ,σ,ς,Τ,τ,Υ,υ,Φ,φ,ϕ,Χ,χ,Ψ,ψ,Ω,ω

template<typename T1>
struct FwdAutoDiff {
private:
	T1 v_;
	T1 d_;
	constexpr FwdAutoDiff(T1 value, T1 derivative) :
			v_(value), d_(derivative) {
	}
public:
	explicit constexpr operator T1() const {
		return v_;
	}
	constexpr T1 D() const {
		return d_;
	}
	friend constexpr FwdAutoDiff operator+(FwdAutoDiff const &A) {
		return A;
	}
	friend constexpr FwdAutoDiff operator-(FwdAutoDiff const &A) {
		FwdAutoDiff C;
		C.v_ = -A.v_;
		C.d_ = -A.d_;
		return C;
	}
	friend constexpr FwdAutoDiff operator+(FwdAutoDiff const &A, FwdAutoDiff const &B) {
		FwdAutoDiff C;
		C.v_ = A.v_ + B.v_;
		C.d_ = A.d_ + B.d_;
		return C;
	}
	friend constexpr FwdAutoDiff operator-(FwdAutoDiff const &A, FwdAutoDiff const &B) {
		FwdAutoDiff C;
		C.v_ = A.v_ - B.v_;
		C.d_ = A.d_ - B.d_;
		return C;
	}
	template<typename T2, typename T3 = decltype(T1() * T2())>
	constexpr FwdAutoDiff<T3> operator*(FwdAutoDiff<T2> const &B) const {
		FwdAutoDiff<T3> A;
		A.v_ = v_ * B.v_;
		A.d_ = d_ * B.v_ + v_ * B.d_;
		return A;
	}
//	template<typename T2, typename T3 = decltype(T1() * T2())>
//	FwdAutoDiff<T3> operator*(T2 const &b) const {
//		FwdAutoDiff<T3> A;
//		A.v_ = b * v_;
//		A.d_ = b * d_;
//		return A;
//	}
//	template<typename T2>
//	friend auto operator*(T2 const &a, FwdAutoDiff const &B) {
//		return B * a;
//	}
	constexpr FwdAutoDiff& operator*=(FwdAutoDiff const &A) {
		*this = *this * A;
		return *this;
	}
	constexpr FwdAutoDiff& operator*=(T1 const &a) {
		*this = *this * a;
		return *this;
	}
	template<typename T2, typename T3 = decltype(T1() / T2())>
	constexpr FwdAutoDiff<T3> operator/(FwdAutoDiff<T2> const &B) const {
		FwdAutoDiff<T3> A;
		A.v_ = v_ / B.v_;
		A.d_ = (d_ * B.v_ - v_ * B.d_) / sqr(B.v_);
		return A;
	}
//	template<typename T2>
//	auto operator/(T2 const &b) const {
//		return *this * (T1(1) / b);
//	}
//	template<typename T2, typename T3 = decltype(T2() / T1())>
//	friend FwdAutoDiff<T3> operator/(T2 const &v, FwdAutoDiff const &B) {
//		FwdAutoDiff<T3> A;
//		A.v_ = v / B.v_;
//		A.d_ = -B.d_ * A.v_ / B.v_;
//		return A;
//	}
	constexpr FwdAutoDiff& operator/=(FwdAutoDiff const &A) {
		*this = *this / A;
		return *this;
	}
	constexpr FwdAutoDiff& operator/=(T1 const &a) {
		*this = *this / a;
		return *this;
	}
	constexpr FwdAutoDiff& operator+=(FwdAutoDiff const &A) {
		*this = *this + A;
		return *this;
	}
	constexpr FwdAutoDiff& operator-=(FwdAutoDiff const &A) {
		*this = *this - A;
		return *this;
	}
	friend constexpr FwdAutoDiff sqr(FwdAutoDiff const &x) {
		return x * x;
	}
	friend constexpr FwdAutoDiff sqrt(FwdAutoDiff const &B) {
		FwdAutoDiff A;
		A.v_ = std::sqrt(B.v_);
		A.d_ = 0.5 * B.d_ / A.v_;
		return A;
	}
	constexpr FwdAutoDiff() {
		v_ = T1(0.0);
		d_ = T1(0.0);
	}
	constexpr FwdAutoDiff(T1 value) {
		v_ = value;
		d_ = T1(0.0);
	}
	template<typename >
	friend class FwdAutoDiff;
	static constexpr FwdAutoDiff independentVariable(T1 const &var) {
		return FwdAutoDiff(var, 1.0);
	}
	static constexpr FwdAutoDiff constant(T1 const &var) {
		return FwdAutoDiff(var, 0.0);
	}
};

namespace Radiation {

static constexpr int ndim = 3;
using DimensionlessType = Quantity<Unit<0, 0, 0, 0>>;
using LengthType = Quantity<Unit<1, 0, 0, 0>>;
using MassType = Quantity<Unit<0, 1, 0, 0>>;
using TimeType = Quantity<Unit<0, 0, 1, 0>>;
using TemperatureType = Quantity<Unit<0, 0, 0, 1>>;
using VelocityType = Quantity<Unit<1, 0, -1, 0>>;
using MomentumType = Quantity<Unit<1, 1, -1, 0>>;
using EnergyType = Quantity<Unit<2, 1, -2, 0>>;
using MassDensityType = Quantity<Unit<-3, 1, 0, 0>>;
using MomentumDensityType = Quantity<Unit<-2, 1, -1, 0>>;
using SpecificEnergyType = Quantity<Unit<2, 0, -2, 0>>;
using EnergyPerKelvinType = Quantity<Unit<2, 1, -2, -1>>;
using EnergyDensityType = Quantity<Unit<-1, 1, -2, 0>>;
using EnergyDensityPerKelvin4Type = Quantity<Unit<-1, 1, -2, -4>>;
using SpecificCrossSectionType = Quantity<Unit<2, -1, 0, 0>>;

#define GAS_CONSERVED_PROPERTIES(arg) \
		XXX(D, MassDensityType, arg, 1) \
		XXX(τ, EnergyDensityType, arg, 1) \
		XXX_FINAL(S, MomentumDensityType, arg, 3)

#define GAS_PRIMITIVE_PROPERTIES(arg) \
		XXX(ρ, MassDensityType, arg, 1) \
		XXX(v, VelocityType, arg, 3) \
		XXX_FINAL(ε, SpecificEnergyType, arg, 1)

#define RADIATION_CONSERVED_PROPERTIES(arg) \
		XXX(Er, EnergyDensityType, arg, 1) \
		XXX_FINAL(F, EnergyDensityType, arg, 3)

#define OPACITY_PROPERTIES(arg) \
		XXX(κₐ, SpecificCrossSectionType, arg, 1) \
		XXX_FINAL(κₛ, SpecificCrossSectionType, arg, 1)

#define MATERIAL_PROPERTIES(arg) \
		XXX(μ, MassType, arg, 1) \
		XXX_FINAL(Γ, DimensionlessType, arg, 1)

#define XXX_FINAL XXX

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
		return SpecificEnergyType(zeroQ) / MassDensityType(oneQ);
	}
	constexpr MassDensityType d_dε(MassDensityType const &ρ, SpecificEnergyType const &ε) const {
		return (Γ - 1) * ρ;
	}
	constexpr auto d2_dε2(MassDensityType const &ρ, SpecificEnergyType const &ε) const {
		return MassDensityType(zeroQ) / SpecificEnergyType(oneQ);
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
static constexpr EnergyPerKelvinType kB = 1.380658e-16 * cgs::g * cgs::cm2 / (cgs::s2 * cgs::K);
static constexpr EnergyDensityPerKelvin4Type aR = 7.5657e-15 * cgs::erg / (cgs::cm3 * cgs::K4);
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
		for (int k = 0; k <= ndim; k++) {
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

// Α,α,Β,β,Γ,γ,Δ,δ,Ε,ε,Ζ,ζ,Η,η,Θ,θ,ϑ,Ι,ι,Κ,κ,ϰ,Λ,λ,Μ,μ,Ν,ν,Ξ,ξ,Ο,ο,Π,π,ϖ,Ρ,ρ,ϱ,Σ,σ,ς,Τ,τ,Υ,υ,Φ,φ,ϕ,Χ,χ,Ψ,ψ,Ω,ω
#define XXX(name, type, arg, dim) copy_const_t<decltype(arg), std::conditional_t<dim == 1, type, ColumnVector<type, dim>>>& name = arg.name;

inline GasPrimitive gasCon2Prim(GasConserved const &con, Material const &mat) {
	using std::abs;
	using std::min;
	using std::max;
	using std::tanh;
	using std::ignore;
	constexpr double toler = 2 * eps;
	static int k = 0;
	GasPrimitive prim;
	GAS_CONSERVED_PROPERTIES(con);
	MATERIAL_PROPERTIES(mat)
	GAS_PRIMITIVE_PROPERTIES(prim);
	auto const S2 = vectorDotProduct(S, S);
	auto const S1 = sqrt(S2);
	if (S2.value() == 0.0) {
		auto const iρ = 1 / D;
		ρ = D;
		ε = τ * iρ;
		v = S * iρ;
	} else {
		auto const S1 = sqrt(S2);
		FwdAutoDiff<EnergyDensityType> W;
		FwdAutoDiff<MassDensityType> _ρ;
		FwdAutoDiff<SpecificEnergyType> _ε;
		auto const F = [Γ, D, τ, S1, S2, &W, &_ε, &_ρ](double β) {
			using AutoType = FwdAutoDiff<double>;
			constexpr FwdAutoDiff<double> one = 1.0;
			auto const _β = AutoType::independentVariable(β);
			auto const p = (FwdAutoDiff(c * S1) / _β - FwdAutoDiff(c2 * D + τ));
			auto const β2 = sqr(_β);
			W = FwdAutoDiff(c2 * D + τ) + p;
			auto const γ = sqrt(one / (one - β2));
			auto const γ2 = sqr(γ);
			auto const iγ = one / γ;
			auto const iγ2 = one / γ2;
			_ρ = FwdAutoDiff(D) * iγ;
			auto const iρ = one / _ρ;
			auto const h = W * iγ2 * iρ;
			_ε = h - c2 - p * iρ;
			auto const f = p - FwdAutoDiff(Γ - 1.0) * _ρ * _ε;
			return f;
		};
		double error = 1.0;
		double βmin = tiny;
		double βmax = 1.0 - eps;
		double β = tanh(double(c * S1 / (τ + c2 * D)));
		β = std::min(β, βmax);
		β = std::max(β, βmin);
		double f0 = 0.0;
		double f1 = 0.0;
		double dβ0 = 1.0;
		for (int j = 0; error > toler; j++, k++) {
			auto const f = F(β);
			double dβ = -double(EnergyDensityType(f) / f.D());
			f1 = EnergyDensityType(f).value();
			double const df_dβ = EnergyDensityType(f.D()).value();
			dβ = max(min(β + dβ, βmax), βmin) - β;
			β += dβ;
			error = double(abs(dβ / β));
			auto df_dβ_test = (f1 - f0) / dβ0;
			f0 = f1;
			dβ0 = dβ;
			printf("%i %i %e %e |  %e  %e  %e\n", j, k, error, β, df_dβ_test, df_dβ, df_dβ_test/df_dβ);
		}
		ρ = MassDensityType(_ρ);
		ε = SpecificEnergyType(_ε);
		v = c2 * S / EnergyDensityType(W);
		printf( "\n");
	}
	return prim;
}

inline GasConserved gasPrim2Con(GasPrimitive const &prim, Material const &mat) {
	GasConserved con;
	GAS_PRIMITIVE_PROPERTIES(prim);
	MATERIAL_PROPERTIES(mat)
	GAS_CONSERVED_PROPERTIES(con);
	EquationOfState pEos(mat);
	auto const v2 = vectorDotProduct(v, v);
	auto const β2 = double(ic2 * v2);
	auto const γ2 = 1 / (1 - β2);
	auto const γ = sqrt(γ2);
	auto const iρ = 1 / ρ;
	auto const p = pEos(ρ, ε);
	auto const h = c2 + ε + p * iρ;
	S = ρ * γ2 * h * ic2 * v;
	D = ρ * γ;
	τ = ρ * γ2 * ε + (γ2 - 1) * p + (D * c2 * (γ - 1));
	return con;
}

inline auto pressureTensor(RadConserved const &rad) {
	using std::max;
	RADIATION_CONSERVED_PROPERTIES(rad);
	auto const magF = vectorMagnitude(F);
	auto const imagF = 1 / max(magF, tiny * cgs::erg / cgs::cm3);
	auto const iEr = 1 / Er;
	auto const f = double(magF * iEr);
	auto const n = F * imagF;
	double const f2 = sqr(f);
	double const ξ = (3 + 4 * f2) / (5 + 2 * sqrt(4 - 3 * f2));
	double const Ωdif = 0.5 * (1 - ξ);
	double const Ωstr = 0.5 * (3 * ξ - 1);
	auto const P = Er * (Ωdif * δ + Ωstr * vectorDyadicProduct(n, n));
	return P;
}

inline double lorentzFactor(GasPrimitive const &prim) {
	GAS_PRIMITIVE_PROPERTIES(prim);
	auto const β = v * ic;
	auto const β2 = vectorDotProduct(β, β);
	auto const γ = sqrt(1 / (1 - β2));
	return double(γ);
}

template<typename QuantityType>
VectorType<QuantityType> velocity(GasPrimitive const &prim) {
	GAS_PRIMITIVE_PROPERTIES(prim);
	auto const γ = lorentzFactor(prim);
	return γ * space2SpaceTime<QuantityType>(c, v);
}

inline TensorType<double> boostTensor(GasPrimitive const &prim) {
	TensorType<double> Λ;
	GAS_PRIMITIVE_PROPERTIES(prim);
	auto const β = ColumnVector<double, ndim>(v * ic);
	double const γ = lorentzFactor(prim);
	auto const Λij = SquareMatrix<double, ndim>(δ + (sqr(γ) / (γ + 1)) * vectorDyadicProduct(β, β));
	auto const Λoi = γ * β;
	double const Λoo = γ;
	return space2SpaceTime<double>(Λoo, Λoi, Λij);
}

auto radiationTensor(RadConserved const &rad) {
	RADIATION_CONSERVED_PROPERTIES(rad);
	auto const P = pressureTensor(rad);
	return space2SpaceTime<typename std::remove_cvref<decltype(Er)>::type>(Er, F, P);
}

inline auto implicitRadiationSolve(GasConserved &gas, RadConserved &rad, Opacity const &opac, Material const &mat, cgs::Seconds const &dt) {
	GasPrimitive prim;
	GAS_CONSERVED_PROPERTIES(gas);
	RADIATION_CONSERVED_PROPERTIES(rad);
	GAS_PRIMITIVE_PROPERTIES(prim);
	MATERIAL_PROPERTIES(mat)
	OPACITY_PROPERTIES(opac)
	prim = gasCon2Prim(gas, mat);
	auto const χ = κₐ + κₛ;
	auto const implicitEnergySolve = [Γ, κₐ, μ, ρ, dt](auto &Er, auto &ε) {
		auto const ρε0 = ρ * ε;
		auto const ρκₐ = ρ * κₐ;
		auto const cdt = c * dt;
		double const λ = double(ρκₐ * cdt);
		auto const a0 = (1 - λ) * (Er + ρε0) - Er;
		auto const a1 = ρ * (λ - 1);
		auto const a4 = λ * aR * sqr(sqr((μ * (Γ - 1)) / (kB * nₐ)));
		for (int j = 0; j < 20; j++) {
			auto const ε2 = sqr(ε);
			auto const ε3 = ε * sqr(ε);
			auto const ε4 = sqr(ε2);
			auto const f = a0 + a1 * ε + a4 * ε4;
			auto const dfdε = a1 + 4 * a4 * ε3;
			auto const dε = -f / dfdε;
			ε += dε;
		}
		Er += ρε0 - ρ * ε;
	};
	auto const implicitFluxSolve = [ρ, χ, dt](auto &F) {
		auto const λ = ρ * c * χ * dt;
		F = F / (1 + λ);
	};
	auto const Rlab = radiationTensor(rad);
	auto const iΛ = boostTensor(prim);
	auto const Λ = η * iΛ * η;
	auto const Rco = Λ * Rlab * Λ;
	auto zeroMom = 0.0 * cgs::g / (cgs::cm2 * cgs::s);
	auto Sco = ColumnVector<decltype(zeroMom), ndim>(zeroMom);
	auto εco = ε;
	auto Eco = spaceTime2Scalar(Rco);
	auto Fco = spaceTime2Vector(Rco);
	implicitEnergySolve(Eco, εco);
	implicitFluxSolve(Fco);
	S += Sco;
	ε = εco;
	Er = Eco;
}

}
#endif /* INCLUDE_RADIATION_HPP_ */
