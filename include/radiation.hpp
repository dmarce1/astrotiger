/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once
#include "constants.hpp"
#include "opacity.hpp"
#include "eos.hpp"
#include "fpe.hpp"

#include <iostream>
#include <climits>
#include <functional>
#include <tuple>
#include "autodiff.hpp"
#include "gas_conserved.hpp"
#include "gas_primitive.hpp"
#include "tensor.hpp"

template<typename Type>
struct ImplicitEnergyFailure: public std::runtime_error {
	explicit ImplicitEnergyFailure(std::string const &msg, Type rhoVal, Type epsVal, Type EVal, Type fVal, Type dfVal, Type stepVal, int iter) :
			std::runtime_error(makeMessage(msg, rhoVal, epsVal, EVal, fVal, dfVal, stepVal, iter)) {
	}
private:
	static std::string makeMessage(std::string const &msg, Type rhoVal, Type epsVal, Type EVal, Type fVal, Type dfVal, Type stepVal, int iter) {
		std::ostringstream out;
		out << msg << "\n" << "  ρ       = " << rhoVal << "\n" << "  ε       = " << epsVal << "\n" << "  E_r     = " << EVal << "\n" << "  f       = " << fVal
				<< "\n" << "  df/de   = " << dfVal << "\n" << "  Δe      = " << stepVal << "\n" << "  iter    = " << iter;
		return out.str();
	}
};

template<typename Type, int dimensionCount>
struct RadFlux {
	EnergyFluxType<Type> E;
	Vector<FluxFluxType<Type>, dimensionCount> F;
	// @formatter:off
	DEFINE_VECTOR_OPERATORS(RadFlux, Type, a.E, a.F)
																																																					// @formatter:on
};
template<typename Type, int dimensionCount>
struct RadConserved {
	// @formatter:off
	DEFINE_VECTOR_OPERATORS(RadConserved, Type, a.E, a.F)
																																																					// @formatter:on
	auto temperature() const {
		return pow<Rational(1, 4)>(E / pc.aR);
	}
	auto pressure() const {
		using std::max;
		auto const magF = sqrt(F.dot(F));
		auto const imagF = 1 / max(magF, EnergyFluxType<Type>(tiny));
		auto const iE = 1 / E;
		auto const f = magF * iE / pc.c;
		auto const n = F * imagF;
		auto const f2 = sqr(f);
		auto const ξ = (3 + 4 * f2) / (5 + 2 * sqrt(4 - 3 * f2));
		auto const Ωdif = 0.5 * (1 - ξ);
		auto const Ωstr = 0.5 * (3 * ξ - 1);
		auto const P = E * (Ωdif * δ + Ωstr * sqr(n));
		jacobian(1);
		return P;
	}
	auto stressEnergy() const {
		return space2spaceTime<EnergyDensityType<Type>, dimensionCount>(E, F / pc.c, pressure());
	}
	auto flux(int ni) const {
		enableFPE();
		using std::max;
		constexpr auto c2 = sqr(pc.c);
		RadFlux<Type, dimensionCount> flux;
		auto const un = unitVector<Type, dimensionCount>(ni);
		auto const magF = sqrt(F.dot(F));
		auto const imagF = 1 / max(magF, EnergyFluxType<Type>(tiny));
		auto const iE = 1 / E;
		auto const f = magF * iE / pc.c;
		auto const N = F * imagF;
		auto const f2 = sqr(f);
		auto const ξ = (3 + 4 * f2) / (5 + 2 * sqrt(4 - 3 * f2));
		auto const Ωdif = 0.5 * (1 - ξ);
		auto const Ωstr = 0.5 * (3 * ξ - 1);
		auto const P = E * (Ωstr * N * N[ni] + Ωdif * un);
		flux.E = F[ni];
		for (int k = 0; k < dimensionCount; k++) {
			flux.F = c2 * P[ni][k];
		}
		return flux;
	}
	auto jacobian(int ni) const {
		enableFPE();
		using std::max;
		using Auto = FwdAutoDiff<Type, 1, std::array<EnergyDensityType<Type>, fieldCount>>;
		constexpr auto _1 = Auto(Type(1));
		constexpr auto _2 = Auto(Type(2));
		constexpr auto _3 = Auto(Type(3));
		constexpr auto _4 = Auto(Type(4));
		constexpr auto _5 = Auto(Type(5));
		constexpr auto half = _1 / _2;
		constexpr auto δ = identity<Auto, dimensionCount>();
		constexpr int ti = dimensionCount;
		constexpr auto c = pc.c;
		SquareMatrix<Type, fieldCount> J;
		Vector<Auto, fieldCount> flux;
		Vector<Auto, dimensionCount> F_;
		auto const un = unitVector<Auto, dimensionCount>(ni);
		auto E_ = Auto::independent(Type(E), ti);
		for (int k = 0; k < dimensionCount; k++) {
			F_[k] = Auto::independent(Type(F[k] / c), k);
		}
		auto const magF = sqrt(F_.dot(F_));
		auto const f = magF / E_;
		auto const n = F_ / max(magF, Auto(tiny));
		auto const f2 = sqr(f);
		auto const χ = (_3 + _4 * f2) / (_5 + _2 * sqrt(_4 - _3 * f2));
		auto const P = half * E_ * ((_1 - χ) * δ + (_3 * χ - _1) * n * n);
		auto const fluxF = P * un;
		flux[ti] = Type(c) * F_[ni];
		for (int k = 0; k < dimensionCount; k++) {
			flux[k] = Type(c) * fluxF[k];
		}
		for (int n = 0; n < fieldCount; n++) {
			auto const dFdU = flux[n].gradient();
			for (int k = 0; k < fieldCount; k++) {
				J(n, k) = dFdU[k];
			}
		}
		return J;
	}
	auto eigenSystem(int ni) const {
		enableFPE();
		using std::max;
		constexpr auto c = pc.c;
		auto const tj = transverseIndices[ni];
		constexpr auto fi = 0;
		constexpr auto li = fieldCount - 1;
		constexpr auto ei = dimensionCount;
		auto const ti = transverseIndices[ni];
		SquareMatrix<Type, fieldCount> R;
		Vector<Type, fieldCount> λ;
		Vector<EnergyFluxType<Type>, fieldCount> flux;
		Vector<DimensionlessType<Type>, transverseCount> βᵗ, nᵗ;
		auto const magF = sqrt(F.dot(F));
		auto const f = magF / (c * E);
		auto const f2 = sqr(f);
		auto const ξ = (3 + 4 * f2) / (5 + 2 * sqrt(4 - 3 * f2));
		auto const H = E * (3 - ξ) / 2;
		auto const β = F / (c * H);
		auto const βˣ = β[ni];
		for (int j = 0; j < transverseCount; j++) {
			nᵗ[j] = βᵗ[j] = β[ti[j]];
		}
		auto const βᵗβᵗ = βᵗ.dot(βᵗ);
		auto const βᵗ1 = sqrt(βᵗβᵗ);
		if (βᵗ1 > Type(0)) {
			nᵗ /= βᵗ1;
		}
		auto const βˣβˣ = sqr(βˣ);
		auto const β2 = βᵗβᵗ + βˣβˣ;
		auto const Vᵗ = sqrt(H);
		auto const V = Vᵗ * β;
		auto const Vˣ = V[ni];
		auto const V2 = V.dot(V);
		for (int dim = 0; dim < dimensionCount; dim++) {
			if (ni != dim) {
				flux[dim] = c * V[dim] * Vˣ;
			} else {
				flux[dim] = c * (Type(0.25) * (sqr(Vᵗ) - V2) + sqr(Vˣ));
			}
		}
		auto const λ1 = 2 * βˣ / (3 - β2);
		auto const λ2 = sqrt((1 - β2) * (3 - β2 - 2 * βˣβˣ)) / (3 - β2);
		auto const λm = λ1 - λ2;
		auto const λo = βˣ;
		auto const λp = λ1 + λ2;
		auto const λp2 = sqr(λp);
		auto const λm2 = sqr(λm);
		R(ei, li) = Type(1 - 2 * βˣ * λp + λp2);
		R(ni, li) = Type(2 * λp - (1 + λp2) * βˣ);
		for (int j = 0; j < transverseCount; ++j) {
			R(tj[j], li) = Type((1 - λp2) * βᵗ[j]);
		}
		R(ei, fi) = Type(1 - 2 * βˣ * λm + λm2);
		R(ni, fi) = Type(2 * λm - (1 + λm2) * βˣ);
		for (int j = 0; j < transverseCount; ++j) {
			R(tj[j], fi) = Type((1 - λm2) * βᵗ[j]);
		}
		if constexpr (dimensionCount >= 2) {
			R(ni, 1) = Type(βˣ * βᵗ1);
			R(ei, 1) = Type(βᵗ1);
			if (βᵗ1 > Type(0)) {
				for (int j = 0; j < transverseCount; ++j) {
					R(tj[j], 1) = Type((1 - βˣβˣ) * nᵗ[j]);
				}
			} else {
				for (int j = 0; j < transverseCount; ++j) {
					R(tj[j], 1) = Type(j == 0);
				}
			}

			if constexpr (dimensionCount == 3) {
				R(ni, 2) = Type(0);
				R(ei, 2) = Type(0);
				if (βᵗ1 > Type(0)) {
					R(tj[0], 2) = Type(+nᵗ[1]);
					R(tj[1], 2) = Type(-nᵗ[0]);
				} else {
					R(tj[0], 2) = Type(0);
					R(tj[1], 2) = Type(1);
				}
			}
		}
		SquareMatrix<Type, fieldCount> A(Type(0));
		for (int j = 0; j < dimensionCount; ++j) {
			A(j, j) = Type(Vᵗ);      // ∂U_j/∂V_j = Vt
			A(j, ei) = Type(V[j]);     // ∂U_j/∂Vt = V_j
			A(ei, j) = Type(V[j] / 2); // ∂U_t/∂V_j = V_j/2
		}
		A(ei, ei) = Type(3 * Vᵗ / 2);    // ∂U_t/∂Vt = 3/2 Vt
		R = A * R;
		λ.front() = Type(c * λm);
		λ.back() = Type(c * λp);
		for (int k = 1; k < fieldCount - 1; k++) {
			λ[k] = Type(c * λo);
		}
		return std::tuple(λ, R);
	}
	RadConserved<Type, dimensionCount> implicitRadiationSolve(GasConserved<Type, dimensionCount> &gasCon, Opacity<Type> const &opac,
			EquationOfState<Type> const &eos, TimeType<Type> const &dt) const {
		auto const implicitEnergySolve = [](SpecificEnergyType<Type> &ε, EnergyDensityType<Type> &E, MassDensityType<Type> const &ρ_,
				SpecificAreaType<Type> const &κₐ, EquationOfState<Type> const &eos, TimeType<Type> const &dt) {
			using std::abs;
			using std::copysign;
			using std::max;
			using std::min;
			constexpr Type toler = 16 * eps;
			constexpr int maxIter = 100;
			using Dimensionless = FwdAutoDiff<DimensionlessType<Type>, 2, std::tuple<DimensionlessType<Type>>>;
			using MassDensity = FwdAutoDiff<MassDensityType<Type>, 2, std::tuple<DimensionlessType<Type>>>;
			using EnergyDensity = FwdAutoDiff<EnergyDensityType<Type>, 2, std::tuple<DimensionlessType<Type>>>;
			using SpecificEnergy = FwdAutoDiff<SpecificEnergyType<Type>, 2, std::tuple<DimensionlessType<Type>>>;
			constexpr auto one = Dimensionless(DimensionlessType(1.0));
			auto const ρ = MassDensity(MassDensityType<Type>(ρ_));
			auto const Eg0 = EnergyDensity(ρ * ε);
			auto const Er0 = EnergyDensity(E);
			auto const Etot0 = EnergyDensity(Er0 + Eg0);
			auto const λ = Dimensionless(ρ * pc.c * κₐ * dt);
			auto const Fgas = [eos, Etot0, Er0, Eg0, λ, ρ, one](SpecificEnergyType<Type> ε_) {
				auto const ε = SpecificEnergy::template independent<0>(ε_);
				auto const T = eos.energy2temperature(ρ, ε);
				auto const T2 = sqr(T);
				auto const T4 = sqr(T2);
				auto const B = pc.aR * T4;
				auto const E = Etot0 - ρ * ε;
				auto const f = ((one + λ) * E - (Er0 + λ * B)) / (Er0 + Eg0);
				return f;
			};
			auto const Frad = [eos, Etot0, Er0, Eg0, λ, ρ, one](EnergyDensityType<Type> E_) {
				auto const E = EnergyDensity::template independent<0>(E_);
				auto const ρε = Etot0 - E;
				auto const ε = ρε / ρ;
				auto const T = eos.energy2temperature(ρ, ε);
				auto const T2 = sqr(T);
				auto const T4 = sqr(T2);
				auto const B = pc.aR * T4;
				auto const f = (λ * B - ((one + λ) * E - Er0)) / (Er0 + Eg0);
				return f;
			};
			Type error = 1.0;
			int iter = 0;
			for (; error > toler; ++iter) {
				if (iter > maxIter) {
					throw ImplicitEnergyFailure("implicitEnergySolve failed to converge", ρ_.value(), ε.value(), E.value(),
							std::numeric_limits<Type>::quiet_NaN(), std::numeric_limits<Type>::quiet_NaN(), std::numeric_limits<Type>::quiet_NaN(), iter);
				}
				bool gasForm = (ρ_ * ε) < E;
				Dimensionless fRoot = gasForm ? Fgas(ε) : Frad(E);
				auto const f = fRoot.template get<0>().value();
				auto const dfde = fRoot.template get<1>().value();
				auto const d2fde2 = fRoot.template get<2>().value();
				if (!std::isfinite(f) || !std::isfinite(dfde) || !std::isfinite(d2fde2)) {
					throw ImplicitEnergyFailure("NaN or Inf detected in implicitEnergySolve derivatives", ρ_.value(), ε.value(), E.value(), f, dfde, 0.0, iter);
				}
				auto const den = Type(max(sqr(dfde), 2 * sqr(dfde) - f * d2fde2));
				auto const num = -2 * f * dfde;
				auto const de = num / (den + copysign(tiny, den));
				if (!std::isfinite(de)) {
					throw ImplicitEnergyFailure("Non-finite update step in implicitEnergySolve", ρ_.value(), ε.value(), E.value(), f, dfde, de, iter);
				}
				auto ε1 = ε;
				auto E1 = E;
				if (gasForm) {
					ε1 += SpecificEnergyType<Type>(de);
					E1 = (Etot0 - ρ * ε1).template get<0>();
				} else {
					E1 += EnergyDensityType<Type>(de);
					ε1 = ((Etot0 - E1) / ρ_).template get<0>();
				}
				ε = max(0.5 * ε, ε1);
				E = max(0.5 * E, E1);
				Type const e = gasForm ? ε.value() : E.value();
				error = abs((de / e));
				if (!std::isfinite(error) || error > 1e10) {
					throw ImplicitEnergyFailure("Divergence detected in implicitEnergySolve", ρ_.value(), ε.value(), E.value(), f, dfde, de, iter);
				}
			}
		};
		using std::max;
		using std::min;
		constexpr Type toler = 4 * eps;
		enableFPE();
		auto radCon = *this;
		GasPrimitive<Type, dimensionCount> gasPrim;
		EnergyDensityType<Type> &τ = gasCon.τ;
		EnergyDensityType<Type> &E = radCon.E;
		EntropyDensityType<Type> &K = gasCon.K;
		MassDensityType<Type> &D = gasCon.D;
		MassDensityType<Type> &ρ = gasPrim.ρ;
		SpecificEnergyType<Type> &ε = gasPrim.ε;
		Vector<VelocityType<Type>, dimensionCount> &v = gasPrim.v;
		Vector<EnergyFluxType<Type>, dimensionCount> &F = radCon.F;
		Vector<MomentumDensityType<Type>, dimensionCount> &S = gasCon.S;
		auto const norm = sqrt(sqr(sqrt(F.dot(F))) / sqr(pc.c) + sqr(E) + sqr(pc.c) * sqr(sqrt(S.dot(S))) + sqr(τ)).value();
		Type error = 1.0;
		auto const F0 = F;
		auto const E0 = E;
		auto const S0 = S;
		auto const τ0 = τ;
		auto R0 = radCon.stressEnergy();
		for (int j = 0; error > toler; j++) {
			auto const F1 = F;
			auto const E1 = E;
			gasPrim = gasCon.toPrimitive(eos);
			auto const γ = gasPrim.lorentzFactor();
			auto const dt0 = dt / γ;
			auto const Λ = space2spaceTime<DimensionlessType<Type>>(γ, γ * v / pc.c, δ + (sqr(γ) / (γ + 1)) * sqr(v / pc.c));
			auto const iΛ = space2spaceTime<DimensionlessType<Type>>(γ, -γ * v / pc.c, δ + (sqr(γ) / (γ + 1)) * sqr(v / pc.c));
			auto const λ = ρ * pc.c * (opac.κₐ + opac.κₛ) * dt0;
			auto R = symmetric(Λ * R0 * Λ);
			E = spaceTime2Tensor0(R);
			F = pc.c * spaceTime2Tensor1(R);
			implicitEnergySolve(ε, E, ρ, opac.κₛ, eos, dt0);
			F = F / (1 + λ);
			R = radCon.stressEnergy();
			R = symmetric(iΛ * R * iΛ);
			E = spaceTime2Tensor0(R);
			F = pc.c * spaceTime2Tensor1(R);
			S = S0 + (F - F0) / sqr(pc.c);
			auto const τo = sqr(pc.c) * (γ - 1.0) * D;
			auto const τmin = 0.5 * (τo + τ);
			τ = τ0 + E0 - E;
			τ = max(τmin, τ);
			auto const T = eos.energy2temperature(ρ, ε);
			K = D * eos.temperature2entropy(ρ, T);
			error = Type(sqrt(sqr(sqrt((F - F1).dot(F - F1))) / sqr(pc.c) + sqr(E - E1)).value() / norm);
		}
		return radCon;
	}
	void setEnergy(EnergyDensityType<Type> const &E_) {
		E = E_;
	}
	void setFlux(Vector<EnergyFluxType<Type>, dimensionCount> const &F_) {
		F = F_;
	}
	void setFlux(int k, EnergyFluxType<Type> const &Fk) {
		F[k] = Fk;
	}
	EnergyDensityType<Type> getEnergy() const {
		return E;
	}
	Vector<EnergyFluxType<Type>, dimensionCount> getFlux() const {
		return F;
	}
	EnergyFluxType<Type> getFlux(int k) const {
		return F[k];
	}
private:
	static constexpr auto δ = identity<DimensionlessType<Type>, dimensionCount>();
	static constexpr Type tiny = sqrt(std::numeric_limits < Type > ::min());
	static constexpr Type eps = std::numeric_limits < Type > ::epsilon();
	static constexpr Type toler = 2 * eps;
	static constexpr int fieldCount = 1 + dimensionCount;
	static constexpr int transverseCount = dimensionCount - 1;
	static constexpr PhysicalConstants<Type> pc {};
	// @formatter:off
	static constexpr auto transverseIndices =
	(dimensionCount == 1) ? std::array<std::array<int, 2>, 3> { { {-1,-1}, {-1,-1}, {-1,-1}}}:
	((dimensionCount == 2) ? std::array<std::array<int, 2>, 3> { { {1,-1}, {0,-1}, {-1,-1}}}:
			/*dimensionCount == 3)*/std::array<std::array<int, 2>, 3> { { {1, 2}, {0, 2}, {0, 1}}});
	// @formatter:on
	EnergyDensityType<Type> E;
	Vector<EnergyFluxType<Type>, dimensionCount> F;
};

template<typename Type>
inline EnergyDensityType<Type> temperature2radiationEnergy(TemperatureType<Type> const &T) {
	static constexpr PhysicalConstants<Type> pc { };
	return pc.aR * sqr(sqr(T));
}

#undef XXX

