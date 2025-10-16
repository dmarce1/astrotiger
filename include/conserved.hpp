/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_SRHD_CONSERVED_HPP_
#define INCLUDE_SRHD_CONSERVED_HPP_

#include <functional>
#include "autodiff.hpp"

#include "matrix.hpp"
#include "vector.hpp"
#include "eos.hpp"
#include "constants.hpp"
#include "math.hpp"
#include "srhd.hpp"

struct PrimitiveRecoveryFailure: public std::runtime_error {
	explicit PrimitiveRecoveryFailure(std::string const &msg, double Dval, double tauVal, double S2val, double fval, double betaVal, int iter) :
			std::runtime_error(makeMessage(msg, Dval, tauVal, S2val, fval, betaVal, iter)) {
	}
private:
	static std::string makeMessage(std::string const &msg, double Dval, double tauVal, double S2val, double fval, double betaVal, int iter) {
		std::ostringstream out;
		out << msg << "\n" << "  D      = " << Dval << "\n" << "  τ      = " << tauVal << "\n" << "  |S|²   = " << S2val << "\n" << "  f(β)   = " << fval
				<< "\n" << "  β      = " << betaVal << "\n" << "  iter   = " << iter;
		return out.str();
	}
};

template<typename Type, int dimensionCount>
struct GasConserved {
	MassDensityType<Type> D;
	EnergyDensityType<Type> τ;
	EntropyDensityType<Type> K;
	Vector<MomentumDensityType<Type>, dimensionCount> S;
	GasPrimitive<Type, dimensionCount> toPrimitive(EquationOfState<Type> const &eos) const {
		using namespace Constants;
		static constexpr FwdAutoDiff<DimensionlessType<Type>, 2, std::tuple<DimensionlessType<Type>>> one(DimensionlessType<Type>(1.0));
		static constexpr double toler = 2 * eps;
		static constexpr DimensionlessType<Type> Φo(1e-3);
		static constexpr int maxIter = 100;
		using std::abs;
		using std::copysign;
		using std::ignore;
		using std::max;
		using std::min;
		using std::tanh;
		GasPrimitive<Type, dimensionCount> prim;
		MassDensityType<Type> &ρ = prim.ρ;
		SpecificEnergyType<Type> &ε = prim.ε;
		Vector<DimensionlessType<Type>, dimensionCount> &β = prim.β;
		auto const S2 = S.dot(S);
		auto const D2 = sqr(D);
		auto const S1 = sqrt(S2);
		auto const τo = c * (sqrt(c2 * D2 + S2) - c * D);
		if (S2.value() < tiny) {
			auto const iρ = 1 / D;
			ρ = D;
			ε = max(τ, τo) * iρ;
			β = S * iρ * ic;
		} else {
			using AutoEnergyDensity = FwdAutoDiff<EnergyDensityType<Type>, 2, std::tuple<DimensionlessType<Type>>>;
			using AutoMassDensity = FwdAutoDiff<MassDensityType<Type>, 2, std::tuple<DimensionlessType<Type>>>;
			using AutoDimensionless = FwdAutoDiff<DimensionlessType<Type>, 2, std::tuple<DimensionlessType<Type>>>;
			using AutoSpecificEnergy = FwdAutoDiff<SpecificEnergyType<Type>, 2, std::tuple<DimensionlessType<Type>>>;
			using Function = std::function<AutoDimensionless(AutoDimensionless)>;
			AutoEnergyDensity W;
			AutoMassDensity ρ_;
			AutoSpecificEnergy ε_;
			AutoDimensionless ϰ_;
			Function const fEnergy = [this, eos, S1, S2, τo, &W, &ϰ_, &ε_, &ρ_](AutoDimensionless β) {
				auto const p = (AutoEnergyDensity(c * S1) / β - AutoEnergyDensity(c2 * D + max(τ, τo)));
				AutoDimensionless const β2 = β * β;
				W = AutoEnergyDensity(c2 * D + max(τ, τo)) + p;
				auto const iγ2 = one - β2;
				auto const γ2 = one / iγ2;
				auto const γ = sqrt(γ2);
				auto const iγ = one / γ;
				ρ_ = AutoMassDensity(D) * iγ;
				auto const iρ = one / ρ_;
				auto const h = W * iγ2 * iρ;
				ε_ = h - c2 - p * iρ;
				AutoEnergyDensity const peos = eos.energy2pressure(ρ_, ε_);
				return AutoDimensionless((p - peos) / W);
			};
			Function const fEntropy = [this, eos, S2, D2](AutoDimensionless β) {
				AutoDimensionless const one(DimensionlessType(1.0));
				auto const β2 = β;
				auto const γ = one / sqrt(one - β2);
				auto const ρ = AutoMassDensity(D) / γ;
				auto const ϰ(K / D);
				auto const Tg(eos.entropy2temperature(ρ, ϰ));
				auto const e = eos.temperature2energy(ρ, Tg);
				auto const ε = AutoSpecificEnergy(e);
				auto const p = eos.energy2pressure(ρ, ε);
				auto const h = c2 + ε + p / ρ;
				auto const h2 = h * h;
				auto const γ2 = γ * γ;
				return AutoDimensionless(γ2 - (c2 * S2) / (D2 * h2) - one);
			};
			bool useEntropy = false;
			REDO_PRIMITIVE: double β1, βmin, βmax, f;
			Function F;
			β1 = tanh(double(c * S1 / (τo + c2 * D)));
			βmin = β1 * eps;
			βmax = 0.9;
			β1 = min(β1, βmax);
			β1 = max(β1, βmin);
			F = useEntropy ? fEntropy : fEnergy;
			auto fRoot = F(AutoDimensionless::template independent<0>(DimensionlessType(β1)));
			f = fRoot.template get<0>().value();
			int iter = 0;
			while (abs(f) > toler) {
				if (++iter > maxIter || !std::isfinite(f)) {
					throw PrimitiveRecoveryFailure("GasConserved::toPrimitive failed to converge", D.value(), τ.value(), S2.value(), f, β1, iter);
				}
				double const dfdβ = fRoot.template get<1>().value();
				double dβ = -f / dfdβ;
				βmin = 0.5 * β1;
				βmax = 0.5 + 0.5 * β1;
				dβ = max(min(β1 + dβ, βmax), βmin) - β1;
				fRoot = F(AutoDimensionless::template independent<0>(DimensionlessType(β1 + dβ)));
				auto const f0 = f;
				f = fRoot.template get<0>().value();
				β1 += dβ;
				while (abs(f) > abs(f0)) {
					β1 -= dβ;
					dβ *= 0.5;
					β1 += dβ;
					fRoot = F(AutoDimensionless::template independent<0>(DimensionlessType(β1)));
					f = fRoot.template get<0>().value();
				}
			}
			ρ = MassDensityType<Type>(ρ_);
			β = c * S / EnergyDensityType<Type>(W);
			ε = SpecificEnergyType<Type>(ε_);
			if (!useEntropy) {
				if (prim.dualEnergySwitch(eos) < Φo) {
					useEntropy = true;
					goto REDO_PRIMITIVE;
				}
			}
		}
		return prim;
	}
};

#endif /* INCLUDE_SRHD_CONSERVED_HPP_ */
