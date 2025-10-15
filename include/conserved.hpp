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



template<int NDIM>
struct GasConserved {
#define XXX(name, type, arg, dim) std::conditional_t<dim == 1, Tensor0<type, dim>, Tensor1<type, dim>> name;
	GAS_CONSERVED_PROPERTIES(0, NDIM)
#undef XXX
	GasPrimitive<NDIM> toPrimitive(EquationOfState const &eos) const {
		using namespace Constants;
		static constexpr FwdAutoDiff<DimensionlessType, NDIM, std::tuple<DimensionlessType>> one(DimensionlessType(1.0));
		static constexpr double toler = 2 * eps;
		static constexpr DimensionlessType Φo(1e-3);
		static constexpr int maxIter = 100;
		using std::abs;
		using std::copysign;
		using std::ignore;
		using std::max;
		using std::min;
		using std::tanh;
		GasPrimitive<NDIM> gasPrim;
#define XXX(name, type, arg, dim) CopyConstType<decltype(arg), std::conditional_t<dim == 1, Tensor0<type, dim>, Tensor1<type, dim>>>& name = arg.name;
		GAS_PRIMITIVE_PROPERTIES(gasPrim, NDIM);
#undef XXX
		auto const S2 = S.dot(S);
		auto const D2 = sqr(D);
		auto const S1 = sqrt(S2);
		auto const τo = c * (sqrt(c2 * D2 + S2) - c * D);
		if (S2.value() < tiny) {
			auto const iρ = 1 / D;
			γ = 1.0;
			ρ = D;
			ε = max(τ, τo) * iρ;
			β = S * iρ * ic;
		} else {
			FwdAutoDiff<EnergyDensityType, NDIM, std::tuple<DimensionlessType>> W;
			FwdAutoDiff<MassDensityType, NDIM, std::tuple<DimensionlessType>> ρ_;
			FwdAutoDiff<SpecificEnergyType, NDIM, std::tuple<DimensionlessType>> ε_;
			FwdAutoDiff<DimensionlessType, NDIM, std::tuple<DimensionlessType>> ϰ_;
			using Function = std::function<FwdAutoDiff<DimensionlessType, NDIM, std::tuple<DimensionlessType>>(FwdAutoDiff<DimensionlessType, NDIM, std::tuple<DimensionlessType>>)>;
			Function const fEnergy = [this, eos, S1, S2, τo, &W, &ϰ_, &ε_, &ρ_](FwdAutoDiff<DimensionlessType, NDIM, std::tuple<DimensionlessType>> β) {
				auto const p = (FwdAutoDiff<EnergyDensityType, NDIM, std::tuple<DimensionlessType>>(c * S1) / β - FwdAutoDiff<EnergyDensityType, NDIM, std::tuple<DimensionlessType>>(c2 * D + max(τ, τo)));
				FwdAutoDiff<DimensionlessType, NDIM, std::tuple<DimensionlessType>> const β2 = β * β;
				W = FwdAutoDiff<EnergyDensityType, NDIM, std::tuple<DimensionlessType>>(c2 * D + max(τ, τo)) + p;
				auto const iγ2 = one - β2;
				auto const γ2 = one / iγ2;
				auto const γ = sqrt(γ2);
				auto const iγ = one / γ;
				ρ_ = FwdAutoDiff<MassDensityType, NDIM, std::tuple<DimensionlessType>>(D) * iγ;
				auto const iρ = one / ρ_;
				auto const h = W * iγ2 * iρ;
				ε_ = h - c2 - p * iρ;
				auto const peos = FwdAutoDiff<EnergyDensityType, NDIM, std::tuple<DimensionlessType>>(eos.energy2pressure(ρ_, ε_));
				auto const rc = FwdAutoDiff<DimensionlessType, NDIM, std::tuple<DimensionlessType>>((p - peos) / W);
				return rc;
			};
			Function const fEntropy = [this, eos, S2, D2](FwdAutoDiff<DimensionlessType, NDIM, std::tuple<DimensionlessType>> β) {
				FwdAutoDiff<DimensionlessType, NDIM, std::tuple<DimensionlessType>> const one(DimensionlessType(1.0));
				auto const β2 = β;
				auto const γ = one / sqrt(one - β2);
				auto const ρ = FwdAutoDiff<MassDensityType, NDIM, std::tuple<DimensionlessType>>(D) / γ;
				auto const ϰ(K / D);
				auto const Tg(eos.entropy2temperature(ρ, ϰ));
				auto const e = eos.temperature2energy(ρ, Tg);
				auto const ε = FwdAutoDiff<SpecificEnergyType, NDIM, std::tuple<DimensionlessType>>(e);
				auto const p = eos.energy2pressure(ρ, ε);
				auto const h = c2 + ε + p / ρ;
				auto const h2 = h * h;
				auto const γ2 = γ * γ;
				auto const rc = FwdAutoDiff<DimensionlessType, NDIM, std::tuple<DimensionlessType>>(γ2 - (c2 * S2) / (D2 * h2) - one);
				return rc;
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
			auto fRoot = F(FwdAutoDiff<DimensionlessType, NDIM, std::tuple<DimensionlessType>>::template independent<0>(DimensionlessType(β1)));
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
				fRoot = F(FwdAutoDiff<DimensionlessType, NDIM, std::tuple<DimensionlessType>>::template independent<0>(DimensionlessType(β1 + dβ)));
				auto const f0 = f;
				f = fRoot.template get<0>().value();
				β1 += dβ;
				while (abs(f) > abs(f0)) {
					β1 -= dβ;
					dβ *= 0.5;
					β1 += dβ;
					fRoot = F(FwdAutoDiff<DimensionlessType, NDIM, std::tuple<DimensionlessType>>::template independent<0>(DimensionlessType(β1)));
					f = fRoot.template get<0>().value();
				}
			}
			ρ = MassDensityType(ρ_);
			β = c * S / EnergyDensityType(W);
			ε = SpecificEnergyType(ε_);
			γ = sqrt(1.0 / (1.0 - β1 * β1));
			if (!useEntropy) {
				if (gasPrim.dualEnergySwitch(eos) < Φo) {
					useEntropy = true;
					goto REDO_PRIMITIVE;
				}
			}
		}
		return gasPrim;
	}
};


#endif /* INCLUDE_SRHD_CONSERVED_HPP_ */
