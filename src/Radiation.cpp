#include "Radiation.hpp"

namespace Radiation {

#define XXX(name, type, arg, dim) copy_const_t<decltype(arg), std::conditional_t<dim == 1, type, ColumnVector<type, dim>>>& name = arg.name;

inline auto pressureTensor(RadConserved const &rad) {
	using std::max;
	RADIATION_CONSERVED_PROPERTIES(rad);
	auto const magF = sqrt(F | F);
	auto const imagF = 1 / max(magF, EnergyFluxType(tiny));
	auto const iEr = 1 / Er;
	auto const f = magF * ic * iEr;
	auto const n = F * imagF;
	auto const f2 = sqr(f);
	auto const ξ = (3 + 4 * f2) / (5 + 2 * sqrt(4 - 3 * f2));
	auto const Ωdif = 0.5 * (1 - ξ);
	auto const Ωstr = 0.5 * (3 * ξ - 1);
	auto const P = Er * (Ωdif * δ < ndim > +Ωstr * (n ^ n));
	return P;
}

inline auto velocity(GasPrimitive const &gasPrim) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	GAS_PRIMITIVE_PROPERTIES(gasPrim);
#pragma GCC diagnostic pop
	return c * γ * space2SpaceTime<double>(1, β);
}

inline auto radiationTensor(RadConserved const &rad) {
	RADIATION_CONSERVED_PROPERTIES(rad);
	auto const P = pressureTensor(rad);
	return space2SpaceTime<typename std::remove_cvref<decltype(Er)>::type>(Er, ic * F, P);
}

template<typename T>
using AD = FwdAutoDiff<T, 2, DimensionlessType>;

auto GasPrimitive::dualEnergySwitch(EquationOfState const &eos) const {
	auto const h = c2 + ε + eos.energy2pressure(ρ, ε) / ρ;
	auto const num = ρ * ε;
	auto const den = ρ * γ * (γ * h - c2);
	return num / den;
}

GasPrimitive GasConserved::toPrimitive(EquationOfState const &eos) const {
	constexpr double toler = 2 * eps;
	constexpr DimensionlessType Φo(1e-3);
	using std::abs;
	using std::copysign;
	using std::ignore;
	using std::max;
	using std::min;
	using std::tanh;
	GasPrimitive gasPrim;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	GAS_PRIMITIVE_PROPERTIES(gasPrim);
#pragma GCC diagnostic pop
	auto const S2 = (S | S);
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
		printf("<<<\n");
		AD<EnergyDensityType> W;
		AD<MassDensityType> ρ_;
		AD<SpecificEnergyType> ε_;
		AD<DimensionlessType> ϰ_;
		using Function = std::function<AD<DimensionlessType>(AD<DimensionlessType>)>;
		Function const fEnergy = [this, eos, S1, S2, τo, &W, &ϰ_, &ε_, &ρ_](AD<DimensionlessType> β) {
			auto const one = makeConstantLike(β, DimensionlessType(1));
			auto const p = (AD<EnergyDensityType>(c * S1) / β - AD<EnergyDensityType>(c2 * D + max(τ, τo)));
			AD<DimensionlessType> const β2 = β * β;
			W = AD<EnergyDensityType>(c2 * D + max(τ, τo)) + p;
			AD<DimensionlessType> const γ = sqrt(one / (one - β2));
			auto const iγ = one / γ;
			ρ_ = AD(D) * iγ;
			auto const γ2 = sqr(γ);
			auto const iγ2 = one / γ2;
			auto const iρ = one / ρ_;
			auto const h = W * iγ2 * iρ;
			ε_ = h - c2 - p * iρ;
			auto const peos = AD(eos.energy2pressure(ρ_, ε_));
			auto const rc = AD<DimensionlessType>((p - peos) / W);
			return rc;
		};
		//ϰ D = P / ρ^(Γ-1)  g^0 cm^2 s^-2 cm^-3 g^(Γ-1) cm^(1 - Γ)
		Function const fEntropy = [this, eos, S2, D2](AD<DimensionlessType> β) {
			AD<DimensionlessType> const one(DimensionlessType(1.0));
			auto const β2 = β;
			auto const γ = one / sqrt(one - β2);
			auto const ρ = D / γ;
			auto const ϰ(K / D);
			auto const Tg(eos.entropy2temperature(ρ, ϰ));
			auto const e = eos.temperature2energy(ρ, Tg);
			auto const ε = AD<SpecificEnergyType>(e);
			auto const p = eos.energy2pressure(ρ, ε);
			auto const h = ε + p / ρ;
			auto const h2 = h * h;
			auto const γ2 = γ * γ;
			auto const rc = AD<DimensionlessType>(γ2 - c2 * S2 / (D2 * h2) - one);
			return rc;
		};
		bool useEntropy = false;
		REDO_PRIMITIVE: double β1, βmin, βmax, f;
		Function F;
		AD<DimensionlessType> fRoot;
		β1 = tanh(double(c * S1 / (τo + c2 * D)));
		βmin = β1 * eps;
		βmax = 0.9;
		β1 = min(β1, βmax);
		β1 = max(β1, βmin);
		F = useEntropy ? fEntropy : fEnergy;
		fRoot = F(AD<DimensionlessType>::independentVariable(DimensionlessType(β1)));
		f = fRoot[0].value();
		printf("%e %e %e\n", fRoot[0].value(), fRoot[1].value(), fRoot[2].value());
		for (int j = 0; abs(f) > toler; j++) {
			double const dfdβ = fRoot[1].value();
			double dβ = -f / dfdβ;
			βmin = 0.5 * β1;
			βmax = 0.5 + 0.5 * β1;
			dβ = max(min(β1 + dβ, βmax), βmin) - β1;
			fRoot = F(AD<DimensionlessType>::independentVariable(DimensionlessType(β1 + dβ)));
			auto const f0 = f;
			f = fRoot[0].value();
			β1 += dβ;
			while (abs(f) > abs(f0)) {
				β1 -= dβ; //	printf("Tr0 = %e | Tg0 = %e   %e  %e\n", (double) Tr.value(), (double) Tg.value(), (double) (ρ * ε).value(), (double) (Er).value());

				dβ *= 0.5;
				β1 += dβ;
				fRoot = F(AD<DimensionlessType>::independentVariable(DimensionlessType(β1)));
				f = fRoot[0].value();
			}
		}
		ρ = MassDensityType(ρ_);
		β = c * S / EnergyDensityType(W);
		ε = SpecificEnergyType(ε_);
		γ = sqrt(1.0 / (1.0 - β1 * β1));
		if (!useEntropy) {
			if (gasPrim.dualEnergySwitch(eos) < Φo) {
				useEntropy = true;
				printf("ENTROPY\n");
				goto REDO_PRIMITIVE;
			}
		}
		printf(">>>\n");
	}
	return gasPrim;
}

void GasPrimitive::updateLorentz() {
	auto const β2 = (β | β);
	γ = DimensionlessType(1.0) / sqrt(DimensionlessType(1.0) - β2);
}

GasConserved GasPrimitive::toConserved(EquationOfState const &eos) const {
	GasConserved con;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	GAS_CONSERVED_PROPERTIES(con);
#pragma GCC diagnostic pop
	auto const β2 = double(β | β);
	auto const γ2 = 1 / (1 - β2);
	auto const γ = sqrt(γ2);
	auto const iρ = 1 / ρ;
	auto const p = eos.energy2pressure(ρ, ε);
	auto const T = eos.energy2temperature(ρ, ε);
	auto const ϰ = eos.temperature2entropy(ρ, T);
	auto const h = c2 + ε + p * iρ;
	S = ρ * γ2 * h * β * ic;
	D = ρ * γ;
	τ = ρ * γ2 * ε + (γ2 - 1) * p + c2 * D * (γ - 1);
	K = D * ϰ;
	return con;
}

auto GasPrimitive::temperature(EquationOfState const &eos) const {
	return eos.energy2temperature(ρ, ε);
}

auto RadConserved::temperature() const {
	return pow<Rational(1, 4)>(Er / aR);
}

void implicitEnergySolve(GasPrimitive &prim, RadConserved &rad, Opacity const &opac, EquationOfState const &eos, TimeType const &dt) {
	using std::abs;
	using std::copysign;
	using std::min;
	using std::max;
	constexpr double toler = 16 * eps;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	RADIATION_CONSERVED_PROPERTIES(rad);
	OPACITY_PROPERTIES(opac);
	GAS_PRIMITIVE_PROPERTIES(prim);
#pragma GCC diagnostic pop
	constexpr auto one = AD<DimensionlessType>(DimensionlessType(1.0));
	auto const ρ_ = AD<MassDensityType>(ρ);
	auto const Eg0 = AD<EnergyDensityType>(ρ * ε);
	auto const Er0 = AD<EnergyDensityType>(Er);
	auto const Etot0 = AD<EnergyDensityType>(Er0 + Eg0);
	auto const λ = AD<DimensionlessType>(ρ * c * κₐ * dt);
	auto const Fgas = [eos, prim, Etot0, Er0, Eg0, λ, ρ_, one](SpecificEnergyType ε_) {
		auto const ε = AD<SpecificEnergyType>::independentVariable(ε_);
		auto const Tg = prim.temperature(eos);
		auto const Tg2 = Tg * Tg;
		auto const Tg4 = Tg2 * Tg2;
		auto const Bp = aR * Tg4;
		auto const Er = Etot0 - ρ_ * ε;
		auto const f = ((one + λ) * Er - (Er0 + λ * Bp)) / (Er0 + Eg0);
		return f;
	};
	auto const Frad = [eos, prim, Etot0, Er0, Eg0, λ, ρ_, one](EnergyDensityType Er_) {
		auto const Er = AD<EnergyDensityType>::independentVariable(Er_);
		auto const Tg = prim.temperature(eos);
		auto const Bp = aR * sqr(sqr(Tg));
		auto const f = (λ * Bp - ((one + λ) * Er - Er0)) / (Er0 + Eg0);
		return f;
	};
	double error = 1.0;
	for (int j = 0; error > toler; j++) {
		bool gasForm;
		double e;
		AD<DimensionlessType> fRoot;
		if ((ρ * ε) < Er) {
			gasForm = true;
			fRoot = Fgas(ε);
		} else {
			gasForm = false;
			fRoot = Frad(Er);
		}
		auto const f = fRoot[0];
		auto const dfde = fRoot[1];
		auto const d2fde2 = fRoot[2];
		auto const den = double(max(sqr(dfde), 2 * sqr(dfde) - f * d2fde2));
		auto const num = -2 * f * dfde;
		auto de = num / (den + copysign(tiny, den));
		auto ε1 = ε;
		auto Er1 = Er;
		if (gasForm) {
			ε1 += SpecificEnergyType(de.value());
			Er1 = (Etot0 - ρ * ε1)[0];
		} else {
			Er1 += EnergyDensityType(de.value());
			ε1 = ((Etot0 - Er) / ρ_)[0];
		}
		ε = max(0.5 * ε, ε1);
		Er = max(0.5 * Er, Er1);
		e = gasForm ? ε.value() : Er.value();
		error = abs((de / e).value());
	}
}

void implicitRadiationSolve(GasConserved &gasCon, RadConserved &radCon, Opacity const &opac, EquationOfState const &eos, TimeType const &dt,
		DimensionlessType ξo) {
	using std::max;
	using std::min;
	enableFPE();
	constexpr double toler = 4 * eps;
	GasPrimitive gasPrim;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	GAS_CONSERVED_PROPERTIES(gasCon);
	RADIATION_CONSERVED_PROPERTIES(radCon);
	GAS_PRIMITIVE_PROPERTIES(gasPrim);
	OPACITY_PROPERTIES(opac)
	auto const norm = sqrt(sqr(sqrt(F | F)) * ic2 + sqr(Er) + c2 * sqr(sqrt(S | S)) + sqr(τ)).value();
#pragma GCC diagnostic pop
	double error = 1.0;
	auto const F0 = F;
	auto const E0 = Er;
	auto const S0 = S;
	auto const τ0 = τ;
	auto R0 = radiationTensor(radCon);
	for (int j = 0; error > toler; j++) {
		auto const F1 = F;
		auto const E1 = Er;
		gasPrim = gasCon.toPrimitive(eos);
		auto const Λ = space2SpaceTime<DimensionlessType>(γ, γ * β, δ < ndim > +(sqr(γ) / (γ + 1)) * (β ^ β));
		auto const dt0 = dt / γ;
		auto const iΛ = space2SpaceTime<DimensionlessType>(γ, -γ * β, δ < ndim > +(sqr(γ) / (γ + 1)) * (β ^ β));
		auto const λ = ρ * c * (κₐ + κₛ) * dt0;
		auto R = Λ * R0 * Λ;
		Er = spaceTime2Scalar(R);
		F = c * spaceTime2Vector(R);
		implicitEnergySolve(gasPrim, radCon, opac, eos, dt0);
		F = F / (1 + λ);
		R = radiationTensor(radCon);
		R = iΛ * R * iΛ;
		Er = spaceTime2Scalar(R);
		F = c * spaceTime2Vector(R);
		S = S0 + (F - F0) * ic2;
		auto const τo = c2 * (γ - 1.0) * D;
		auto const τmin = 0.5 * (τo + τ);
		τ = τ0 + E0 - Er;
		τ = max(τmin, τ);
		auto const T = eos.energy2temperature(ρ, ε);
		K = D * eos.temperature2entropy(ρ, T);
		error = double(sqrt(sqr(sqrt((F - F1) | (F - F1))) * ic2 + sqr(Er - E1)).value() / norm);
	}
}

}

#include <iostream>
#include <vector>
#include <cmath>
using namespace Radiation;

// Helper to compare two primitive states
void compare(GasPrimitive const &dustPrim, GasPrimitive const &fullPrim, DimensionlessType xi, DimensionlessType xi0, int caseId) {
	auto relErr = [](auto a, auto b) {
		double denom = std::max( { 1e-14, std::abs(b) });
		return std::abs(a - b) / denom;
	};

	double errRho = relErr(dustPrim.ρ.value(), fullPrim.ρ.value());
	double errEps = relErr(dustPrim.ε.value(), fullPrim.ε.value());
	double errGam = relErr(dustPrim.γ.value(), fullPrim.γ.value());

	std::cout << "Case " << caseId << " xi=" << xi.value() << " xi0=" << xi0.value() << "  ρ_err=" << errRho << "  ε_err=" << errEps << "  γ_err=" << errGam
			<< "\n";
}

void radiation_test() {
}

