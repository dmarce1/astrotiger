#include "Radiation.hpp"

namespace Radiation {

#define XXX(name, type, arg, dim) copy_const_t<decltype(arg), std::conditional_t<dim == 1, type, ColumnVector<type, dim>>>& name = arg.name;

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

inline auto velocity(GasPrimitive const &gasPrim) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"v
	GAS_PRIMITIVE_PROPERTIES(gasPrim);
#pragma GCC diagnostic pop
	return c * γ * space2SpaceTime<double>(1, β);
}

inline auto radiationTensor(RadConserved const &rad) {
	RADIATION_CONSERVED_PROPERTIES(rad);
	auto const P = pressureTensor(rad);
	return space2SpaceTime<typename std::remove_cvref<decltype(Er)>::type>(Er, F, P);
}

template<typename T>
using AD = FwdAutoDiff<T, 2, 1>;

GasPrimitive GasConserved::toPrimitive(EquationOfState const &eos, DimensionlessType ξo) const {
	constexpr double toler = 2 * eps;
	using std::abs;
	using std::copysign;
	using std::ignore;
	using std::max;
	using std::min;
	using std::tanh;
	GasPrimitive gasPrim;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	EOS_PROPERTIES(eos)
	GAS_PRIMITIVE_PROPERTIES(gasPrim);
#pragma GCC diagnostic pop
	auto const S2 = vectorDotProduct(S, S);
	auto const D2 = sqr(D);
	auto const S1 = sqrt(S2);
	auto const τo = (sqrt(D2 + S2) - D);
	if (S2.value() < tiny) {
		auto const iρ = c2 / D;
		γ = 1.0;
		ρ = D * ic2;
		ε = max(τ, τo) * iρ;
		β = S * iρ * ic2;
	} else {
		FwdAutoDiff<EnergyDensityType> W;
		FwdAutoDiff<MassDensityType> ρ_;
		FwdAutoDiff<SpecificEnergyType> ε_;
		FwdAutoDiff<DimensionlessType> ϰ_;
		auto const F = [this, Γ, eos, S1, S2, τo, &W, &ϰ_, &ε_, &ρ_](AD<DimensionlessType> β) {
			auto const one = AD<DimensionlessType>(DimensionlessType(1.0));
			auto const p = (AD<EnergyDensityType>(S1) / β - AD<EnergyDensityType>(D + max(τ, τo)));
			AD<DimensionlessType> const β2 = β * β;
			W = AD<EnergyDensityType>(D + max(τ, τo)) + p;
			auto const γ = sqrt(one / (one - β2));
			auto const iγ = one / γ;
			ρ_ = ic2 * FwdAutoDiff(D) * iγ;
			auto const γ2 = sqr(γ);
			auto const iγ2 = one / γ2;
			auto const iρ = one / ρ_;
			auto const h = W * iγ2 * iρ;
			ε_ = h - c2 - p * iρ;
			auto const peos = FwdAutoDiff(eos.pressure(ρ_, ε_));
			return (p - peos) / W;
		};
		double β1 = tanh(double(S1 / (τo + D)));
		double βmin = β1 * eps;
		double βmax = 0.9;
		β1 = std::min(β1, βmax);
		β1 = std::max(β1, βmin);
		auto fRoot = F(AD<DimensionlessType>::independentVariable(DimensionlessType(β1)));
		double f = fRoot.D(0).value();
		printf("%e %e %e\n", fRoot.D(1).value(), fRoot.D(1).value(), fRoot.D(1).value());
		for (int j = 0; abs(f) > toler; j++) {
			double const dfdβ = fRoot.D(1).value();
			double dβ = -f / dfdβ;
			βmin = 0.5 * β1;
			βmax = 0.5 + 0.5 * β1;
			dβ = max(min(β1 + dβ, βmax), βmin) - β1;
			fRoot = F(AD<DimensionlessType>::independentVariable(DimensionlessType(β1 + dβ)));
			auto const f0 = f;
			f = fRoot.D(0).value();
			β1 += dβ;
			while (abs(f) > abs(f0)) {
				β1 -= dβ;
				dβ *= 0.5;
				β1 += dβ;
				fRoot = F(AD<DimensionlessType>::independentVariable(DimensionlessType(β1)));
				f = fRoot.D(0).value();
			}
		}
		ρ = MassDensityType(ρ_);
		β = S / EnergyDensityType(W);
		ε = SpecificEnergyType(ε_);
		γ = sqrt(1.0 / (1.0 - β1 * β1));
	}
	return gasPrim;
}

GasConserved GasPrimitive::toConserved(EquationOfState const &eos) const {
	GasConserved con;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	EOS_PROPERTIES(eos)
	GAS_CONSERVED_PROPERTIES(con);
#pragma GCC diagnostic pop
	auto const β2 = double(vectorDotProduct(β, β));
	auto const γ2 = 1 / (1 - β2);
	auto const γ = sqrt(γ2);
	auto const iρ = 1 / ρ;
	auto const p = eos.pressure(ρ, ε);
	auto const ϰ = eos.entropy(ρ, ε);
	auto const h = c2 + ε + p * iρ;
	S = ρ * γ2 * h * β;
	D = ρ * γ * c2;
	τ = ρ * γ2 * ε + (γ2 - 1) * p + (D * (γ - 1));
	K = D * ϰ;
	return con;
}

void implicitEnergySolve(GasPrimitive &prim, RadConserved &rad, Opacity const &opac, EquationOfState const &eos, cgs::Seconds const &dt) {
	using std::abs;
	using std::copysign;
	using std::min;
	using std::max;
	constexpr double toler = 16 * eps;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	RADIATION_CONSERVED_PROPERTIES(rad);
	EOS_PROPERTIES(eos);
	OPACITY_PROPERTIES(opac);
	GAS_PRIMITIVE_PROPERTIES(prim);
#pragma GCC diagnostic pop
	constexpr auto one = AD<DimensionlessType>(DimensionlessType(1.0));
	auto const ρ_ = AD<MassDensityType>(ρ);
	auto const Eg0 = AD<EnergyDensityType>(ρ * ε);
	auto const Er0 = AD<EnergyDensityType>(Er);
	auto const Etot0 = AD<EnergyDensityType>(Er0 + Eg0);
	auto const λ = AD<DimensionlessType>(ρ * c * κₐ * dt);
	auto const Fgas = [eos, Etot0, Er0, Eg0, λ, ρ_, one](SpecificEnergyType ε_) {
		auto const ε = AD<SpecificEnergyType>::independentVariable(ε_);
		auto const Tg = gasEnergy2temperature(eos, ε);
		auto const Tg2 = Tg * Tg;
		auto const Tg4 = Tg2 * Tg2;
		auto const Bp = aR * Tg4;
		auto const Er = Etot0 - ρ_ * ε;
		auto const f = ((one + λ) * Er - (Er0 + λ * Bp)) / (Er0 + Eg0);
		return f;
	};
	auto const Frad = [eos, Etot0, Er0, Eg0, λ, ρ_, one](EnergyDensityType Er_) {
		auto const Er = FwdAutoDiff<EnergyDensityType>::independentVariable(Er_);
		auto const ε = (Etot0 - Er) / ρ_;
		auto const Tg = gasEnergy2temperature(eos, ε);
		auto const Bp =aR * sqr(sqr(Tg));
		auto const f = (λ * Bp - ((one + λ) * Er - Er0)) / (Er0 + Eg0);
		return f;
	};
	double error = 1.0;
	auto Tr = radiationEnergy2temperature(Er);
	auto Tg = gasEnergy2temperature(eos, ε);
//	printf("Tr0 = %e | Tg0 = %e   %e  %e\n", (double) Tr.value(), (double) Tg.value(), (double) (ρ * ε).value(), (double) (Er).value());
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
		auto const f = fRoot.D(0);
		auto const dfde = fRoot.D(1);
		auto const d2fde2 = fRoot.D(2);
		auto const den = double(max(sqr(dfde), 2 * sqr(dfde) - f * d2fde2));
		auto const num = -2 * f * dfde;
		auto de = num / (den + copysign(tiny, den));
		auto ε1 = ε;
		auto Er1 = Er;
		if (gasForm) {
			ε1 += SpecificEnergyType(de.value());
			Er1 = (Etot0 - ρ * ε1).D(0);
		} else {
			Er1 += EnergyDensityType(de.value());
			ε1 = ((Etot0 - Er) / ρ_).D(0);
		}
		ε = max(0.5 * ε, ε1);
		Er = max(0.5 * Er, Er1);
		e = gasForm ? ε.value() : Er.value();
		error = abs((de / e).value());
		Tr = radiationEnergy2temperature(Er);
		Tg = gasEnergy2temperature(eos, ε);
		printf("%c: %i %e | Tr  = %e | Tg  = %e | %e  %e \n", gasForm ? 'G' : 'R', j, (double) error, (double) Tr.value(), (double) Tg.value(), (ρ * ε).value(),
				(Er).value());
	}
}

void implicitRadiationSolve(GasConserved &gasCon, RadConserved &radCon, Opacity const &opac, EquationOfState const &eos, cgs::Seconds const &dt,
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
	EOS_PROPERTIES(eos)
	OPACITY_PROPERTIES(opac)
	auto const norm = sqrt(sqr(vectorMagnitude(F)) + sqr(Er) + sqr(vectorMagnitude(S)) + sqr(τ)).value();
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
		gasPrim = gasCon.toPrimitive(eos, ξo);
		auto const β2 = vectorDotProduct(β, β);
		auto const Λ = space2SpaceTime<DimensionlessType>(γ, γ * β, δ + (sqr(γ) / (γ + 1)) * vectorDyadicProduct(β, β));
		auto const dt0 = dt / γ;
		auto const iΛ = space2SpaceTime<DimensionlessType>(γ, -γ * β, δ + (sqr(γ) / (γ + 1)) * vectorDyadicProduct(β, β));
		auto const λ = ρ * c * (κₐ + κₛ) * dt0;
		auto R = Λ * R0 * Λ;
		Er = spaceTime2Scalar(R);
		F = spaceTime2Vector(R);
		implicitEnergySolve(gasPrim, radCon, opac, eos, dt0);
		F = F / (1 + λ);
		R = radiationTensor(radCon);
		R = iΛ * R * iΛ;
		Er = spaceTime2Scalar(R);
		F = spaceTime2Vector(R);
		S = S0 + F - F0;
		auto const τo = (γ - 1.0) * D;
		auto const τmin = 0.5 * (τo + τ);
		τ = τ0 + E0 - Er;
		τ = max(τmin, τ);
		K = D * eos.entropy(ρ, ε);
		error = double(sqrt(sqr(vectorMagnitude(F - F1)) + sqr(Er - E1)).value() / norm);
		auto const Tr = radiationEnergy2temperature(Er);
		auto const Tg = gasEnergy2temperature(eos, ε);
		printf("%i %e %e %e %e\n", j, (double) sqrt(β2.value()), (double) Tr.value(), (double) Tg.value(), (double) error);
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

