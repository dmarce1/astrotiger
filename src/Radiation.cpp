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

inline double lorentzFactor(GasPrimitive const &gasPrim) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	GAS_PRIMITIVE_PROPERTIES(gasPrim);
#pragma GCC diagnostic pop
	auto const β = v * ic;
	auto const β2 = vectorDotProduct(β, β);
	auto const γ = sqrt(1 / (1 - β2));
	return double(γ);
}

inline auto velocity(GasPrimitive const &gasPrim) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	GAS_PRIMITIVE_PROPERTIES(gasPrim);
#pragma GCC diagnostic pop
	auto const γ = lorentzFactor(gasPrim);
	return γ * space2SpaceTime<VelocityType>(c, v);
}

inline TensorType<double> boostTensor(GasPrimitive const &gasPrim, int sign = +1) {
	TensorType<double> Λ;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	GAS_PRIMITIVE_PROPERTIES(gasPrim);
#pragma GCC diagnostic pop
	auto const β = ColumnVector<double, ndim>(v * ic);
	double const γ = lorentzFactor(gasPrim);
	auto const Λij = SquareMatrix<double, ndim>(δ + (sqr(γ) / (γ + 1)) * vectorDyadicProduct(β, β));
	auto const Λoi = double(sign) * γ * β;
	double const Λoo = γ;
	return space2SpaceTime<double>(Λoo, Λoi, Λij);
}

inline auto radiationTensor(RadConserved const &rad) {
	RADIATION_CONSERVED_PROPERTIES(rad);
	auto const P = pressureTensor(rad);
	return space2SpaceTime<typename std::remove_cvref<decltype(Er)>::type>(Er, F, P);
}

GasPrimitive GasConserved::toPrimitive(EquationOfState const &eos) const {
	constexpr double toler = 2 * eps;
	constexpr auto ξo = EnergyDensityType(1e-3);
	static int k = 0;
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
	auto const S1 = sqrt(S2);
	auto const ξi = τ - S2 / (2 * D);
	auto const D2 = sqr(D);
	if (ξi < ξo) {
		auto const β2 = S2 / (S2 + c2 * D2);
		auto const iγ2 = 1.0 - β2;
		auto const iγ = sqrt(iγ2);
		auto const iD = 1 / D;
		ρ = D * iγ;
		v = S * iγ * iD;
		ϰ = K * iD;
		ε = eos.entropy2energy(ρ, ϰ);
	} else if (abs(S1.value()) < tiny) {
		auto const iρ = 1 / D;
		ρ = D;
		ε = τ * iρ;
		v = S * iρ;
		ϰ = K * iρ;
	} else {
		FwdAutoDiff<EnergyDensityType> W;
		FwdAutoDiff<MassDensityType> ρ_;
		FwdAutoDiff<SpecificEnergyType> ε_;
		FwdAutoDiff<DimensionlessType> ϰ_;
		auto const F = [this, Γ, eos, S1, S2, &W, &ϰ_, &ε_, &ρ_](double β) {
			using AutoType = FwdAutoDiff<double>;
			constexpr FwdAutoDiff<double> one = 1.0;
			auto const β_ = AutoType::independentVariable(β);
			auto const p = (FwdAutoDiff(c * S1) / β_ - FwdAutoDiff(c2 * D + τ));
			auto const β2 = sqr(β_);
			W = FwdAutoDiff(c2 * D + τ) + p;
			auto const γ = sqrt(one / (one - β2));
			auto const iγ = one / γ;
			ρ_ = FwdAutoDiff(D) * iγ;
			auto const γ2 = sqr(γ);
			auto const iγ2 = one / γ2;
			auto const iρ = one / ρ_;
			auto const h = W * iγ2 * iρ;
			ε_ = h - c2 - p * iρ;
			return (p - FwdAutoDiff(eos.pressure(ρ_, ε_))) / W;
		};
		double β = tanh(double(c * S1 / (τ + c2 * D)));
		double βmin = β * eps;
		double βmax = 1.0 - eps;
		β = std::min(β, βmax);
		β = std::max(β, βmin);
		auto fRoot = F(β);
		double f = fRoot.D(0).value();
		for (int j = 0; abs(f) > toler; j++, k++) {
			printf("%e %i\n", abs(f), j);
			double const dfdβ = fRoot.D(1).value();
			double const d2fdβ2 = fRoot.D(2).value();
			double const den = 2 * sqr(dfdβ) - f * d2fdβ2;
			double const num = -2 * f * dfdβ;
			double dβHalley = num / (den + copysign(tiny, den));
			double dβNewton = -f / dfdβ;
			βmin = β * eps;
			βmax = 0.99 + (0.015 - 0.005 * β) * β;
			dβNewton = max(min(β + dβNewton, βmax), βmin) - β;
			dβHalley = max(min(β + dβHalley, βmax), βmin) - β;
			auto const fNewtonRoot = F(β + dβNewton);
			auto const fHalleyRoot = F(β + dβHalley);
			double const fNewton = fNewtonRoot.D(0).value();
			double const fHalley = fHalleyRoot.D(0).value();
			double const fLast = f;
			double dβ;
			if (abs(fHalley) < abs(fNewton)) {
				fRoot = fHalleyRoot;
				f = fHalley;
				dβ = dβHalley;
			} else {
				fRoot = fNewtonRoot;
				f = fNewton;
				dβ = dβNewton;
			}
			while (abs(f) > abs(fLast)) {
				dβ *= 0.5;
				auto const fRoot = F(β + dβ);
				f = fRoot.D(0).value();
			}
			β += dβ;
		}
		ρ = MassDensityType(ρ_);
		v = c2 * S / EnergyDensityType(W);
		ε = SpecificEnergyType(ε_);
		ϰ = eos.energy2entropy(ρ, ε);
	}
	return gasPrim;
}

void GasPrimitive::updateEntropy(EquationOfState const &eos) {
	ϰ = eos.energy2entropy(ρ, ε);
}

GasConserved GasPrimitive::toConserved(EquationOfState const &eos) const {
	GasConserved con;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	EOS_PROPERTIES(eos)
	GAS_CONSERVED_PROPERTIES(con);
#pragma GCC diagnostic pop
	auto const v2 = vectorDotProduct(v, v);
	auto const β2 = double(ic2 * v2);
	auto const γ2 = 1 / (1 - β2);
	auto const γ = sqrt(γ2);
	auto const iρ = 1 / ρ;
	auto const p = eos.pressure(ρ, ε);
	auto const h = c2 + ε + p * iρ;
	S = ρ * γ2 * h * ic2 * v;
	D = ρ * γ;
	τ = ρ * γ2 * ε + (γ2 - 1) * p + (D * c2 * (γ - 1));
	K = D * ϰ;
	return con;
}

void implicitEnergySolve(GasPrimitive &prim, RadConserved &rad, Opacity const &opac, EquationOfState const &eos, cgs::Seconds const &dt) {
	using std::abs;
	using std::copysign;
	using std::max;
	constexpr double toler = 16 * eps;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	RADIATION_CONSERVED_PROPERTIES(rad);
	EOS_PROPERTIES(eos);
	OPACITY_PROPERTIES(opac);
	GAS_PRIMITIVE_PROPERTIES(prim);
#pragma GCC diagnostic pop
	auto const ρ_ = FwdAutoDiff(ρ);
	auto const Eg0 = FwdAutoDiff(ρ * ε);
	auto const Er0 = FwdAutoDiff(Er);
	auto const Etot0 = Er0 + Eg0;
	auto const λ = FwdAutoDiff(ρ * c * κₐ * dt);
	auto const Fgas = [eos, Etot0, Er0, λ, ρ_](SpecificEnergyType ε_) {
		constexpr auto one = FwdAutoDiff(DimensionlessType(1.0));
		auto const ε = FwdAutoDiff<SpecificEnergyType>::independentVariable(ε_);
		auto const Tg = gasEnergy2temperature(eos, ε);
		auto const Bp = FwdAutoDiff(aR) * sqr(sqr(Tg));
		auto const Er = Etot0 - ρ_ * ε;
		auto const f = ((one + λ) * Er - (Er0 + λ * Bp)) / Er0;
		return f;
	};
	auto const Frad = [eos, Etot0, Er0, Eg0, λ, ρ_](EnergyDensityType Er_) {
		constexpr auto one = FwdAutoDiff(DimensionlessType(1.0));
		auto const Er = FwdAutoDiff<EnergyDensityType>::independentVariable(Er_);
		auto const ε = (Etot0 - Er) / ρ_;
		auto const Tg = gasEnergy2temperature(eos, ε);
		auto const Bp = FwdAutoDiff(aR) * sqr(sqr(Tg));
		auto const f = (λ * Bp - ((one + λ) * Er - Er0)) / Eg0;
		return f;
	};
	double error = 1.0;
	auto Tr = radiationEnergy2temperature(Er);
	auto Tg = gasEnergy2temperature(eos, ε);
	printf("Tr0 = %e | Tg0 = %e   %e  %e\n", (double) Tr.value(), (double) Tg.value(), (double) (ρ * ε).value(), (double) (Er).value());
	for (int j = 0; error > toler; j++) {
		bool gasForm;
		double e;
		FwdAutoDiff<DimensionlessType> fRoot;
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
		if (gasForm) {
			ε += SpecificEnergyType(de.value());
			Er = (Etot0 - ρ * ε).D(0);
			e = ε.value();
		} else {
			Er += EnergyDensityType(de.value());
			ε = ((Etot0 - Er) / ρ_).D(0);
			e = Er.value();
		}
		error = abs((de / e).value());
		printf("Er = %e\n", Er.value());
		Tr = radiationEnergy2temperature(Er);
		Tg = gasEnergy2temperature(eos, ε);
		printf("%c: %i %e | Tr  = %e | Tg  = %e | %e  %e \n", gasForm ? 'G' : 'R', j, (double) error, (double) Tr.value(), (double) Tg.value(), (ρ * ε).value(),
				(Er).value());
	}
}

void implicitRadiationSolve(GasConserved &gasCon, RadConserved &radCon, Opacity const &opac, EquationOfState const &eos, cgs::Seconds const &dt) {
	constexpr double toler = 4 * eps;
	GasPrimitive gasPrim;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	GAS_CONSERVED_PROPERTIES(gasCon);
	RADIATION_CONSERVED_PROPERTIES(radCon);
	GAS_PRIMITIVE_PROPERTIES(gasPrim);
	EOS_PROPERTIES(eos)
	OPACITY_PROPERTIES(opac)
	auto const norm = sqrt(sqr(vectorMagnitude(F)) + sqr(Er) + sqr(vectorMagnitude(c * S)) + sqr(τ)).value();
#pragma GCC diagnostic pop
	double error = 1.0;
	auto const F0 = F;
	auto const E0 = Er;
	auto const S0 = S;
	auto const τ0 = τ;
	for (int j = 0; error > toler; j++) {
		auto const F1 = F;
		auto const E1 = Er;
		gasPrim = gasCon.toPrimitive(eos);
		auto const v2 = vectorDotProduct(v, v);
		auto const β = sqrt(ic2 * v2);
		auto const λ = ρ * c * (κₐ + κₛ) * dt;
		auto R = radiationTensor(radCon);
		auto Λ = boostTensor(gasPrim, -1);
		auto iΛ = boostTensor(gasPrim, +1);
		R = Λ * R * Λ;
		Er = spaceTime2Scalar(R);
		F = spaceTime2Vector(R);
		implicitEnergySolve(gasPrim, radCon, opac, eos, dt);
		F = F / (1 + λ);
		R = radiationTensor(radCon);
		R = iΛ * R * iΛ;
		Er = spaceTime2Scalar(R);
		F = spaceTime2Vector(R);
		auto const dF = F - F0;
		auto const dE = Er - E0;
		S = S0 - ic * dF;
		τ = τ0 - dE;
		F = F0 + dF;
		Er = E0 + dE;
		error = double(sqrt(sqr(vectorMagnitude(F - F1)) + sqr(Er - E1)).value() / norm);
		printf("%i %e %e %e %e\n", j, (double) β.value(), (double) vectorMagnitude(F).value(), (double) λ, (double) error);
	}
}

}

void radiation_test() {
}
