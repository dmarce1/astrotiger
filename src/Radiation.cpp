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

inline TensorType<double> boostTensor(GasPrimitive const &gasPrim) {
	TensorType<double> Λ;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	GAS_PRIMITIVE_PROPERTIES(gasPrim);
#pragma GCC diagnostic pop
	auto const β = ColumnVector<double, ndim>(v * ic);
	double const γ = lorentzFactor(gasPrim);
	auto const Λij = SquareMatrix<double, ndim>(δ + (sqr(γ) / (γ + 1)) * vectorDyadicProduct(β, β));
	auto const Λoi = γ * β;
	double const Λoo = γ;
	return space2SpaceTime<double>(Λoo, Λoi, Λij);
}

inline auto radiationTensor(RadConserved const &rad) {
	RADIATION_CONSERVED_PROPERTIES(rad);
	auto const P = pressureTensor(rad);
	return space2SpaceTime<typename std::remove_cvref<decltype(Er)>::type>(Er, F, P);
}

GasPrimitive gasCon2Prim(GasConserved const &con, Material const &mat) {
	constexpr double toler = eps;
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
	GAS_CONSERVED_PROPERTIES(con);
	MATERIAL_PROPERTIES(mat)
	GAS_PRIMITIVE_PROPERTIES(gasPrim);
#pragma GCC diagnostic pop
	auto const S2 = vectorDotProduct(S, S);
	auto const S1 = sqrt(S2);
	if (S2.value() == 0.0) {
		auto const iρ = 1 / D;
		ρ = D;
		ε = τ * iρ;
		v = S * iρ;
	} else {
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
			auto const f = (p - FwdAutoDiff(Γ - 1.0) * _ρ * _ε) / W;
			return f;
		};
		double β = tanh(double(c * S1 / (τ + c2 * D)));
		double βmin = β * eps;
		double βmax = 1.0 - eps;
		β = std::min(β, βmax);
		β = std::max(β, βmin);
		auto fRoot = F(β);
		double f = fRoot.D(0).value();
		for (int j = 0; abs(f) > toler; j++, k++) {
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
		ρ = MassDensityType(_ρ);
		ε = SpecificEnergyType(_ε);
		v = c2 * S / EnergyDensityType(W);
	}
	return gasPrim;
}

GasConserved gasPrim2Con(GasPrimitive const &gasPrim, Material const &mat) {
	GasConserved con;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	GAS_PRIMITIVE_PROPERTIES(gasPrim);
	MATERIAL_PROPERTIES(mat)
	GAS_CONSERVED_PROPERTIES(con);
#pragma GCC diagnostic pop
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
// μ (Γ - 1) ε = p = kB nₐ T
template<typename Type>
auto gasEnergy2temperature(Material const &mat, Type ε) {
	MATERIAL_PROPERTIES(mat);
	return ε * ((μ * (Γ - 1)) / (kB * nₐ));
}

FwdAutoDiff<TemperatureType> gasEnergy2temperature(Material const &mat, FwdAutoDiff<SpecificEnergyType> const &ε) {
	MATERIAL_PROPERTIES(mat);
	return ε * FwdAutoDiff(gasEnergy2temperature(mat, ε.D(0)) / ε.D(0));
}

template<typename Type>
auto radiationEnergy2temperature(Type const &Er) {
	return sqrt(sqrt(Er / aR));
}

void implicitEnergySolve(GasPrimitive &prim, RadConserved &rad, Opacity const &opac, Material const &mat, cgs::Seconds const &dt) {
	using std::abs;
	using std::copysign;
	using std::max;
	constexpr double toler = 16 * eps;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	RADIATION_CONSERVED_PROPERTIES(rad);
	MATERIAL_PROPERTIES(mat);
	OPACITY_PROPERTIES(opac);
	GAS_PRIMITIVE_PROPERTIES(prim);
#pragma GCC diagnostic pop
	auto const ρ_ = FwdAutoDiff(ρ);
	auto const Eg0 = FwdAutoDiff(ρ * ε);
	auto const Er0 = FwdAutoDiff(Er);
	auto const Etot0 = Er0 + Eg0;
	auto const λ = FwdAutoDiff(ρ * c * κₐ * dt);
	auto const Fgas = [mat, Etot0, Er0, λ, ρ_](SpecificEnergyType ε_) {
		constexpr auto one = FwdAutoDiff(DimensionlessType(1.0));
		auto const ε = FwdAutoDiff<SpecificEnergyType>::independentVariable(ε_);
		auto const Tg = gasEnergy2temperature(mat, ε);
		auto const Bp = FwdAutoDiff(aR) * sqr(sqr(Tg));
		auto const Er = Etot0 - ρ_ * ε;
		auto const f = ((one + λ) * Er - (Er0 + λ * Bp)) / Er0;
		return f;
	};
	auto const Frad = [mat, Etot0, Er0, Eg0, λ, ρ_](EnergyDensityType Er_) {
		constexpr auto one = FwdAutoDiff(DimensionlessType(1.0));
		auto const Er = FwdAutoDiff<EnergyDensityType>::independentVariable(Er_);
		auto const ε = (Etot0 - Er) / ρ_;
		auto const Tg = gasEnergy2temperature(mat, ε);
		auto const Bp = FwdAutoDiff(aR) * sqr(sqr(Tg));
		auto const f = (λ * Bp - ((one + λ) * Er - Er0)) / Eg0;
		return f;
	};
	double error = 1.0;
	auto Tr = radiationEnergy2temperature(Er);
	auto Tg = gasEnergy2temperature(mat, ε);
	printf("Tr0 = %e | Tg0 = %e   %e  %e\n", Tr, Tg, (ρ * ε).value(), (Er).value());
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
		Tr = radiationEnergy2temperature(Er);
		Tg = gasEnergy2temperature(mat, ε);
		printf("%c: %i %e | Tr  = %e | Tg  = %e | %e  %e \n", gasForm ? 'G' : 'R', j, error, Tr, Tg, (ρ * ε).value(), (Er).value());
	}
}

void implicitRadiationSolve(GasConserved &gas, RadConserved &rad, Opacity const &opac, Material const &mat, cgs::Seconds const &dt) {
	GasPrimitive gasPrim;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	GAS_CONSERVED_PROPERTIES(gas);
	RADIATION_CONSERVED_PROPERTIES(rad);
	GAS_PRIMITIVE_PROPERTIES(gasPrim);
	MATERIAL_PROPERTIES(mat)
	OPACITY_PROPERTIES(opac)
#pragma GCC diagnostic pop
	gasPrim = gasCon2Prim(gas, mat);
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
	auto const iΛ = boostTensor(gasPrim);
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

void radiation_test() {
}
