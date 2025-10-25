/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/
#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "rad_conserved.hpp"
#include "fpe.hpp"

constexpr int ndim = 3;
using Catch::Approx;
static constexpr PhysicalConstants<double> pc { };

int main(int argc, char *argv[]) {
	enableFPE();
	return Catch::Session().run(argc, argv);
}

EquationOfState<double> makeEquationOfState(double Γ, double μ);

struct RadTestCase {
	double rho;
	double kappaA;
	double kappaS;
	double beta[3];   // fluid velocity / c
	double E;         // radiation energy density
	double Ffrac[3];  // flux fractions relative to Er
	double Tg;        // gas temperature for initial ε
	double dt;
	const char *name;
};

TEST_CASE("Implicit energy solve conserves total energy in all regimes", "[implicitEnergySolve]") {
	EquationOfState<double> eos(MolarMassType(1.0));

	// @formatter:off
	std::vector<RadTestCase> cases = {
		{1.0,  1.0,   0.0,     {0.0,0.0,0.0},   1.0,     {0.01,0.0,0.0},                1e6,   1.0, "Thick diffusion"},
		{1.0,  1.0,   0.0,     {0.5,0.0,0.0},   1.0,     {0.01,0.0,0.0},                1e6,   1.0, "Thick diffusion, boosted"},
		{1.0,  0.8,   0.2,     {0.0,0.0,0.0},   1.0,     {0.3,0.0,0.0},                 1e6,   1.0, "Semi-transparent"},
		{1.0,  0.0,   1e-6,    {0.0,0.0,0.0},   1.0,     {0.999,0.0,0.0},               1e6,   1.0, "Free-streaming"},
		{1.0,  1e-3,  0.0,     {0.0,0.0,0.0},   1.0,     {0.9,0.0,0.0},                 1e6,   1.0, "Thin with weak absorption"},
		{1.0,  100.0, 0.0,     {0.0,0.0,0.0},  10.0,     {0.1,0.0,0.0},                 1e6,   1.0, "Stiff absorption dominated"},
		{1.0,  0.01,  100.0,   {0.0,0.0,0.0},   1.0,     {0.6,0.0,0.0},                 1e6,   1.0, "Stiff scattering dominated"},
		{1.0,  1.0,   0.0,     {0.0,0.0,0.0}, 100.0,     {0.1,0.0,0.0},                 1e5,   1.0, "Radiation-dominated gas"},
		{1.0,  1.0,   0.0,     {0.0,0.0,0.0},  0.01,     {0.1,0.0,0.0},                 1e6,   1.0, "Matter-dominated"},
		{1.0,  1.0,   0.0,     {0.95,0.0,0.0},  1.0,     {0.0,0.0,0.0},                 1e6,   1.0, "Ultra-relativistic flow"},
		{1.0,  0.0,   1e-6,    {0.0,0.0,0.0},   1.0,     {+1.0,0.0,0.0},                1e6,   1.0, "Counter-beam +x"},
		{1.0,  0.0,   1e-6,    {0.0,0.0,0.0},   1.0,     {-1.0,0.0,0.0},                1e6,   1.0, "Counter-beam -x"},
		{1.0,  0.5,   0.5,     {0.3,0.0,0.0},   1.0,     {0.7/sqrt(2),0.7/sqrt(2),0.0},1e6,   1.0, "Oblique anisotropy"},
		{1.0,  0.05,  0.05,    {0.7,0.0,0.0},   1.0,     {0.999,0.0,0.0},               1e6,   1.0, "Flux cap edge"},
		{1.0,  1.0,   0.0,     {0.0,0.0,0.0}, 1e-12,     {0.0,0.0,0.0},                 1e7,   1.0, "Emission-only"},
		{1.0,  1.0,   0.0,     {0.0,0.0,0.0},   1.0,     {0.2,0.0,0.0},                 1e3,   1.0, "Absorption-only"},
		{1.0,  0.0,   1.0,     {0.8,0.0,0.0},   1.0,     {0.6,0.0,0.0},                 1e6,   1.0, "Scattering-only"}
	};
		// @formatter:on

	for (size_t i = 0; i < cases.size(); i++) {
		auto const &tc = cases[i];
		Opacity<double> opac;
		GasPrimitive<double, ndim> gasPrim;
		RadConserved<double, ndim> radCon;
		gasPrim.setMassDensity(MassDensityType<double>(tc.rho));
		auto const ε = eos.temperature2energy(gasPrim.getMassDensity(), TemperatureType<double>(tc.Tg));
		gasPrim.setSpecificEnergy(ε);
		gasPrim.setVelocity(
				{ pc.c * DimensionlessType<double>(tc.beta[0]), pc.c * DimensionlessType<double>(tc.beta[1]), pc.c * DimensionlessType<double>(tc.beta[2]) });
		auto const Er = temperature2radiationEnergy<double>(TemperatureType<double>(tc.E));
		radCon.setEnergy(Er);
		radCon.setFlux(
				{ pc.c * DimensionlessType<double>(tc.Ffrac[0]) * Er, pc.c * DimensionlessType<double>(tc.Ffrac[1]) * Er, pc.c
						* DimensionlessType<double>(tc.Ffrac[2]) * Er });
		opac.κₐ = SpecificAreaType(tc.kappaA);
		opac.κₛ = SpecificAreaType(tc.kappaS);
		auto const radCon0 = radCon;
		auto gasCon = gasPrim.toConserved(eos);
		radCon = radCon.implicitRadiationSolve(gasCon, opac, eos, TimeType<double>(tc.dt));
		auto const deltaE = radCon0.getEnergy() - radCon.getEnergy();
		double const Eg0 = (gasPrim.getMassDensity() * gasPrim.getSpecificEnergy()).value();
		double const Er0 = radCon0.getEnergy().value();
		double const Eg1 = Eg0 - deltaE.value();
		double const Er1 = Er0 + deltaE.value();
		INFO("Case: " << tc.name);
		INFO(
				"rho=" << tc.rho << " kappaA=" << tc.kappaA << " kappaS=" << tc.kappaS << " beta=(" << tc.beta[0] << "," << tc.beta[1] << "," << tc.beta[2] << ")" << " E=" << tc.E << " dt=" << tc.dt);

		REQUIRE(Eg0 + Er0 == Approx(Eg1 + Er1).epsilon(1e-10));
	}
}
