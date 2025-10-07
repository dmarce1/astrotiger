#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>   // <-- this is the key include
#include "Radiation.hpp"

using Catch::Approx;
using namespace Radiation;
int main(int argc, char *argv[]) {
	installFpeHandler();
	return Catch::Session().run(argc, argv);
}

EquationOfState makeEquationOfState(double Γ, double μ);

struct RadTestCase {
	double rho;
	double kappaA;
	double kappaS;
	double beta[3]; // fluid velocity / c
	double E;       // radiation energy density
	double Ffrac[3];       // flux fractions relative to Er
	double Tg;      // gas temperature for initial ε
	double dt;
	const char *name;
};

TEST_CASE("Implicit energy solve conserves total energy in all regimes", "[implicitEnergySolve]") {
	EquationOfState eos(MolarMassType(1.0));

	// @formatter:off
    std::vector<RadTestCase> cases = {
        // idx | ρ     κₐ     κₛ       β={vx,vy,vz}     Er       Ffrac={Fx,Fy,Fz}              Tg     dt    description
        /* 0 */ {1.0,  1.0,   0.0,     {0.0,0.0,0.0},   1.0,     {0.01,0.0,0.0},                1e6,   1.0, "Thick diffusion"},
        /* 1 */ {1.0,  1.0,   0.0,     {0.5,0.0,0.0},   1.0,     {0.01,0.0,0.0},                1e6,   1.0, "Thick diffusion, boosted"},
        /* 2 */ {1.0,  0.8,   0.2,     {0.0,0.0,0.0},   1.0,     {0.3,0.0,0.0},                 1e6,   1.0, "Semi-transparent"},
        /* 3 */ {1.0,  0.0,   1e-6,    {0.0,0.0,0.0},   1.0,     {0.999,0.0,0.0},               1e6,   1.0, "Free-streaming"},
        /* 4 */ {1.0,  1e-3,  0.0,     {0.0,0.0,0.0},   1.0,     {0.9,0.0,0.0},                 1e6,   1.0, "Thin with weak absorption"},
        /* 5 */ {1.0,  100.0, 0.0,     {0.0,0.0,0.0},  10.0,     {0.1,0.0,0.0},                 1e6,   1.0, "Stiff absorption dominated"},
        /* 6 */ {1.0,  0.01,  100.0,   {0.0,0.0,0.0},   1.0,     {0.6,0.0,0.0},                 1e6,   1.0, "Stiff scattering dominated"},
        /* 7 */ {1.0,  1.0,   0.0,     {0.0,0.0,0.0}, 100.0,     {0.1,0.0,0.0},                 1e5,   1.0, "Radiation-dominated gas"},
        /* 8 */ {1.0,  1.0,   0.0,     {0.0,0.0,0.0},  0.01,     {0.1,0.0,0.0},                 1e6,   1.0, "Matter-dominated"},
        /* 9 */ {1.0,  1.0,   0.0,     {0.95,0.0,0.0},  1.0,     {0.0,0.0,0.0},                 1e6,   1.0, "Ultra-relativistic flow"},
        /*10 */ {1.0,  0.0,   1e-6,    {0.0,0.0,0.0},   1.0,     {+1.0,0.0,0.0},                1e6,   1.0, "Counter-beam +x"},
        /*11 */ {1.0,  0.0,   1e-6,    {0.0,0.0,0.0},   1.0,     {-1.0,0.0,0.0},                1e6,   1.0, "Counter-beam -x"},
        /*12 */ {1.0,  0.5,   0.5,     {0.3,0.0,0.0},   1.0,     {0.7/sqrt(2),0.7/sqrt(2),0.0},1e6,   1.0, "Oblique anisotropy"},
        /*13 */ {1.0,  0.05,  0.05,    {0.7,0.0,0.0},   1.0,     {0.999,0.0,0.0},               1e6,   1.0, "Flux cap edge"},
        /*14 */ {1.0,  1.0,   0.0,     {0.0,0.0,0.0}, 1e-12,     {0.0,0.0,0.0},                 1e7,   1.0, "Emission-only"},
        /*15 */ {1.0,  1.0,   0.0,     {0.0,0.0,0.0},   1.0,     {0.2,0.0,0.0},                 1e3,   1.0, "Absorption-only"},
        /*16 */ {1.0,  0.0,   1.0,     {0.8,0.0,0.0},   1.0,     {0.6,0.0,0.0},                 1e6,   1.0, "Scattering-only"}
    };
    	// @formatter:on

	for (int i = 0; i < cases.size(); i++) {
		auto const &tc = cases[i];
		Opacity opac;
		GasPrimitive gasPrim;
		RadConserved radCon;

		gasPrim.ρ = MassDensityType(tc.rho);
		gasPrim.ε = eos.temperature2energy(gasPrim.ρ, TemperatureType(tc.Tg));
		gasPrim.β[0] = tc.beta[0];
		gasPrim.β[1] = tc.beta[1];
		gasPrim.β[2] = tc.beta[2];

		radCon.Er = temperature2radiationEnergy(TemperatureType(tc.E)); // treat E as equivalent radiation T^4
		radCon.F[0] = c * tc.Ffrac[0] * radCon.Er;
		radCon.F[1] = c * tc.Ffrac[1] * radCon.Er;
		radCon.F[2] = c * tc.Ffrac[2] * radCon.Er;

		opac.κₐ = SpecificAreaType(tc.kappaA);
		opac.κₛ = SpecificAreaType(tc.kappaS);

		auto const radCon0 = radCon;
		auto gasCon = gasPrim.toConserved(eos);

		implicitRadiationSolve(gasCon, radCon, opac, eos, TimeType(tc.dt));

		auto deltaE = radCon0.Er - radCon.Er;

		double Eg0 = (gasPrim.ρ * gasPrim.ε).value();
		double Er0 = radCon0.Er.value();
		double Eg1 = Eg0 - deltaE.value();
		double Er1 = Er0 + deltaE.value();

		INFO("Case: " << tc.name);
		INFO(
				"rho=" << tc.rho << " kappaA=" << tc.kappaA << " kappaS=" << tc.kappaS << " beta=(" << tc.beta[0] << "," << tc.beta[1] << "," << tc.beta[2] << ")" << " E=" << tc.E << " dt=" << tc.dt);
		REQUIRE(Eg0 + Er0 == Approx(Eg1 + Er1).epsilon(1e-10));
	}
}

