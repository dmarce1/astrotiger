/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/
#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>   // <-- this is the key include
#include "Radiation.hpp"

using namespace Radiation;



EquationOfState makeEquationOfState(double Γ, double μ) {
	EquationOfState eos;
	eos.Γ = Γ;
	eos.μ = μ * cgs::g/cgs::mol;
	return eos;
}

//static void roundTripCheck(GasPrimitive const &primIn, EquationOfState const &eosIn, double tol = 1e-10) {
//	enableFpeTrapsThisThread();
//	// forward to conserved
//	auto con = primIn.toConserved(eosIn);
//	// back to primitive
//	auto primOut = con.toPrimitive(eosIn);
//	auto units1 = cgs::cm / cgs::s;
//	auto units2 = cgs::g / cgs::cm3;
//	auto units3 = cgs::cm2 / cgs::s2;
//
//	printf( "%e %e\n", (double) primIn.ρ.value(), (double) primOut.ρ.value());
//	printf( "%e %e\n", (double) primIn.ε.value(), (double) primOut.ε.value());
//	for (int i = 0; i < 3; i++) {
//		printf( "%e %e\n", (double) primIn.v[i].value(), (double) primOut.v[i].value());
//	}
//	REQUIRE(primOut.ρ / units2 == Catch::Approx(primIn.ρ / units2).epsilon(tol));
//	REQUIRE(primOut.ε / units3 == Catch::Approx(primIn.ε / units3).epsilon(tol));
//	for (int i = 0; i < 3; i++) {
//		REQUIRE(primOut.v[i] / units1 == Catch::Approx(primIn.v[i] / units1).epsilon(tol));
//	}
//}
//
//TEST_CASE("Ultra-relativistic γ~100", "[con2prim]") {
//	auto units = cgs::cm / cgs::s;
//	GasPrimitive prim;
//	EquationOfState const eos = makeEquationOfState(5. / 3., 1.);
//	prim.ρ = 1.0 * cgs::g / cgs::cm3;
//	prim.v = { 0.99995 * c, 0.0 * units, 0.0 * units };
//	prim.ε = 0.00001 * c2;
//	prim.updateEntropy(eos);
//
//	roundTripCheck(prim, eos, 1e-6);
//}
//
//
//TEST_CASE("Normal gas ", "[con2prim]") {
//	auto units = cgs::cm / cgs::s;
//	GasPrimitive prim;
//	EquationOfState const eos = makeEquationOfState(5. / 3., 1.);
//	prim.ρ = 1.0 * cgs::g / cgs::cm3;
//	prim.v = { 0.01 * c, 0.0 * c, 0.0 * units };
//	prim.ε = 1.0e-1 * c2;
//	prim.updateEntropy(eos);
//
//	roundTripCheck(prim, eos, 1e-5);
//}
//
//TEST_CASE("Cold gas near rest", "[con2prim]") {
//	auto units = cgs::cm / cgs::s;
//	GasPrimitive prim;
//	EquationOfState const eos = makeEquationOfState(5. / 3., 1.);
//	prim.ρ = 1.0 * cgs::g / cgs::cm3;
//	prim.v = { 0.0 * units, 0.0 * units, 0.0 * units };
//	prim.ε = 1e-6 * cgs::erg / cgs::g;
//	prim.updateEntropy(eos);
//
//	roundTripCheck(prim, eos);
//}
//
//TEST_CASE("Mildly relativistic γ~2", "[con2prim]") {
//	auto units = cgs::cm / cgs::s;
//	GasPrimitive prim;
//	EquationOfState const eos = makeEquationOfState(5. / 3., 1.);
//	prim.ρ = 1.0 * cgs::g / cgs::cm3;
//	prim.v = { 0.866 * c, 0.0 * units, 0.0 * units }; // gamma ~ 2
//	prim.ε = 0.1 * c2;
//	prim.updateEntropy(eos);
//
//	roundTripCheck(prim, eos, 1e-8);
//}
//
//TEST_CASE("Ultra-relativistic γ~1000", "[con2prim]") {
//	auto units = cgs::cm / cgs::s;
//	GasPrimitive prim;
//	EquationOfState const eos = makeEquationOfState(5. / 3., 1.);
//	prim.ρ = 1.0 * cgs::g / cgs::cm3;
//	prim.v = { 0.9999995 * c, 0.0 * units, 0.0 * units };
//	prim.ε = 0.5 * c2;
//	prim.updateEntropy(eos);
//
//	roundTripCheck(prim, eos, 1e-6);
//}
//
//
//TEST_CASE("Ultra-relativistic γ~10", "[con2prim]") {
//	auto units = cgs::cm / cgs::s;
//	GasPrimitive prim;
//	EquationOfState const eos = makeEquationOfState(5. / 3., 1.);
//	prim.ρ = 1.0 * cgs::g / cgs::cm3;
//	prim.v = { 0.995 * c, 0.0 * units, 0.0 * units }; // gamma ~ 10
//	prim.ε = 0.5 * c2;
//	prim.updateEntropy(eos);
//
//	roundTripCheck(prim, eos, 1e-6); // allow looser tolerance
//}
//
//TEST_CASE("Hot gas (ε >> c^2)", "[con2prim]") {
//	auto units = cgs::cm / cgs::s;
//	GasPrimitive prim;
//	EquationOfState const eos = makeEquationOfState(5. / 3., 1.);
//	prim.ρ = 1.0 * cgs::g / cgs::cm3;
//	prim.v = { 0.2 * c, 0.1 * c, 0.0 * units };
//	prim.ε = 50.0 * c2;
//	prim.updateEntropy(eos);
//
//	roundTripCheck(prim, eos, 1e-5);
//}
//
//TEST_CASE("Low-density / vacuum-like", "[con2prim]") {
//	auto units = cgs::cm / cgs::s;
//	GasPrimitive prim;
//	EquationOfState const eos = makeEquationOfState(5. / 3., 1.);
//	prim.ρ = 1e-12 * cgs::g / cgs::cm3;
//	prim.v = { 0.0 * units, 0.0 * units, 0.0 * units };
//	prim.ε = 1.0 * cgs::erg / cgs::g;
//	prim.updateEntropy(eos);
//
//	roundTripCheck(prim, eos, 1e-6);
//}
