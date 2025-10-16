/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/
#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "radiation.hpp"

constexpr int ndim = 3;

EquationOfState<double> makeEquationOfState(double μ) {
	return EquationOfState<double>(MolarMassType<double>(μ));
}

static void roundTripCheck(GasPrimitive<double, ndim> const &primIn, EquationOfState<double> const &eosIn, double tol = 1e-10) {
	enableFPE();
	auto con = primIn.toConserved(eosIn);
	auto primOut = con.toPrimitive(eosIn);
	printf("%e %e\n", (double) primIn.ρ.value(), (double) primOut.ρ.value());
	printf("%e %e\n", (double) primIn.ε.value(), (double) primOut.ε.value());
	for (int i = 0; i < 3; i++) {
		printf("%e %e\n", (double) primIn.β[i].value(), (double) primOut.β[i].value());
	}
	REQUIRE(primOut.ρ.value() == Catch::Approx(primIn.ρ.value()).epsilon(tol));
	REQUIRE(primOut.ε.value() == Catch::Approx(primIn.ε.value()).epsilon(tol));
	for (int i = 0; i < 3; i++) {
		REQUIRE(primOut.β[i].value() == Catch::Approx(primIn.β[i].value()).epsilon(tol));
	}
}

TEST_CASE("Ultra-relativistic γ~100", "[con2prim]") {
	using namespace Constants;
	GasPrimitive<double, ndim> prim;
	EquationOfState<double> const eos = makeEquationOfState(1);
	prim.ρ = MassDensityType(1.0);
	prim.β = { (DimensionlessType<double>) 0.99995, (DimensionlessType<double>) 0.0, (DimensionlessType<double>) 0.0 };
	prim.ε = 0.00001 * c2;

	roundTripCheck(prim, eos, 1e-6);
}

TEST_CASE("Normal gas ", "[con2prim]") {
	using namespace Constants;
	GasPrimitive<double, ndim> prim;
	EquationOfState<double> const eos = makeEquationOfState(1);
	prim.ρ = MassDensityType<double>(1.0);
	prim.β = { (DimensionlessType<double>) 0.01, (DimensionlessType<double>) 0.0, (DimensionlessType<double>) 0.0 };
	prim.ε = 1.0e-1 * c2;

	roundTripCheck(prim, eos, 1e-5);
}

TEST_CASE("Cold gas near rest", "[con2prim]") {
	GasPrimitive<double, ndim> prim;
	EquationOfState<double> const eos = makeEquationOfState(1);
	prim.ρ = MassDensityType(1.0);
	prim.β = { (DimensionlessType<double>) 0.0, (DimensionlessType<double>) 0.0, (DimensionlessType<double>) 0.0 };
	prim.ε = SpecificEnergyType<double>(1e-6);

	roundTripCheck(prim, eos);
}

TEST_CASE("Mildly relativistic γ~2", "[con2prim]") {
	using namespace Constants;
	GasPrimitive<double, ndim> prim;
	EquationOfState<double> const eos = makeEquationOfState(1);
	prim.ρ = MassDensityType<double>(1.0);
	prim.β = { (DimensionlessType<double>) 0.866, (DimensionlessType<double>) 0.0, (DimensionlessType<double>) 0.0 }; // gamma ~ 2
	prim.ε = 0.1 * c2;

	roundTripCheck(prim, eos, 1e-8);
}

TEST_CASE("Ultra-relativistic γ~1000", "[con2prim]") {
	using namespace Constants;
	GasPrimitive<double, ndim> prim;
	EquationOfState<double> const eos = makeEquationOfState(1);
	prim.ρ = MassDensityType(1.0);
	prim.β = { (DimensionlessType<double>) 0.9999995, (DimensionlessType<double>) 0.0, (DimensionlessType<double>) 0.0 };
	prim.ε = 0.5 * c2;

	roundTripCheck(prim, eos, 1e-6);
}

TEST_CASE("Ultra-relativistic γ~10", "[con2prim]") {
	using namespace Constants;
	GasPrimitive<double, ndim> prim;
	EquationOfState<double> const eos = makeEquationOfState(1);
	prim.ρ = MassDensityType(1.0);
	prim.β = { DimensionlessType(0.995), DimensionlessType(0.0), DimensionlessType(0.0) }; // gamma ~ 10
	prim.ε = 0.5 * c2;

	roundTripCheck(prim, eos, 1e-6); // allow looser tolerance
}

TEST_CASE("Hot gas (ε >> c^2)", "[con2prim]") {
	using namespace Constants;
	GasPrimitive<double, ndim> prim;
	EquationOfState<double> const eos = makeEquationOfState(1);
	prim.ρ = MassDensityType(1.0);
	prim.β = { DimensionlessType(0.2), DimensionlessType(0.1), DimensionlessType(0.0) };
	prim.ε = 50.0 * c2;

	roundTripCheck(prim, eos, 1e-5);
}

TEST_CASE("Low-density / vacuum-like", "[con2prim]") {
	GasPrimitive<double, ndim> prim;
	EquationOfState<double> const eos = makeEquationOfState(1);
	prim.ρ = MassDensityType(1e-12);
	prim.β = { DimensionlessType(0.0), DimensionlessType(0.0), DimensionlessType(0.0) };
	prim.ε = SpecificEnergyType(1.0);

	roundTripCheck(prim, eos, 1e-6);
}
