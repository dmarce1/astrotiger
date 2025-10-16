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

	REQUIRE(primOut.getMassDensity().value() == Catch::Approx(primIn.getMassDensity().value()).epsilon(tol));
	REQUIRE(primOut.getSpecificEnergy().value() == Catch::Approx(primIn.getSpecificEnergy().value()).epsilon(tol));
	for (int i = 0; i < ndim; i++) {
		REQUIRE(primOut.getVelocity(i).value() == Catch::Approx(primIn.getVelocity(i).value()).epsilon(tol));
	}
}

TEST_CASE("Ultra-relativistic γ~100", "[con2prim]") {
	using namespace Constants;
	GasPrimitive<double, ndim> prim;
	auto const eos = makeEquationOfState(1);

	prim.setMassDensity(MassDensityType(1.0));
	prim.setVelocity({ c * DimensionlessType(0.99995), c * DimensionlessType(0.0), c * DimensionlessType(0.0) });
	prim.setSpecificEnergy(0.00001 * c2);

	roundTripCheck(prim, eos, 1e-6);
}

TEST_CASE("Normal gas", "[con2prim]") {
	using namespace Constants;
	GasPrimitive<double, ndim> prim;
	auto const eos = makeEquationOfState(1);

	prim.setMassDensity(MassDensityType(1.0));
	prim.setVelocity({ c * DimensionlessType(0.01), c * DimensionlessType(0.0), c * DimensionlessType(0.0) });
	prim.setSpecificEnergy(0.1 * c2);

	roundTripCheck(prim, eos, 1e-5);
}

TEST_CASE("Cold gas near rest", "[con2prim]") {
	using namespace Constants;
	GasPrimitive<double, ndim> prim;
	auto const eos = makeEquationOfState(1);

	prim.setMassDensity(MassDensityType(1.0));
	prim.setVelocity({ c * DimensionlessType(0.0), c * DimensionlessType(0.0), c * DimensionlessType(0.0) });
	prim.setSpecificEnergy(1e-6);

	roundTripCheck(prim, eos);
}

TEST_CASE("Mildly relativistic γ~2", "[con2prim]") {
	using namespace Constants;
	GasPrimitive<double, ndim> prim;
	auto const eos = makeEquationOfState(1);

	prim.setMassDensity(MassDensityType(1.0));
	prim.setVelocity({ c * DimensionlessType(0.866), c * DimensionlessType(0.0), c * DimensionlessType(0.0) }); // γ ≈ 2
	prim.setSpecificEnergy(0.1 * c2);

	roundTripCheck(prim, eos, 1e-8);
}

TEST_CASE("Ultra-relativistic γ~1000", "[con2prim]") {
	using namespace Constants;
	GasPrimitive<double, ndim> prim;
	auto const eos = makeEquationOfState(1);

	prim.setMassDensity(MassDensityType(1.0));
	prim.setVelocity({ c * DimensionlessType(0.9999995), c * DimensionlessType(0.0), c * DimensionlessType(0.0) });
	prim.setSpecificEnergy(0.5 * c2);

	roundTripCheck(prim, eos, 1e-6);
}

TEST_CASE("Ultra-relativistic γ~10", "[con2prim]") {
	using namespace Constants;
	GasPrimitive<double, ndim> prim;
	auto const eos = makeEquationOfState(1);

	prim.setMassDensity(MassDensityType(1.0));
	prim.setVelocity({ c * DimensionlessType(0.995), c * DimensionlessType(0.0), c * DimensionlessType(0.0) }); // γ ≈ 10
	prim.setSpecificEnergy(0.5 * c2);

	roundTripCheck(prim, eos, 1e-6);
}

TEST_CASE("Hot gas (ε >> c²)", "[con2prim]") {
	using namespace Constants;
	GasPrimitive<double, ndim> prim;
	auto const eos = makeEquationOfState(1);

	prim.setMassDensity(MassDensityType(1.0));
	prim.setVelocity({ c * DimensionlessType(0.2), c * DimensionlessType(0.1), c * DimensionlessType(0.0) });
	prim.setSpecificEnergy(50.0 * c2);

	roundTripCheck(prim, eos, 1e-5);
}

TEST_CASE("Low-density / vacuum-like", "[con2prim]") {
	using namespace Constants;
	GasPrimitive<double, ndim> prim;
	auto const eos = makeEquationOfState(1);

	prim.setMassDensity(MassDensityType(1e-12));
	prim.setVelocity({ c * DimensionlessType(0.0), c * DimensionlessType(0.0), c * DimensionlessType(0.0) });
	prim.setSpecificEnergy(1.0);

	roundTripCheck(prim, eos, 1e-6);
}
