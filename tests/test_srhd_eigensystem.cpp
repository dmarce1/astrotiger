#include "gas_primitive.hpp"
#include "catch2/catch_test_macros.hpp"
#include "eos.hpp"
#include "matrix.hpp"
#include "vector.hpp"
#include <cmath>
#include <iostream>

TEST_CASE("SRHD primitive eigenSystem consistent with Jacobian", "[srhd][eigen][jacobian]") {
	using std::abs;
	using Type = double;
	constexpr Type tol = 16 * std::numeric_limits < Type > ::epsilon();
	constexpr int dimensionCount = 3;
	constexpr int fieldCount = 2 + dimensionCount;
	static constexpr PhysicalConstants<Type> pc { };
	constexpr auto c = pc.c;

	EquationOfState<Type> eos(Type(5.0 / 3.0));
	GasPrimitive<Type, dimensionCount> prim;
	prim.setMassDensity(MassDensityType<Type>(1.0));
	prim.setSpecificEnergy(SpecificEnergyType<Type>(1.0e-1 * c * c));
	prim.setVelocity( { 7e-1 * c, 1e-1 * c, -1e-2 * c });
	for (int dim = 0; dim < dimensionCount; dim++) {
		SquareMatrix<Type, fieldCount> J = prim.jacobian(eos, dim);
		auto [λ, R] = prim.eigenSystem(eos, dim);
		SquareMatrix<Type, fieldCount> Λ;
		for (int i = 0; i < fieldCount; i++) {
			for (int j = 0; j < fieldCount; j++) {
				Λ(i, j) = ((i == j) ? Type(λ[i]) : Type(0.0));
			}
		}
		auto JR = J * R;
		auto RL = R * Λ;
//		std::cout << "Λ = \n";
//		std::cout << (Λ);
//		std::cout << "R = \n";
//		std::cout << (R);
//		std::cout << "J = \n";
//		std::cout << (J);
//		std::cout << "J.R = \n";
//		std::cout << (J * R);
//		std::cout << "R * Λ = \n";
//		std::cout << (R * Λ);
//		std::cout << "J.R - R * Λ = \n";
//		std::cout << (JR - RL) / frobeniusNorm(RL);
		for (int i = 0; i < fieldCount; i++) {
			for (int j = 0; j < fieldCount; j++) {
				CHECK((abs((JR(i, j) - RL(i, j))) / c) < tol);
			}
		}
		for (int i = 0; i < fieldCount; i++) {
			auto Ri = R.getColumn(i);
			auto JRi = J * Ri;
			auto λRi = Ri * λ[i];
			for (int k = 0; k < fieldCount; k++) {
				CHECK((abs((JRi[k] - λRi[k])) / c) < tol);
			}
		}
		Type detR = determinant(R);
		CHECK(abs(detR) > 1e-10);
	}
}
