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
	constexpr Type tol = 8 * sqrt(std::numeric_limits<Type>::epsilon());
	constexpr int dimensionCount = 3;
	constexpr int fieldCount = 2 + dimensionCount;
	static constexpr PhysicalConstants<Type> pc { };
	constexpr auto c = pc.c;

	// Construct an ideal-gas EOS
	EquationOfState<Type> eos(Type(5.0 / 3.0));

	// --- Step 1: Build primitive state using setters ---
	GasPrimitive<Type, dimensionCount> prim;
	prim.setMassDensity(MassDensityType<Type>(1.0));
	prim.setSpecificEnergy(SpecificEnergyType<Type>(1.0e-3 * c * c));
	prim.setVelocity(Vector<VelocityType<Type>, dimensionCount> { 0.05 * c, -0.03 * c, 0.1 * c });

	// --- Step 2: Loop over each spatial direction ---
	for (int dim = 0; dim < dimensionCount; dim++) {
		SquareMatrix<Type, fieldCount> J = prim.jacobian(eos, dim);
		auto [λ, R] = prim.eigenSystem(eos, dim);

		// --- Step 3: Build Λ (diagonal matrix of eigenvalues) ---
		SquareMatrix<Type, fieldCount> Λ;
		for (int i = 0; i < fieldCount; i++) {
			for (int j = 0; j < fieldCount; j++) {
				Λ(i, j) = ((i == j) ? Type(λ[i]) : Type(0.0));
			}
		}
//		std::cout << "R = \n";
//		std::cout << toString(R);
//		std::cout << "J = \n";
//		std::cout << toString(J);
//		std::cout << "J.R = \n";
//		std::cout << toString(J * R);
//		std::cout << "R * Λ = \n";
//		std::cout << toString(R * Λ);
//		std::cout << "J.R - R * Λ = \n";
//		std::cout << toString(J * R - R * Λ);

//		std::cout << "R = \n";
//		std::cout << toMathematica(R);
//		std::cout << "J = \n";
//		std::cout << toMathematica(TypeJ);

// --- Step 4: Check J * R ≈ R * Λ ---
		auto JR = J * R;
		auto RL = R * Λ;
		for (int i = 0; i < fieldCount; i++) {
			for (int j = 0; j < fieldCount; j++) {
				CHECK((abs((JR(i, j) - RL(i, j))) / c).value() < tol);
			}
		}

		// --- Step 5: Check each eigenpair individually ---
		for (int i = 0; i < fieldCount; i++) {
			auto Ri = R.getColumn(i);
			auto JRi = J * Ri;
			auto λRi = Ri * λ[i];
			for (int k = 0; k < fieldCount; k++) {
				CHECK((abs((JRi[k] - λRi[k])) / c).value() < tol);
			}
		}

		// --- Step 6: Check eigenvector independence ---
		Type detR = determinant(R);
		CHECK(abs(detR) > 1e-10);
		printf("\n");
	}
}
