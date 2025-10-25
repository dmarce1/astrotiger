#include "rad_conserved.hpp"
#include "catch2/catch_test_macros.hpp"
#include "matrix.hpp"
#include "vector.hpp"
#include <cmath>
#include <iostream>

TEST_CASE("M1 eigenSystem consistent with Jacobian", "[srhd][eigen][jacobian]") {
	using std::abs;
	using Type = double;
	constexpr Type tol = 4 * (std::numeric_limits < Type > ::epsilon());
	constexpr int dimensionCount = 3;
	constexpr int fieldCount = 1 + dimensionCount;
	static constexpr PhysicalConstants<Type> pc { };
	constexpr auto c = pc.c.value();
	RadConserved<Type, dimensionCount> rad;
	rad.setEnergy(1.0);
	rad.setFlux(0, c * 0.49);
	rad.setFlux(1, c * 0.09);
	rad.setFlux(2, -c * 0.04);
	for (int dim = 0; dim < dimensionCount; dim++) {
		SquareMatrix<Type, fieldCount> J = rad.jacobian(dim);
		auto [λ, R] = rad.eigenSystem(dim);
		SquareMatrix<Type, fieldCount> Λ;
		for (int i = 0; i < fieldCount; i++) {
			for (int j = 0; j < fieldCount; j++) {
				Λ(i, j) = ((i == j) ? Type(λ[i]) : Type(0.0));
			}
		}
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
//		std::cout << (J * R - R * Λ);
		auto JR = J * R;
		auto RL = R * Λ;
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
		printf("\n");
	}
}
