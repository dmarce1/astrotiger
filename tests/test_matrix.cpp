/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

/******************************************************************************
 *  Catch2 Unit Tests for Matrix.hpp
 *  Thorough coverage of constructors, operators, and mathematical utilities.
 ******************************************************************************/
#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "rad_conserved.hpp"
#include "fpe.hpp"

#include "matrix.hpp"

using Catch::Approx;
using Catch::Matchers::WithinRel;

TEST_CASE("Matrix basic construction and indexing", "[matrix]") {
	using M2 = Matrix<double, 2, 2>;
	M2 A { 1.0, 2.0, 3.0, 4.0 };

	REQUIRE(A(0, 0) == 1.0);
	REQUIRE(A(0, 1) == 2.0);
	REQUIRE(A(1, 0) == 3.0);
	REQUIRE(A(1, 1) == 4.0);

	A(0, 1) = 5.0;
	REQUIRE(A(0, 1) == 5.0);
}

TEST_CASE("Matrix symmetric indexing and assignment", "[matrix]") {
	using MSym = SquareMatrix<double, 3, SymmetryType::symmetric>;
	MSym S(0.0);
	S(0, 0) = 1.0;
	S(1, 0) = 2.0;
	S(1, 1) = 3.0;
	S(2, 0) = 4.0;
	S(2, 1) = 5.0;
	S(2, 2) = 6.0;

	REQUIRE(S(0, 1) == 2.0);
	REQUIRE(S(1, 0) == 2.0);
	REQUIRE(S(2, 0) == 4.0);
	REQUIRE(S(0, 2) == 4.0);
}

TEST_CASE("Matrix antisymmetric indexing", "[matrix]") {
	using MAnti = SquareMatrix<double, 3, SymmetryType::antisymmetric>;
	MAnti A(0.0);
	A(1, 0) = 5.0;

	REQUIRE(A(1, 0) == 5.0);
	REQUIRE(A(0, 1) == -5.0);
	REQUIRE(A(0, 0) == 0.0);
}

TEST_CASE("Matrix addition and subtraction", "[matrix]") {
	using M2 = Matrix<double, 2, 2>;
	M2 A { 1, 2, 3, 4 };
	M2 B { 5, 6, 7, 8 };

	auto C = A + B;
	REQUIRE(C(0, 0) == 6.0);
	REQUIRE(C(1, 1) == 12.0);

	auto D = B - A;
	REQUIRE(D(0, 0) == 4.0);
	REQUIRE(D(1, 1) == 4.0);

	A += B;
	REQUIRE(A(0, 0) == 6.0);
	A -= B;
	REQUIRE(A(0, 0) == 1.0);
}

TEST_CASE("Matrix scalar arithmetic", "[matrix]") {
	using M2 = Matrix<double, 2, 2>;
	M2 A { 1, 2, 3, 4 };

	auto B = A * 2.0;
	REQUIRE(B(1, 1) == 8.0);

	auto C = B / 2.0;
	REQUIRE(C(1, 1) == 4.0);

	A *= 3.0;
	REQUIRE(A(0, 0) == 3.0);

	A /= 3.0;
	REQUIRE(A(0, 0) == 1.0);
}

TEST_CASE("Matrix transpose", "[matrix]") {
	using M2 = Matrix<double, 2, 3>;
	M2 A { 1, 2, 3, 4, 5, 6 };
//	M2 A{1, 4,
//		 2, 5,
//		 3, 6};
	auto T = transpose(A);

	REQUIRE(T(0, 0) == 1.0);
	REQUIRE(T(1, 0) == 2.0);
	REQUIRE(T(0, 1) == 4.0);
	REQUIRE(T(2, 1) == 6.0);
}

TEST_CASE("Matrix trace", "[matrix]") {
	using M2 = SquareMatrix<double, 2>;
	M2 A { 1, 2, 3, 4 };
	double tr = trace(A);
	REQUIRE(tr == Approx(5.0));
}

TEST_CASE("Matrix determinant 2x2", "[matrix]") {
	using M2 = SquareMatrix<double, 2>;
	M2 A { 1, 2, 3, 4 };
	double det = determinant(A);
	REQUIRE(det == Approx(-2.0));
}

TEST_CASE("Matrix determinant 3x3", "[matrix]") {
	using M3 = SquareMatrix<double, 3>;
	M3 A { 6, 1, 1, 4, -2, 5, 2, 8, 7 };
	double det = determinant(A);
	REQUIRE(det == Approx(-306.0));
}

TEST_CASE("Matrix comatrix and adjoint", "[matrix]") {
	using M2 = SquareMatrix<double, 2>;
	M2 A { 1, 2, 3, 4 };
//	M2 A{4, -2,
//		 -3, 1};
	auto coA = comatrix(A);
	auto adjA = adjoint(A);

	REQUIRE(coA(0, 0) == 4.0);
	REQUIRE(coA(1, 1) == 1.0);
	REQUIRE(adjA(0, 1) == -2.0);
	REQUIRE(adjA(1, 0) == -3.0);
}

TEST_CASE("Matrix inverse", "[matrix]") {
	using M2 = SquareMatrix<double, 2>;
	M2 A { 1, 2, 3, 4 };
	auto invA = inverse(A);
	auto I = A * invA;

	REQUIRE(I(0, 0) == Approx(1.0).margin(1e-12));
	REQUIRE(I(1, 1) == Approx(1.0).margin(1e-12));
	REQUIRE(I(0, 1) == Approx(0.0).margin(1e-12));
	REQUIRE(I(1, 0) == Approx(0.0).margin(1e-12));
}

TEST_CASE("Matrix inversion failure throws", "[matrix][exception]") {
	using M2 = SquareMatrix<double, 2>;
	M2 A { 1, 2, 2, 4 }; // singular

	REQUIRE_THROWS_AS(inverse(A), MatrixException);
}

TEST_CASE("Identity matrix creation", "[matrix]") {
	auto I = identity<double, 3>();
	REQUIRE(I(0, 0) == 1.0);
	REQUIRE(I(1, 1) == 1.0);
	REQUIRE(I(2, 2) == 1.0);
	REQUIRE(I(0, 1) == 0.0);
	REQUIRE(I(1, 2) == 0.0);
}

TEST_CASE("Vector outer product", "[matrix]") {
	Vector<double, 2> a { 1.0, 2.0 };
	Vector<double, 3> b { 3.0, 4.0, 5.0 };
	auto M = a * b;

	REQUIRE(M(0, 0) == 3.0);
	REQUIRE(M(1, 2) == 10.0);
}

TEST_CASE("sqr(Vector) produces symmetric matrix", "[matrix]") {
	Vector<double, 3> a { 1.0, 2.0, 3.0 };
	auto M = sqr(a);

	REQUIRE(M(0, 0) == 1.0);
	REQUIRE(M(2, 1) == Approx(6.0));
	REQUIRE(M(1, 2) == Approx(6.0));
}

TEST_CASE("Symmetric and antisymmetric decomposition", "[matrix]") {
	using M3 = SquareMatrix<double, 3>;
	M3 A { 0, 1, 2, 3, 4, 5, 6, 7, 8 };

	auto symA = symmetrize(A);
	auto antiA = antisymmetrize(A);

	REQUIRE(symA(1, 2) == Approx((A(1, 2) + A(2, 1)) / 2));
	REQUIRE(antiA(1, 2) == Approx((A(1, 2) - A(2, 1)) / 2));
}

TEST_CASE("toString produces formatted output", "[matrix]") {
	using M2 = Matrix<double, 2, 2>;
	M2 A { 1.0, 2.0, 3.0, 4.0 };
	std::ostringstream os;
	os << A;
	std::string s = os.str();

	REQUIRE(s.find("1.000e+00") != std::string::npos);
	REQUIRE(s.find("4.000e+00") != std::string::npos);
	REQUIRE(s.find('|') != std::string::npos);
}

TEST_CASE("SignedReference assignment and arithmetic", "[signedref]") {
	double val = 10.0;
	SignedReference<double> ref(+1, val);
	ref = 5.0;
	REQUIRE(val == 5.0);
	ref += 2.0;
	REQUIRE(val == 7.0);
	ref -= 3.0;
	REQUIRE(val == 4.0);

	SignedReference<double> refNeg(-1, val);
	refNeg = 6.0;
	REQUIRE(val == -6.0);
}

TEST_CASE("multSymmetry combinations", "[symmetry]") {
using enum SymmetryType;
		REQUIRE(multSymmetry(symmetric, symmetric) == symmetric);
	REQUIRE(multSymmetry(symmetric, antisymmetric) == antisymmetric);
	REQUIRE(multSymmetry(antisymmetric, antisymmetric) == symmetric);
	REQUIRE(multSymmetry(asymmetric, symmetric) == asymmetric);
}
