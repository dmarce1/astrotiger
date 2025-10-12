#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "numbers/autodiff.hpp" // your header with the class template

using Catch::Approx;

// ----------------- Helpers -----------------
template<int O, int D>
static void forAllMultiIndices(std::function<void(std::array<int, D> const&)> const &f) {
	std::array<int, D> α { };
	α.fill(0);
	while (true) {
		f(α);
		// increment base-O counter over D dims: 0..O-1 in each slot
		int d = D - 1;
		while (d >= 0) {
			α[d] += 1;
			if (α[d] < O) {
				break;
			}
			α[d] = 0;
			--d;
		}
		if (d < 0) {
			break;
		}
	}
}

template<int O, int D>
static double getCoeff(FwdAutoDiffDouble<O, D> const &a, int n) {
	static_assert(D == 1, "D==1 overload only");
	return a[n];
}

template<int O, int D>
static double getCoeff(FwdAutoDiffDouble<O, D> const &a, std::array<int, D> const &α) {
	if constexpr (D == 1) {
		return a[α[0]];
	} else {
		return a[α];
	}
}

template<int O, int D>
static void expectSeriesEqual(FwdAutoDiffDouble<O, D> const &A, FwdAutoDiffDouble<O, D> const &B, double tol = 1e-12) {
	if constexpr (D == 1) {
		for (int n = 0; n < O; ++n) {
			REQUIRE(getCoeff<O, D>(A, n) == Approx(getCoeff<O, D>(B, n)).margin(tol));
		}
	} else {
		forAllMultiIndices<O, D>([&](std::array<int, D> const &α) {
			REQUIRE(getCoeff<O, D>(A, α) == Approx(getCoeff<O, D>(B, α)).margin(tol));
		});
	}
}

//// ----------------- Main -----------------
//int main(int argc, char* argv[]) {
//	return Catch::Session().run(argc, argv);
//}

TEST_CASE("FwdAutoDiffDouble<O=6,D=1>: basic identities", "[FwdAutoDiffDouble][univariate]") {
	constexpr int O = 6;
	constexpr int D = 1;
	using AD = FwdAutoDiffDouble<O, D>;

	// choose a positive expansion point so log/sqrt/cbrt are well-defined
	double const a = 1.3;

	AD x = AD::independent(a, 0);
	AD y;
//	SECTION("exp(log(x)) == x (series equality)") {
//		y = exp(x);
//	}
	SECTION("exp(log(x)) == x (series equality)") {
		FwdAutoDiffDouble<O, D> z = log(x);
		FwdAutoDiffDouble<O, D> y = exp(z);
		expectSeriesEqual<O, D>(y, x, 1e-12);
	}
	SECTION("x * (1/x) == 1 (series equality)") {
		auto y = x * (1.0 / x);
		AD one(1.0);
		expectSeriesEqual<O, D>(y, one, 1e-12);
	}

	SECTION("sqrt(x)^2 == x (within truncation)") {
		auto s = sqrt(x);
		auto y = s * s;
		expectSeriesEqual<O, D>(y, x, 1e-12);
	}

	SECTION("cbrt(x)^3 == x (within truncation)") {
		auto r = cbrt(x);
		auto y = r * r * r;
		expectSeriesEqual<O, D>(y, x, 1e-12);
	}

	SECTION("scalar ops preserve series structure") {
		auto y = (x * 2.5) / 5.0; // net scale 0.5
		AD half = x / 2.0;
		expectSeriesEqual<O, D>(y, half, 1e-12);
	}
}

TEST_CASE("FwdAutoDiffDouble<O=5,D=2>: multivariate composition/product identities", "[FwdAutoDiffDouble][bivariate]") {
	constexpr int O = 5;
	constexpr int D = 2;
	using AD = FwdAutoDiffDouble<O, D>;

	// positive expansion point for both dims
	double const ax = 0.9;
	double const ay = 1.2;

	AD x = AD::independent(ax, 0);
	AD y = AD::independent(ay, 1);

	SECTION("exp(log(x*y)) == x*y (series equality)") {
		AD xy = x * y;
		AD f  = exp(log(xy));
		std::cout << xy << "\n" << f << "\n";
		expectSeriesEqual<O, D>(f, xy, 1e-12);
	}

	SECTION("(sqrt(x*y))^2 == x*y (within truncation)") {
		auto xy = x * y;
		auto s  = sqrt(xy);
		auto back = s * s;
		expectSeriesEqual<O, D>(back, xy, 1e-12);
	}

	SECTION("(x / y) * y == x (series equality)") {
		auto q = x / y;
		auto back = q * y;
		expectSeriesEqual<O, D>(back, x, 1e-12);
	}
}

// Optional: quick numerical sanity via evaluated series at a point (uses operator())
TEST_CASE("Evaluation at a point agrees with constructed identities", "[FwdAutoDiffDouble][eval]") {
	constexpr int O = 6;
	constexpr int D = 1;
	using AD = FwdAutoDiffDouble<O, D>;

	double const a = 1.1;
	AD x = AD::independent(a, 0);

	auto id1 = exp(log(x));
	auto id2 = (sqrt(x) * sqrt(x));
	auto id3 = x * (1.0 / x);

	std::array<double, D> evalPoint{ a };

	REQUIRE(id1(evalPoint) == Approx(x(evalPoint)).epsilon(1e-12));
	REQUIRE(id2(evalPoint) == Approx(x(evalPoint)).epsilon(1e-12));
	REQUIRE(id3(evalPoint) == Approx(1.0).epsilon(1e-12));
}
