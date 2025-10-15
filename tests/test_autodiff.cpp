
#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <iostream>

#include "autodiff.hpp"

using Catch::Approx;

template<int O, int D>
static void forAllMultiIndices(std::function<void(Indices<O, D, int> const&)> const &f) {
	auto α = Indices<O, D, int>::zero();
	while (true) {
		f(α);
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
static double getCoeff(FwdAutoDiff< double, O, std::array<double, D> > const &a, int n) {
	static_assert(D == 1, "D==1 overload only");
	return a[n];
}

template<int O, int D>
static double getCoeff(FwdAutoDiff< double, O, std::array<double, D>> const &a, Indices<O, D, int> const &α) {
	if constexpr (D == 1) {
		return a[α[0]];
	} else {
		return a[α];
	}
}

template<int O, int D>
static void expectSeriesEqual(FwdAutoDiff< double, O, std::array<double, D> > const &A, FwdAutoDiff< double, O, std::array<double, D> > const &B, double tol = 1e-12) {
	if constexpr (D == 1) {
		for (int n = 0; n < O; ++n) {
			REQUIRE(getCoeff<O, D>(A, n) == Approx(getCoeff<O, D>(B, n)).margin(tol));
		}
	} else {
		forAllMultiIndices<O, D>([&](Indices<O, D, int> const &α) {
			REQUIRE(getCoeff<O, D>(A, α) == Approx(getCoeff<O, D>(B, α)).margin(tol));
		});
	}
}


TEST_CASE("FwdAutoDiffDouble<O=6,D=1>: basic identities", "[FwdAutoDiffDouble][univariate]") {
	constexpr int O = 6;
	constexpr int D = 1;
	using AD = FwdAutoDiff< double, O, std::array<double, D> >;

	double const a = 1.3;

	AD x = AD::independent(a, 0);
	AD y;
	SECTION("exp(log(x)) == x (series equality)") {
		FwdAutoDiff< double, O, std::array<double, D> > z = log(x);
		FwdAutoDiff< double, O, std::array<double, D> > y = exp(z);
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
	using AD = FwdAutoDiff< double, O, std::array<double, D> >;

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

TEST_CASE("Evaluation at a point agrees with constructed identities", "[FwdAutoDiffDouble][eval]") {
	constexpr int O = 6;
	constexpr int D = 1;
	using AD = FwdAutoDiff< double, O, std::array<double, D> >;

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
