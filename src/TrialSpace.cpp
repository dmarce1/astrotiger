#include "Definitions.hpp"
#include "Quadrature.hpp"
#include "TrialSpace.hpp"

#include <algorithm>
#include <vector>

using namespace Constants;

std::vector<std::vector<Real>> getLagrangePolynomials(int N) {
	static Memoize<std::vector<std::vector<Real>>> memoiy([](int N) {
		std::vector<std::vector<Real>> polynomials;
		auto const quadraturePoints = getLobattoQuadraturePoints(N);
		for (int j = 0; j < N; j++) {
			std::vector<Real> P(N, Real(0));
			auto const myPoint = quadraturePoints[j];
			P[0] = one;
			for (int i = 0; i < N; i++) {
				if (i == j) {
					continue;
				}
				auto const quadPoint = quadraturePoints[i];
				Real const root = quadPoint.x;
				Real const den = myPoint.x - root;
				if (den == zero) {
					continue;
				}
				auto const Pold = P;
				std::shift_right(P.begin(), P.end(), 1);
				P[0] = zero;
				for (int i = 0; i < (int) P.size(); i++) {
					P[i] -= root * Pold[i];
					P[i] /= den;
				}
			}
			polynomials.push_back(P);
		}
		return polynomials;
	});
	return memoiy(N);
}

