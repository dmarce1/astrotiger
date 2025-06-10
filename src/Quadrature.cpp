#include "Definitions.hpp"
#include "Quadrature.hpp"

#include <algorithm>

std::vector<QuadraturePoint> getLobattoQuadraturePoints(int N) {
	static Memoize<std::vector<QuadraturePoint>> memoiy([](int N) {
		using namespace Constants;
		std::vector<QuadraturePoint> points;
		Real const NNm1 = Real(N * (N - 1));
		Real x, w;
		points.push_back( { Real(-1), two / NNm1 });
		for (int j = 1; j < N - 1; j++) {
			Real const Np1 = Real(N + 1);
			Real const Nm1 = Real(N - 1);
			x = std::cos(pi * (one - Real(j) / Nm1));
			Real prev;
			do {
				prev = x;
				Real const Pnm1 = std::legendre(N - 1, x);
				Real const Pn = std::legendre(N, x);
				Real const x2 = squared(x);
				Real const num = (one - x2) * (Pn - x * Pnm1);
				Real const den = (Np1 * x2 - Nm1) * Pnm1 - two * x * Pn;
				x += num / den;
			} while (Real(x) != Real(prev));
			w = two / (NNm1 * squared(std::legendre(N - 1, x)));
			points.push_back( { Real(x), Real(w) });
		}
		points.push_back( { Real(+1), two / NNm1 });
		return points;
	});
	return memoiy(N);
}

