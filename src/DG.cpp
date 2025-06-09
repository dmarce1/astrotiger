#include "Quadrature.hpp"
#include "Util.hpp"

// x=0   y=0  z=0

// x=1   y=0  z=0
// x=1   y=0  z=1
// x=1   y=1  z=0

// x=0   y=0  z=0
// x=0   y=0  z=1
// x=0   y=1  z=0
// x=0   y=0  z=2
// x=0   y=1  z=1
// x=0   y=2  z=0

// x=1   y=0  z=0
// x=1   y=0  z=1
// x=1   y=1  z=0
// x=1   y=0  z=2
// x=1   y=1  z=1
// x=1   y=2  z=0

// x=2   y=0  z=0
// x=2   y=0  z=1
// x=2   y=1  z=0
// x=2   y=0  z=2
// x=2   y=1  z=1
// x=2   y=2  z=0

template<typename Type, int M, int D>
auto nodalToModal(std::array<Type, power(M, D)> const &x) {
	static auto const Q = getLegendreQuadraturePoints(M);
	constexpr int sqrSize = power(M, D);
	constexpr int triSize = binomialCoefficient(M + D - 1, D);
	if constexpr (D > 1) {
		constexpr int dSize = binomialCoefficient(M + D - 2, D - 1);
		std::array<std::array<Type, M>, dSize> z;
		for (int i = 0; i < M; i++) {
			z[i] = nodalToModal(x + i * power(M, D - 1));
		}
		std::array<Type, triSize> w;
		int ii = 0;
		int mm = 1;
		int d = 0;
		do {
			int jj = 0;
			for (int i = 0; i < mm; i++) {
				w[ii] = 0.0;
				for (int j = 0; j < M; j++) {
					w[ii] += Q[j].w * (2 * i + 1) * 0.5 * std::legendre(i, Q[j].x) * z[j][jj];
				}
				ii++;
				jj++;
			}
			mm *= M + d;
			mm /= d + 1;
			d++;
		} while (ii < triSize);
	} else {
		std::array<Type, M> y;
		for (int n = 0; n < M; n++) {
			y[n] = 0.0;
			for (int m = 0; m < M; m++) {
				y[n] += Q[m].w * (2 * n + 1) * 0.5 * std::legendre(n, Q[m].x) * x[m];
			}
		}
		return y;
	}
}
