/*
 * LinearGravity.hpp
 *
 *  Created on: Feb 22, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_LINEARGRAVITY_HPP_
#define INCLUDE_LINEARGRAVITY_HPP_

#include <functional>

#include <silo.h>

#include "GridAttributes.hpp"
#include "Real.hpp"
#include "Vector.hpp"
#include "Kernels.hpp"

using namespace Math;

template<typename T>
SymmetricMatrix<T, DIM4> minkowskiMetric() {
	static constexpr T zero = T(0), one = T(1);
	using return_type = SymmetricMatrix<T, DIM4>;
	return_type nu = return_type(zero);
	nu[0, 0] = -one;
	nu[1, 1] = +one;
	nu[2, 2] = +one;
	nu[3, 3] = +one;
	return nu;
}

static void siloErrorHandler(char *errorString_) {
	std::string const errorString(errorString_);
	std::cout << "SILO returned an error." << "\n";
	std::cout << errorString << "\n";
	std::cout << "Aborting..." << "\n";
	abort();
}

template<typename T>
SymmetricMatrix<T, DIM4> reverseTrace(SymmetricMatrix<T, DIM4> Tij) {
	static constexpr T zero = T(0), half = T(0.5);
	static auto const nu = minkowskiMetric<T>();
	T trT = zero;
	for (int k = 0; k < DIM4; k++) {
		trT += nu[k, k] * Tij[k, k];
	}
	for (int k = 0; k < DIM4; k++) {
		Tij[k, k] -= half * nu[k, k] * trT;
	}
	return Tij;
}

template<typename RealType>
struct LinearGravity {
	using real_type = RealType;
	using grid_type = std::valarray<real_type>;
	static constexpr real_type zero = real_type(0), half = real_type(0.5), one = real_type(1);
	static constexpr int BW = 2;
	LinearGravity(int N_) :
			gAttr(Vector<int, NDIM>( { N_, N_, N_ }), BW), Ni(gAttr.intSizes[0]), Ne(gAttr.extSizes[0]), sizes(gAttr.intSizes), strides(gAttr.extStrides), N3(
					Ne * Ne * Ne), start(Ne * (Ne * BW + BW) + BW) {
		for (int i = 0; i < DIM4; i++) {
			for (int j = 0; j <= i; j++) {
				for (int k = 0; k < DIM4; k++) {
					U[k][i, j] = grid_type(zero, N3);
				}
			}
		}
		for (int k = 0; k < NDIM; k++) {
			X[k] = grid_type(N3);
			Vector<size_t, NDIM> x;
			for (x[XDIM] = 0; x[XDIM] != Ne; x[XDIM]++) {
				for (x[YDIM] = 0; x[YDIM] != Ne; x[YDIM]++) {
					for (x[ZDIM] = 0; x[ZDIM] != Ne; x[ZDIM]++) {
						int index = 0;
						for (int n = 0; n < NDIM; n++) {
							index += strides[n] * x[n];
						}
						X[k][index] = (real_type(((int) x[k]) - BW) + half);
					}
				}
			}
		}
		return;
	}
	void step(SymmetricMatrix<grid_type, DIM4> const &T, real_type dt) {
		auto const slopeLimiter = [](real_type const &a, real_type const &b) -> real_type {
			using namespace std;
			return (copysign(half, a) + copysign(half, b)) * (min(abs(a), abs(b)));
		};
		constexpr real_type G = real_type(-8.0 * M_PI);
		real_type const dx = one / real_type(Ni);
		real_type const dxinv = real_type(Ni);
		Vector<grid_type, DIM4> D(grid_type(zero, N3));
		auto const &strides = gAttr.extStrides;
		auto const dV = gAttr.d3R;
		auto const N3 = gAttr.extSize;
		grid_type UxR = grid_type(N3);
		grid_type UtR = grid_type(N3);
		grid_type UxL = grid_type(N3);
		grid_type UtL = grid_type(N3);
		grid_type Fx = grid_type(N3);
		grid_type Ft = grid_type(N3);
		decltype(U) U0 = U;
		decltype(U) dU(grid_type(zero, N3));
		for (int rkSubstep = 0; rkSubstep < 2; rkSubstep++) {
			for (int k = 0; k < DIM4; k++) {
				for (int m = 0; m < DIM4; m++) {
					for (int j = 0; j <= m; j++) {
						dU[k][m, j] = (G * dV) * T[m, j];
					}
				}
			}
			for (int k = 1; k < DIM4; k++) {
				size_t const dk = strides[k - 1];
				for (int m = 0; m < DIM4; m++) {
					for (int j = 0; j <= m; j++) {
						for (size_t i0 = dk; i0 < N3 - dk; i0++) {
							size_t const ip = i0 + dk;
							size_t const im = i0 - dk;
							real_type const uxp = U[k][m, j][ip];
							real_type const utp = U[0][m, j][ip];
							real_type const ux0 = U[k][m, j][i0];
							real_type const ut0 = U[0][m, j][i0];
							real_type const uxm = U[k][m, j][im];
							real_type const utm = U[0][m, j][im];
							real_type const slp_x = slopeLimiter(uxp - ux0, ux0 - uxm);
							real_type const slp_t = slopeLimiter(utp - ut0, ut0 - utm);
							UxR[i0] = ux0 - half * slp_x;
							UtR[i0] = ut0 - half * slp_t;
							UxL[ip] = ux0 + half * slp_x;
							UtL[ip] = ut0 + half * slp_t;
						}
						for (size_t i0 = 2 * dk; i0 < N3 - dk; i0++) {
							Fx[i0] = half * (-UtL[i0] - UtR[i0]) - half * (UxR[i0] - UxL[i0]);
							Ft[i0] = half * (-UxL[i0] - UxR[i0]) - half * (UtR[i0] - UtL[i0]);
						}
						for (size_t i0 = 2 * dk; i0 < N3 - 2 * dk; i0++) {
							size_t const ip = i0 + dk;
							size_t const im = i0;
							dU[k][m, j][i0] -= (Fx[ip] - Fx[im]) * dxinv;
							dU[0][m, j][i0] -= (Ft[ip] - Ft[im]) * dxinv;
						}
					}
				}
			}
			for (int k = 0; k < DIM4; k++) {
				size_t const dk = strides[k - 1];
				for (int m = 0; m < DIM4; m++) {
					for (int j = 0; j <= m; j++) {
						for (size_t i0 = dk; i0 < N3 - dk; i0++) {
							U[k][m, j][i0] += dt * dU[k][m, j][i0];
						}
						if (rkSubstep == 1) {
							for (size_t i0 = dk; i0 < N3 - dk; i0++) {
								U[k][m, j][i0] += half * (U0[k][m, j][i0] - U[k][m, j][i0]);
							}
						}
					}
				}
			}
		}
	}
	void output(DBfile *db, DBoptlist *optList, SymmetricMatrix<grid_type, DIM4> const *Tptr = nullptr) const {
		using namespace Math;
		static constexpr int silo_data_type = DB_DOUBLE;
		char const *const coordnames[NDIM] = { "x", "y", "z" };
		Vector<std::vector<real_type>, NDIM> xCoordinates;
		for (int dim = 0; dim < NDIM; dim++) {
			xCoordinates[dim].resize(Ni + 1);
			for (int n = 0; n <= int(Ni); n++) {
				xCoordinates[dim][n] = real_type(n) / real_type(Ni);
			}
		}
		void const *const coords[NDIM] = { xCoordinates[XDIM].data(), xCoordinates[YDIM].data(), xCoordinates[ZDIM].data() };
		int const dims2[NDIM] = { (int) Ni, (int) Ni, (int) Ni };
		int dims1[NDIM] = { dims2[0] + 1, dims2[1] + 1, dims2[2] + 1 };
		DBPutQuadmesh(db, "quadMesh", coordnames, coords, dims1, NDIM, silo_data_type, DB_COLLINEAR, optList);
		auto const gSlice = std::gslice(start, sizes, strides);
		for (int i = 0; i < DIM4; i++) {
			for (int j = 0; j <= i; j++) {
				for (int k = 0; k < DIM4; k++) {
					std::valarray<real_type> const u = U[k][i, j][gSlice];
					std::string fieldname = std::string("g") + std::to_string(i) + std::to_string(j) + std::string("_") + std::to_string(k);
					DBPutQuadvar1(db, fieldname.c_str(), "quadMesh", std::begin(u), dims2, NDIM, NULL, 0, silo_data_type, DB_ZONECENT, optList);
				}
			}
		}
		if (Tptr) {
			auto const &T = *Tptr;
			for (int i = 0; i < DIM4; i++) {
				for (int j = 0; j <= i; j++) {
					std::valarray<real_type> const u = T[i, j][gSlice];
					std::string fieldname = std::string("T") + std::to_string(i) + std::to_string(j);
					DBPutQuadvar1(db, fieldname.c_str(), "quadMesh", std::begin(u), dims2, NDIM, NULL, 0, silo_data_type, DB_ZONECENT, optList);
				}
			}
		}
	}
	void enforceBoundaryConditions() {
		for (int n = 0; n < 2 * NDIM; n++) {
			for (int i = 0; i < DIM4; i++) {
				for (int j = 0; j <= i; j++) {
					for (int k = 0; k < DIM4; k++) {
						U[k][i, j][gAttr.dstBoundarySlices[n]] = U[k][i, j][gAttr.srcBoundarySlices[n]];
					}
				}
			}
		}
	}
	auto metricDerivatives(real_type x0, real_type y0, real_type z0) const {
		using return_type = Vector<SymmetricMatrix<real_type, DIM4>, DIM4>;
		return_type u = return_type(zero);
		int const i0 = x0 * real_type(Ni) + BW;
		int const j0 = y0 * real_type(Ni) + BW;
		int const k0 = z0 * real_type(Ni) + BW;
		Vector<int, NDIM> originIndices( { i0, j0, k0 });
		Vector<int, NDIM> indices;
		int &i = indices[XDIM];
		int &j = indices[YDIM];
		int &k = indices[ZDIM];
		Vector<SymmetricMatrix<real_type, DIM4>, DIM4> V;
		Vector<real_type, NDIM> const X0( { x0, y0, z0 });
		size_t const start = (i0 - 1) * strides[XDIM] + (j0 - 1) * strides[YDIM] + (k0 - 1) * strides[ZDIM];
		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				for (k = 0; k < 3; k++) {
					size_t index = start;
					for (int dim = 0; dim < NDIM; dim++) {
						index += indices[dim] * strides[dim];
					}
					real_type weight = one;
					for (int dim = 0; dim < NDIM; dim++) {
						weight *= kernelTSC(abs(X[dim][index]) - (real_type(Ni) * X0[dim]));
					}
					for (int l = 0; l < DIM4; l++) {
						for (int m = 0; m < DIM4; m++) {
							for (int n = 0; n <= m; n++) {
								u[l][m, n] += weight * U[l][m, n][index];
							}
						}
					}
				}
			}
		}
		for (int k = 0; k < DIM4; k++) {
			u[k] = reverseTrace(u[k]);
		}
		return u;
	}
	auto getStateVars(int xi, int yi, int zi) const {
		Vector<SymmetricMatrix<real_type, DIM4>, DIM4> v;
		int const index = xi * gAttr.extSizes[XDIM] + yi * gAttr.extSizes[YDIM] + zi * gAttr.extSizes[ZDIM];
		for (int k = 0; k < DIM4; k++) {
			for (int i = 0; i < DIM4; i++) {
				for (int j = 0; j <= i; j++) {
					v[k][i, j] = U[k][i, j][index];
				}
			}
		}
		return v;
	}
private:
	GridAttributes<real_type> gAttr;
	size_t const &Ni;
	size_t const &Ne;
	std::valarray<size_t> const &sizes;
	std::valarray<size_t> const &strides;
	size_t N3;
	size_t start;
	Vector<grid_type, NDIM> X;
	Vector<SymmetricMatrix<grid_type, DIM4>, DIM4> U;
}
;

#endif /* INCLUDE_LINEARGRAVITY_HPP_ */
