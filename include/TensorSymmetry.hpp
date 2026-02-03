#pragma once
#include "Matrix.hpp"
#include "SparseMatrix.hpp"
#include "Young.hpp"

struct SymmetryElement {
	Index row;
	Index col;
	double value;
	friend std::ostream &operator<<(std::ostream &os, SymmetryElement const &ele) {
		os << "(" << ele.row << ", " << ele.col << "): " << ele.value;
		return os;
	}
};

namespace detail {
template <std::integral auto O, std::integral auto D, IntegralArray auto... Ts>
constexpr auto countDetracerConditions() {
	int count = 0;
	((count +=
	  [](unsigned n) {
		  return n * (n - 1) / 2;
	  }(Ts.size())),
	 ...);
	return ipow(D, O - 2) * count;
}

template <std::integral auto O, std::integral auto D, std::integral auto Ndof, IntegralArray auto... Ts>
constexpr auto detracer(MatrixType auto Q) {
	constexpr int N = Q.rowCount();
	constexpr int R = Q.colCount();
	constexpr int L = countDetracerConditions<O, D, Ts...>();
	SparseMatrix<double> Ω(L, N);
	using IndexType = Indices<O, D>;
	int n = 0;
	(([&n, &Ω](IntegralArray auto T) {
		 for (unsigned i1 = 0; i1 < T.size(); i1++) {
			 auto const k1 = T[i1];
			 for (unsigned i2 = i1 + 1; i2 < T.size(); i2++) {
				 auto const k2 = T[i2];
				 IndexType ms{};
				 int const count = ipow(D, O - 2);
				 for (int ci = 0; ci < count; ci++) {
					 for (ms[k1] = ms[k2] = 0; ms[k1] < D; ms[k1] = ++ms[k2]) {
						 auto const m = Index(ms);
						 assert(n <= L);
						 assert(m <= N);
						 Ω(n, m) = 1;
					 }
					 if (ci + 1 < count) {
						 for (unsigned k = O - 1; k >= 0; k--) {
							 if (k == k1) continue;
							 if (k == k2) continue;
							 if (++ms[k] != D) break;
							 ms[k] = 0;
						 }
					 }
					 n++;
				 }
			 }
		 }
	 }(Ts)),
	 ...);
	auto V = (Ω * Q).reducedRowEchelonForm();
	if constexpr (Ndof < 0) {
		return R - V.rowCount();
	} else {
		SparseMatrix<double> A(N, Ndof);
		std::bitset<R> isPivot{};
		for (int r = 0; r < V.rowCount(); r++) {
			isPivot[V[r].begin()->first] = true;
		}
		auto W = SparseMatrix<double>(Q);
		Matrix<double, N, Ndof> X{};
		int k = 0;
		for (int m = 0; m < R; m++) {
			if (isPivot[m]) {
				int j;
				for (j = 0; W(j, m) != 1.0; j++) {
				}
				for (int n = 0; n < N; n++) {
					if (W(n, m) != 0.0) {
						W[n] -= W(n, m) * W[j];
					}
				}
			} else {
				k++;
			}
		}
		k = 0;
		for (int m = 0; m < R; m++) {
			if (!isPivot[m]) {
				for (int n = 0; n < N; n++) {
					X(n, k) = W(n, m);
				}
				k++;
			}
		}
		return X;
	}
}

} // namespace detail

// template <std::integral auto O, std::integral auto D, IntegralArray auto... Ts>
// constexpr auto detracer(MatrixType auto Q) {
//	constexpr auto N1 = ipow(D, O);
//	constexpr auto N2 = detail::detracer<O, D, -1, Ts...>(Q);
//	return detail::detracer<O, D, N2, Ts...>(Q);
// }

template <std::integral auto D, IntegerPartitionType auto Λ, IntegralArray auto... Ts>
constexpr auto symmetrizer() {
	if constexpr (sizeof...(Ts)) {
		constexpr int O = Λ.size();
		constexpr auto Q = symmetrizer<D, Λ>();
		constexpr auto N1 = ipow(D, O);
		constexpr auto N2 = detail::detracer<O, D, -1, Ts...>(Q);
		return detail::detracer<O, D, N2, Ts...>(Q);
	} else {
		constexpr auto Y = createYoungTableau(Λ);
		constexpr auto basis = genSemistandardTableau<Λ, D>();
		constexpr auto ps = genYoungPermutations(Λ);
		constexpr int R = basis.size();
		constexpr int O = Y.size();
		constexpr int N = ipow(D, O);
		using IndexType = Indices<O, D>;
		Matrix<double, R, N> Ψ{};
		int i = 0;
		std::bitset<N> isBasis{};
		for (auto const &b : basis) {
			isBasis[b.template index<D>()] = true;
		}
		for (auto ns = IndexType::ibegin(); ns != IndexType ::iend(); ns++) {
			auto const n = Index(ns);
			if (isBasis[n]) {
				for (auto p : ps) {
					auto const ms = apply<O>(p.second, ns);
					auto const m = Index(ms);
					Ψ(i, m) += p.first;
				}
				for (auto p : ps) {
					auto const ms = apply<O>(p.second, ns);
					auto const m = Index(ms);
					auto &v = Ψ(i, m);
					v = std::max(std::min(v, 1.0), -1.0);
				}
				i++;
			}
		}
		return transpose(Ψ);
	}
}
