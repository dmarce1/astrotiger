#pragma once

#include "Indices.hpp"
#include "IntegerPartition.hpp"
#include "LinearCombination.hpp"
#include "Matrix.hpp"
#include "TensorSymmetry.hpp"
#include "SparseMatrix.hpp"
#include "Young.hpp"
#include <type_traits>

//template<int R, int D>
//class TensorSymmetry {
//	constexpr int static N = ipow(D, R);
//	std::array<std::pair<std::pair<int, int>, Rational>, N> transform_;
//
//public:
//	using LiteralType = std::array<std::array<Rational, N>, N>;
//	template<size_t N>
//	constexpr TensorSymmetry(std::array<std::pair<std::pair<int, int>, Rational>, N> const &t) :
//			transform_ { } {
//		transform_ = t;
//	}
//	constexpr auto unpack() const {
//		auto tmp = pack();
//		return rightPseudoinverse(tmp);
//	}
//	constexpr auto pack() const {
//		auto tmp = SparseMatrix<Rational>(N, N, transform_);
//		rankReduce(tmp);
//		return tmp;
//	}
//	constexpr auto packedSize() const {
//		for (int n = 0; n < N; n++) {
//			bool z = true;
//			for (int m = 0; m < N; m++) {
//				if (transform_(n, m) != Rational(0)) {
//					z = false;
//					break;
//				}
//			}
//			if (z) {
//				return n;
//			}
//		}
//		return N;
//	}
//
//	constexpr auto unpackedSize() const {
//		return N;
//	}
//};

template <IntegerPartitionType auto Λ>
constexpr auto youngSymmetrizer() {
	constexpr auto R = Λ.size();
	constexpr YoungTableau ν(Λ);
	using TableauType = std::remove_cvref_t<decltype(ν)>;
	constexpr auto nRow = ν.rowPermutationCount();
	constexpr auto nCol = ν.colPermutationCount();
	constexpr auto rows = [R, ν]<Permutation<R> σ, auto rowRecur, int iR = 0>() {
		constexpr auto cols = [R]<Permutation<R> τ, auto colRecur, int iC = 0>() {
			constexpr auto τNext = TableauType(τ).nextColPermutation();
			constexpr auto ele = std::array<std::pair<Rational, Permutation<R>>, 1> {
				std::pair<Rational, Permutation<R>> {Rational(parity<R>(τ) * parity<R>(σ)), τ}};
			if constexpr (iC + 1 == nCol) {
				return ele;
			} else {
				constexpr auto next = colRecur.template operator()<(Permutation<R>)τNext, colRecur, iC + 1>();
				return plus<Rational, Permutation<R>, ele.size(), next.size(), ele, next>();
			}
		};
		constexpr auto σNext = TableauType(σ).nextRowPermutation();
		constexpr auto ele = cols.template operator()<(Permutation<R>)σ, cols>();
		if constexpr (iR + 1 == nRow) {
			return ele;
		} else {
			constexpr auto next = rowRecur.template operator()<(Permutation<R>)σNext, rowRecur, iR + 1>();
			return plus<Rational, Permutation<R>, ele.size(), next.size(), ele, next>();
		}
	};
	auto const A = rows.template operator()<(Permutation<R>)ν, rows>();
	return A;
}

//template <IntegerPartitionType auto Λ, std::integral auto D>
//constexpr auto createTensorSymmetrySize() {
//	constexpr auto R = Λ.size();
//	constexpr auto N = ipow(D, R);
//	using IndexType = Indices<R, D>;
//	constexpr auto ySym = youngSymmetrizer<Λ>();
//	SparseMatrix<Rational> A{};
//	std::map<int, std::map<int, Rational>> data{};
//	for (unsigned i = 0; i < ySym.size(); i++) {
//		for (auto ns = IndexType::ibegin(); ns != IndexType::iend(); ns++) {
//			auto const c = ySym[i].first;
//			auto const p = ySym[i].second;
//			auto const ms = apply<R>(p, ns);
//			A(Index(ns), Index(ms)) += c;
//		}
//	}
//	return A.density();
//}
//

template <IntegerPartitionType auto Λ, std::integral auto D, bool sizeOnly = false>
constexpr auto createTensorSymmetry() {
	constexpr auto R = Λ.size();
	constexpr auto N = ipow(D, R);
	using IndexType = Indices<R, D>;
	constexpr auto ySym = youngSymmetrizer<Λ>();
	SparseMatrix<Rational> A(N);
	for (unsigned i = 0; i < ySym.size(); i++) {
		for (auto ns = IndexType::ibegin(); ns != IndexType::iend(); ns++) {
			auto const c = ySym[i].first;
			auto const p = ySym[i].second;
			auto const ms = apply<R>(p, ns);
			A(Index(ns), Index(ms)) += c;
		}
	}
	auto const rank = rankReduce(A);
	auto const iA = rightPseudoinverse(A);
	//	auto const iA = A;
	if constexpr(sizeOnly) {
		return std::pair(A.density(), iA.density());
	} else {
		constexpr auto sz = createTensorSymmetry<Λ, D, true>();
		constexpr auto Asz = sz.first;
		constexpr auto iAsz = sz.second;
		std::array<std::pair<std::pair<int, int>, Rational>, Asz> Alit {};
		std::array<std::pair<std::pair<int, int>, Rational>, iAsz> iAlit {};
		auto const Blit = A.literal();
		auto const iBlit = iA.literal();
		for(int i = 0; i < Asz; i++) {
			Alit[i] = Blit[i];
		}
		for(int i = 0; i < iAsz; i++) {
			iAlit[i] = iBlit[i];
		}
		return std::tuple(A.rowCount(), Alit, iAlit);
	}
}

//template <std::integral auto R1, std::integral auto R2, std::integral auto D>
//constexpr auto tensorSymmetryProduct(Matrix<Rational, pow(D, R1)> const &A, Matrix<Rational, pow(D, R1)> const &B) {
//	return kroneckerProduct(A, B);
//}
//
//template <std::integral auto R, std::integral auto D, std::integral auto I, std::integral auto J>
//constexpr auto tensorSymmetryContraction(Matrix<Rational, ipow(D, R)> const &A) {
//	constexpr int N1 = pow(D, R);
//	constexpr int N2 = N1 / (D * D);
//	Matrix<Rational, N2, N2> B{};
//	using Index1 = Indices<R, D>;
//	for (auto I1 = Index1::ibegin(); I1 != Index1::iend(); I1++) {
//		auto const I2 = I1.template contract<I, J>();
//		if (I1[I] == I1[J]) {
//			for (auto J1 = Index1::ibegin(); J1 != Index1::iend(); J1++) {
//				if (J1[I] == J1[J]) {
//					auto const J2 = J1.template contract<I, J>();
//					B(Index(I2), Index(J2)) += A(Index(I1), Index(J1));
//				}
//			}
//		}
//	}
//	return B;
//}

// template <std::integral auto R, std::integral auto D, std::integral auto Np, std::integral auto Nt>
// struct TensorSymmetry {
//	static constexpr int N = ipow(D, R);
//	using PType = LinearCombination<Rational, Permutation<R>, Np>;
//	using TType = std::array<std::pair<int, int>, Nt>;
//	PType Ps;
//	TType Ts;
//	constexpr int rank() const {
//		return R;
//	}
//	constexpr TensorSymmetry() {
//	}
//	constexpr TensorSymmetry(PType p, TType t) :
//	Ps(p), Ts(t) {
//	}
//	constexpr auto weight() const {
//		using std::abs;
//		Rational wt = Rational(0);
//		for (unsigned i = 0; i < Np; i++) {
//			wt += abs(Ps[i].first);
//		}
//		return wt;
//	}
//	constexpr auto toMatrix() const {
//		constexpr int N = ipow(D, R);
//		Matrix<Rational, N> A;
//		for (unsigned n = 0; n < N; n++) {
//			for (unsigned m = 0; m < N; m++) {
//				A(n, m) = Rational(0);
//			}
//		}
//		using IndexType = Indices<R, D>;
//		for (unsigned i = 0; i < Np; i++) {
//			for (auto ns = IndexType::ibegin(); ns != IndexType::iend(); ns++) {
//				auto const c = Ps[i].first;
//				auto const p = Ps[i].second;
//				auto const ms = apply<R>(p, ns);
//				A(Index(ns), Index(ms)) += c;
//			}
//		}
//		constexpr int Nc = Nt * ipow(D, R - 2);
//		if constexpr(Nt) {
//			Matrix<Rational, Nc, N> C {};
//			int ci = 0;
//			for (unsigned i = 0; i < Nt; i++) {
//				int const iA = Ts[i].first;
//				int const iB = Ts[i].second;
//				for (auto ns = IndexType::ibegin(); ns != IndexType::iend(); ns++) {
//					if (ns[iB] == ns[iA] && ns[iB] == 0) {
//						for (int k = 0; k < D; k++) {
//							auto ms = ns;
//							ms[iA] = ms[iB] = k;
//							C(ci, Index(ms)) = Rational(1) / Rational(D);
//						}
//						ci++;
//					}
//				}
//			}
//			Matrix<Rational, N> Pc = Matrix<Rational, N>::identity();
//			auto const trC = transpose(C);
//			auto const CtrC = C * trC;
//			auto const Cinv = trC * inverse(CtrC);
//			Pc = Pc - Cinv * C;
//			A = Pc * A;
//		}
//		return A;
//	};
// };
//
// template <auto sym>
// consteval auto genTransform() {
//	constexpr int N = sym.N;
//	constexpr auto A = sym.toMatrix();
//	constexpr auto B = rankReduce<Rational, N, N, A.literal()>();
//	static_assert(B.rowCount());
//	static_assert(B.columnCount());
//	return Rational(sym.weight()) * pseudoinverse(B);
// };
//
//
////template <IntegerPartitionType auto Λ>
////struct SpechtModule {
////	static constexpr auto tableau() {
////		YoungTableau t(Λ);
////		return t;
////	}
////	using TableauType = std::remove_cvref_t<decltype(tableau())>;
////	static constexpr auto size() {
////		return YoungTableau(Λ).standardCount();
////	}
////	template <int dim>
////	constexpr auto access(TableauType P) const {
////		auto eT = YoungSymmetrizer<Λ>();
////		if constexpr (dim == 0) {
////			constexpr int R = Λ.size();
////			for (auto &σ : eT) {
////				σ.second = apply<R>(P, σ.second);
////			}
////			return eT;
////		} else {
////			do {
////				std::next_permutation(P.begin(), P.end());
////			}while (!P.isStandard());
////			return access<dim - 1>(P);
////		}
////	}
////	template <int dim>
////	constexpr auto operator()() const {
////		auto constexpr tab = TableauType {};
////		return access<dim>(tab);
////	}
////};
//
//
////template <std::integral auto R, std::integral auto D, std::pair<int, int>... Ts>
////constexpr auto createTensorSymmetry() {
////	constexpr std::array<std::pair<Rational, Permutation<R>>, 1> eT = {
////		std::pair<Rational, Permutation<R>>(1, identityPermutation<R>())};
////	constexpr std::array<std::pair<int, int>, sizeof...(Ts)> tr = {Ts...};
////	constexpr TensorSymmetry<R, D, eT.size(), sizeof...(Ts)> sym(eT, tr);
////	return sym;
////}
//
////namespace detail {
////template <int R1, int R2, int N1, int N2, LinearCombination<Rational, Permutation<R1>, N1> linCo1,
////LinearCombination<Rational, Permutation<R2>, N2> linCo2>
////consteval auto rawDirectProduct() {
////	constexpr int N3 = N1 * N2;
////	constexpr int R3 = R1 + R2;
////	LinearCombination<Rational, Permutation<R3>, N3> linCo3;
////	for (int n1 = 0; n1 < N1; n1++) {
////		for (int n2 = 0; n2 < N2; n2++) {
////			linCo3[n1 * N2 + n2].first = linCo1[n1].first * linCo2[n2].first;
////			linCo3[n1 * N2 + n2].second = concatenate<R1, R2>(linCo1[n1].second, linCo2[n2].second);
////		}
////	}
////	return linCo3;
////}
////}
// // namespace detail
////
////template <int R1, int R2, int N1, int N2, LinearCombination<Rational, Permutation<R1>, N1> linCo1,
////LinearCombination<Rational, Permutation<R1>, N1> linCo2>
////consteval auto directProduct() {
////constexpr int N3 = N1 + N2;
////constexpr int R3 = R1 + R2;
////constexpr auto linCo3 = detail::rawDirectProduct<R1, R2, N1, N2, linCo1, linCo2>();
////return canonicalize<Rational, Permutation<R3>, N3, linCo3>();
////}
