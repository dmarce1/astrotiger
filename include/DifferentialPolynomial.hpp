#pragma once

#include <cassert>
#include <vector>
#include <unordered_map>

#include "Multidices.hpp"

//template<int O, int D1, int D2>
//struct MultivariateBellPolynomial {
//	struct Polynomial {
//		int n;
//		int k;
//		Multidices<D1> N;
//		Multidices<D2> K;
//		std::array<Multidices<D2>, N> Kj;
//		std::array<Multidices<D1>, N> J;
//	};
//	std::vector<std::vector<Polynomial>> Bnk_;
//	MultivariateBellPolynomial() {
//		Multidices<O> kj;
//		for (int k = 0; k < O; k++) {
//			indices[k] = 0;
//		}
//		auto const initKset = [](std::array<Multidices<D2>, O> &Kjs, Multidices<O> kj) {
//			for (int n = 0; n < O; n++) {
//				Kjs[n] = Multidices<D2>(kj[n]);
//			}
//		};
//		auto const nextKset = [](std::array<Multidices<D2>, O> &Kjs, Multidices<O> kj) -> bool {
//			int j = 0;
//			while (abs(++Kjs[j]) > kj[j]) {
//				Kjs[j] = Multidices<D2>(kj[j]);
//				j++;
//				if (j == O) {
//					return false;
//				}
//			}
//			return true;
//		};
//		int j;
//		do {
//			Polynomial entry;
//			entry.n = 0;
//			entry.k = 0;
//			for (int j = 0; j < O; j++) {
//				entry.n += (j + 1) * kj[j];
//				entry.k += kj[j];
//			}
//			entry.n = n;
//			std::array<Multidices<D2>, O> Kjs;
//			initKset(Kjs, kj);
//			do {
//
//			} while (nextKset(Kjs, kj));
//			j = 0;
//			do {
//				indices[j]++;
//				entry.n += j + 1;
//				if ((entry.n >= O) || (kj[j] * (j + 1) >= O)) {
//					entry.n -= kj[j] * (j + 1);
//					kj[j] = 0;
//					j++;
//				} else {
//					break;
//				}
//			} while (j < O);
//		} while (j < O);
//	}
//private:
//};

//#include <map>
//#include <hpx/synchronization/once.hpp>
//#include <unordered_map>
//
//template<int N>
//struct std::hash<Multidices<N>> {
//	size_t operator()(Multidices<N> const &indices) const {
//		constexpr size_t n = ((size_t(1) << size_t(31)) - size_t(1));
//		constexpr size_t a = 48271;
//		constexpr size_t b = 16807;
//		std::hash<int> keyGen;
//		size_t key = 42;
//		for (int i = 0; i < N; i++) {
//			key = ((a * key) % n) ^ ((b * keyGen(i)) % n);
//		}
//		return key;
//	}
//};
//
//template<int N, int D>
//using DIndex = std::array<Multidices<D>, N>;
//
//
//
//template<int N, int D>
//struct std::hash<DIndex<N, D>> {
//	size_t operator()(DIndex<N, D> const &dIndices) const {
//		constexpr size_t n = ((size_t(1) << size_t(31)) - size_t(1));
//		constexpr size_t a = 0x7FFFFFED;
//		constexpr size_t c = 0x7FFFFFC3;
//		std::hash<<Multidices<D>> kHash;
//		size_t key = 1234;
//		for (int d2 = 0; d2 < D2; d2++) {
//			for (int n = 0; n < N; n++) {
//				key ^= ((a * kHash(dIndices[n][d]) + c) % n);
//			}
//		}
//		return key;
//	}
//};
//

template<int D>
int tri2flat(std::array<int, D> n) {
	int flat = 0;
	for (int d = 0; d < D - 1; d++) {
		n[d] -= n[d + 1];
		flat += binco(n[d] + D - d - 1, D - d);
	}
	flat += n.back();
	return flat;
}

template<typename Type, int maxDegree, int varCount, int dxCount>
class DifferentialPolynomial {
	static constexpr int jetSize = binco(maxDegree + dxCount, dxCount);
	using JetType = std::array<Multidices<varCount>, jetSize>;
	struct HashKey {
		size_t operator()(JetType const &A) const {
			std::hash<int> keyGen;
			constexpr size_t n = ((size_t(1) << size_t(31)) - size_t(1));
			constexpr size_t a = 0x7FFFFFED;
			constexpr size_t c = 0x7FFFFFC3;
			size_t key = 42;
			for (int i = 0; i < jetSize; i++) {
				for (int d = 0; d < varCount; d++) {
					key ^= (a * keyGen(A[i][d]) + c) % n;
				}
			}
			return key;
		}
	};
	std::unordered_map<JetType, Type, HashKey> C_;
public:
	static DifferentialPolynomial getMonomial(int varNum = 0) {
		DifferentialPolynomial P;
		JetType I;
		I.fill(0);
		I[0][varNum] = 1;
		P.C_[I] = Type(1);
		return P;
	}
	DifferentialPolynomial() = default;
	DifferentialPolynomial(DifferentialPolynomial const&) = default;
	DifferentialPolynomial(Type const &A) {
		JetType I;
		I.fill(0);
		C_.clear();
		C_[I] = A;
	}
	DifferentialPolynomial& operator=(DifferentialPolynomial const &A) {
		C_ = A.C_;
		return *this;
	}
	DifferentialPolynomial& operator=(DifferentialPolynomial &&A) {
		C_ = std::move(A.C_);
		return *this;
	}
	DifferentialPolynomial& operator=(Type &&A) {
		JetType I;
		I.fill(0);
		C_.clear();
		C_[I] = std::move(A);
		return *this;
	}
	DifferentialPolynomial& operator+=(DifferentialPolynomial const &A) {
		*this = *this + A;
		return *this;
	}
	DifferentialPolynomial& operator-=(DifferentialPolynomial const &A) {
		*this = *this - A;
		return *this;
	}
	DifferentialPolynomial& operator*=(DifferentialPolynomial const &A) {
		*this = *this * A;
		return *this;
	}
	DifferentialPolynomial& operator*=(Type const &A) {
		*this = *this * A;
		return *this;
	}
	Type& operator()(JetType const &I) {
		auto it = C_.find(I);
		if (it == C_.end()) {
			C_[I] = Type(0);
		}
		return C_[I];
	}
	Type operator()(JetType const &I) const {
		auto it = C_.find(I);
		if (it == C_.end()) {
			return Type(0);
		}
		return C_[I];
	}
	friend std::ostream& operator<<(std::ostream &os, DifferentialPolynomial const &P) {
		std::vector<std::pair<JetType, Type>> terms(P.C_.begin(), P.C_.end());
		std::sort(terms.begin(), terms.end(), [](std::pair<JetType, Type> const &A, std::pair<JetType, Type> const &B) {
			for (int v = 0; v < varCount; v++) {
				for (int i = 0; i < jetSize; i++) {
					auto const a = A.first[i][v];
					auto const b = B.first[i][v];
					if (a != b) {
						return (b < a);
					}
				}
			}
			return false;
		});
		Multidices<varCount> zeroIndex;
		JetType zeroJet;
		zeroJet.fill(0);
		bool first = true;
		for (auto const &term : terms) {
			Type coefficient = term.second;
			if (!first) {
				if (coefficient > 0) {
					os << print2string("+ ");
				} else {
					coefficient = -coefficient;
					os << print2string("- ");
				}
			}
			if ((term.second != Type(1)) || (term.first == zeroJet)) {
				if constexpr (std::is_floating_point<Type>::value) {
					os << print2string("(%e) ", coefficient);
				} else {
					os << print2string("%i ", coefficient);
				}
			}
			std::array<int, dxCount> jetIndex;
			jetIndex.fill(0);
			for (int flat = 0; flat < jetSize; flat++) {
				if (term.first[flat] != zeroIndex) {
					for (int vdim = 0; vdim < varCount; vdim++) {
						int const pwr = term.first[flat][vdim];
						if (pwr > 0) {
							int const N = std::accumulate(jetIndex.begin(), jetIndex.end(), 0);
							if (N != 0) {
								os << print2string("(");
							}
							if (N > 0) {
								os << print2string("d");
								if (N > 1) {
									os << print2string("%i", N);
								}
							}
							char name = 'f' + vdim;
							os << print2string("%c", name);
							if (N > 0) {
								os << print2string("/");
								for (int di = 0; di < dxCount; di++) {
									int const n = jetIndex[di];
									if (n > 0) {
										os << print2string("d%c", 'x' + di);
										if (n > 1) {
											os << print2string("%i", n);
										}
									}
								}
							}
							if (N != 0) {
								os << print2string(")");
							}
							if (pwr > 1) {
								os << print2string("^%i", pwr);
							}
						}
						os << print2string(" ");
					}
				}
				int i = dxCount - 1;
				while (1) {
					if (i == 0) {
						jetIndex[0]++;
						break;
					}
					if (++jetIndex[i] > jetIndex[i - 1]) {
						jetIndex[i] = 0;
						i--;
					} else {
						break;
					}
				}
			}
			first = false;
		}
		return os;
	}
	friend DifferentialPolynomial operator+(DifferentialPolynomial A) {
		return A;
	}
	friend DifferentialPolynomial operator-(DifferentialPolynomial A) {
		for (auto it = A.C_.begin(); it != A.C_.end(); it++) {
			it->second = -it->second;
		}
		return A;
	}
	friend DifferentialPolynomial operator+(DifferentialPolynomial A, DifferentialPolynomial const &B) {
		for (auto it = B.C_.begin(); it != B.C_.end(); it++) {
			A.C_[it->first] += it->second;
		}
		return A;
	}
	friend DifferentialPolynomial operator+(DifferentialPolynomial A, Type const &B) {
		JetType I;
		I.fill(0);
		A.C_[I] += B;
		return A;
	}
	friend DifferentialPolynomial operator+(Type const &A, DifferentialPolynomial B) {
		return B + A;
	}
	friend DifferentialPolynomial operator-(DifferentialPolynomial A, DifferentialPolynomial const &B) {
		return A + (-B);
	}
	friend DifferentialPolynomial operator-(DifferentialPolynomial A, Type const &B) {
		return A + (-B);
	}
	friend DifferentialPolynomial operator-(Type const &A, DifferentialPolynomial B) {
		return A + (-B);
	}
	friend DifferentialPolynomial operator*(DifferentialPolynomial const &A, DifferentialPolynomial const &B) {
		DifferentialPolynomial C;
		for (auto itA = A.C_.begin(); itA != A.C_.end(); itA++) {
			for (auto itB = B.C_.begin(); itB != B.C_.end(); itB++) {
				JetType iC = itA->first + itB->first;
				C.C_[iC] += itA->second * itB->second;
			}
		}
		return C;
	}
	friend DifferentialPolynomial operator*(DifferentialPolynomial A, Type const &B) {
		for (auto itA = A.C_.begin(); itA != A.C_.end(); itA++) {
			itA->second *= B;
		}
		return A;
	}
	friend DifferentialPolynomial operator*(Type const &A, DifferentialPolynomial B) {
		return B * A;
	}
	friend DifferentialPolynomial differentiate(DifferentialPolynomial const &A, int ddim = 0) {
		DifferentialPolynomial D;
		Multidices<varCount> zeroIndex;
		std::array<Type, varCount> zeroValue;
		zeroValue.fill(Type(0));
		for (auto itA = A.C_.begin(); itA != A.C_.end(); itA++) {
			JetType thisJet = itA->first;
			std::array<int, dxCount> jetIndex;
			jetIndex.fill(0);
			for (int flatIndex = 0; flatIndex < jetSize; flatIndex++) {
				if (thisJet[flatIndex] != zeroIndex) {
					JetType newJet = thisJet;
					auto newCo = itA->second;
					auto newIndex = jetIndex;
					newIndex[ddim]++;
					for (int vdim = 0; vdim < varCount; vdim++) {
						int const pwr = newJet[flatIndex][vdim];
						newJet[flatIndex][vdim]--;
						newJet[tri2flat<dxCount>(newIndex)][vdim]++;
						newCo *= pwr;
					}
					D.C_[newJet] += newCo;
				}
				int i = dxCount - 1;
				while (1) {
					if (i == 0) {
						jetIndex[0]++;
						break;
					}
					if (++jetIndex[i] > jetIndex[i - 1]) {
						jetIndex[i] = 0;
						i--;
					} else {
						break;
					}
				}
			}
		}
		return D;
	}
};
