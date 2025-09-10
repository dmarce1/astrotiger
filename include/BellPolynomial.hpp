/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_BELLPOLYNOMIAL_HPP_
#define INCLUDE_BELLPOLYNOMIAL_HPP_

#include "Multidices.hpp"

template<int D1, int D2>
auto const& multivariateBellPolynomial(Multidices<D1> const &N, Multidices<D2> const &K) {
	using TermType = std::vector<std::pair<Multidices<D1>, Multidices<D2>>>;
	using ResultType = std::vector<TermType>;
	using MBMapType = std::unordered_map<Multidices<D1>, Multidices<D2>>;
	static thread_local std::unordered_map<Multidices<D1>, std::unordered_map<Multidices<D2>, std::shared_ptr<ResultType>>> memory;
	auto &mapN = memory[N];
	auto iterator = mapN.find(K);
	if (iterator == mapN.end()) {
		MBMapType J2K;
		std::vector<MBMapType> J2Ks;
		std::function<void(Multidices<D1>, Multidices<D1> const&, Multidices<D2> const&)> recurse;
		recurse = [&recurse, &J2K, &J2Ks](Multidices<D1> J, Multidices<D1> const &N, Multidices<D2> const &K) {
			for (Multidices<D2> Kj = 0; J * abs(Kj) <= N; Kj++) {
				auto const JabsK = J * abs(Kj);
				if ((JabsK <= N) && (Kj <= K)) {
					auto const N1 = N - JabsK;
					auto const K1 = K - Kj;
					if (abs(Kj) > 0) {
						J2K[J] = Kj;
					}
					if ((abs(N1) == 0) && (abs(K1) == 0)) {
						J2Ks.push_back(J2K);
					} else {
						bool end = false;
						auto J1 = J;
						do {
							J1++;
							if (!(abs(J1) <= abs(N))) {
								end = true;
								break;
							}
						} while (!(J1 <= N));
						if (!end) {
							recurse(J1, N1, K1);
						}
					}
				}
			}
			J2K.erase(J);
		};
		recurse(Multidices<D1>(1), N, K);
		ResultType result;
		for (auto const &j2k : J2Ks) {
			TermType term;
			for (auto const &t : j2k) {
				term.push_back(t);
			}
			result.push_back(std::move(term));
		}
		iterator = mapN.insert(std::make_pair(K, std::make_shared<ResultType>(std::move(result)))).first;
	}
	return *iterator->second;
}

#endif /* INCLUDE_BELLPOLYNOMIAL_HPP_ */
