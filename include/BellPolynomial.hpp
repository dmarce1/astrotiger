/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_BELLPOLYNOMIAL_HPP_
#define INCLUDE_BELLPOLYNOMIAL_HPP_

#include "Multidices.hpp"

template<int D1, int D2>
auto const& multivariateBellPolynomial(Multidices<D1> const &N, size_t const &K) {
	using TermType = std::vector<std::pair<Multidices<D1>, size_t>>;
	using ResultType = std::vector<TermType>;
	using MBMapType = std::unordered_map<Multidices<D1>, size_t>;
	static thread_local std::unordered_map<Multidices<D1>, std::unordered_map<size_t, std::shared_ptr<ResultType>>> memory;
	auto &mapN = memory[N];
	auto iterator = mapN.find(K);
	if (iterator == mapN.end()) {
		MBMapType J2k;
		std::vector<MBMapType> J2ks;
		std::function<void(Multidices<D1>, Multidices<D1> const&, size_t const&)> recurse;
		recurse = [&recurse, &J2k, &J2ks](Multidices<D1> J, Multidices<D1> const &N, int k) {
			for (int kj = 0; J * kj <= N; kj++) {
				auto const JabsK = J * kj;
				if ((JabsK <= N) && (kj <= k)) {
					auto const N1 = N - JabsK;
					auto const k1 = k - kj;
					if (kj > 0) {
						J2k[J] = kj;
					}
					if ((abs(N1) == 0) && (k1 == 0)) {
						J2ks.push_back(J2k);
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
							recurse(J1, N1, k1);
						}
					}
				}
			}
			J2k.erase(J);
		};
		recurse(Multidices<D1>(1), N, K);
		ResultType result;
		for (auto const &j2k : J2ks) {
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
