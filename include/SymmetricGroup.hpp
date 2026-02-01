/*
 * SymmetricGroup.hpp
 *
 *  Created on: Dec 26, 2025
 *      Author: dmarce1
 */

#pragma once

#include "Approximate.hpp"
#include "Math.hpp"
#include "Matrix.hpp"
#include "Indices.hpp"
#include "Permutation.hpp"
#include "Rational.hpp"

#include <algorithm>
#include <numeric>
#include <utility>

template<int N>
struct SymmetricGroup {
	static constexpr int size() {
		return factorial<N>();
	}
	constexpr auto const& operator[](int i) const {
		return elements_[i];
	}
	constexpr auto find(Permutation<N> const &p) const {
		for (int i = 0; i < size(); i++) {
			if (p == elements_[i]) {
				return i;
			}
		}
		return -1;
	}
	friend std::ostream& operator<<(std::ostream &os, SymmetricGroup<N> const &Sn) {
		for (int i = 0; i < size(); i++) {
			os << std::to_string(i + 1) << ". " << Sn[i] << std::endl;
		}
		return os;
	}

private:static consteval auto genElements() {
		std::array<Permutation<N>, size()> elements;
		Permutation<N> p {};
		int i = 0;
		do {
			elements[i++] = p;
		}
		while (std::next_permutation(p.begin(), p.end()));
		return elements;
	}
	static constexpr auto elements_ = genElements();
};

