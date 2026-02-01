/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "Definitions.hpp"

#include <array>

template<Dimension D, Index N>
struct Indices {
	constexpr Indices() = default;
	constexpr Indices(Indices const&) = default;
	constexpr Indices(Indices&&) = default;
	constexpr Indices& operator=(Indices const&) = default;
	constexpr Indices& operator=(Indices&&) = default;
	constexpr Index operator[](Dimension d) const {
		return idx_[d];
	}
	constexpr Index& operator[](Dimension d) {
		return idx_[d];
	}
	constexpr operator Index() const {
		Index i = 0;
		for (Dimension d = 0; d < D; d++) {
			i = N * i + idx_[d];
		}
		return i;
	}
	constexpr Indices& operator++() {
		Dimension d = D - 1;
		while (++idx_[d] == N) {
			idx_[d] = 0;
			if (d == 0) {
				*this = iend();
				break;
			}
			d--;
		}
		return *this;
	}
	constexpr Indices operator++(int) {
		auto const original = *this;
		operator++();
		return original;
	}
	constexpr Indices& operator--() {
		Dimension d = D - 1;
		while (idx_[d]-- == 0) {
			idx_[d] = N - 1;
			if (d == 0) {
				*this = irend();
				break;
			}
			d--;
		}
		return *this;
	}
	constexpr Indices operator--(int) {
		auto const original = *this;
		operator++();
		return original;
	}
	constexpr bool operator<(Indices const &other) const {
		for (Dimension d = 0; d < D; d++) {
			if (idx_[d] < other.idx_[d]) {
				return true;
			}
		}
		return false;
	}
	constexpr bool operator<=(Indices const &other) const {
		return !(*this > other);
	}
	constexpr bool operator==(Indices const &other) const {
		for (Dimension d = 0; d < D; d++) {
			if (idx_[d] != other.idx_[d]) {
				return false;
			}
		}
		return true;
	}
	constexpr bool operator>=(Indices const &other) const {
		return !(*this < other);
	}
	constexpr bool operator>(Indices const &other) const {
		return other < *this;
	}
	constexpr bool operator!=(Indices const &other) const {
		return !(*this == other);
	}
	constexpr auto begin() const {
		return idx_.cbegin();
	}
	constexpr auto end() const {
		return idx_.cend();
	}
	constexpr auto begin() {
		return idx_.begin();
	}
	constexpr auto end() {
		return idx_.end();
	}
	template<int I, int J>
	constexpr Indices<D - 2, N> contract() const {
		Indices<D - 2, N> rc;
		for (int i = 0, j = 0; i < D; i++) {
			if ((i != I) && (i != J)) {
				rc[j++] = (*this)[i];
			}
		}
		return rc;
	}
	constexpr static Indices ibegin() {
		std::array<Index, D> idx;
		idx.fill(0);
		return Indices(idx);
	}
	constexpr static Indices iend() {
		std::array<Index, D> idx;
		idx.fill(N);
		return Indices(idx);
	}
	constexpr static Indices irbegin() {
		std::array<Index, D> idx;
		idx.fill(N - 1);
		return Indices(idx);
	}
	constexpr static Indices irend() {
		std::array<Index, D> idx;
		idx.fill(~Index(0));
		return Indices(idx);
	}
	constexpr Indices(std::array<Index, D> const &idx) :
			idx_(idx) {
	}
private:
	std::array<Index, D> idx_;
};
