/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "Definitions.hpp"

#include <array>

template<Dimension D, Index N>
struct Multidices {
	constexpr Multidices() = default;
	constexpr Multidices(Multidices const&) = default;
	constexpr Multidices(Multidices&&) = default;
	constexpr Multidices& operator=(Multidices const&) = default;
	constexpr Multidices& operator=(Multidices&&) = default;
	constexpr Index operator[](Dimension d) const {
		return idx_[d];
	}
	constexpr Index& operator[](Dimension d) {
		return idx_[d];
	}
	constexpr explicit operator Index() const {
		Index i = 0;
		for (Dimension d = 0; d < D; d++) {
			i = N * i + idx_[d];
		}
		return i;
	}
	constexpr Multidices& operator++() {
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
	constexpr Multidices operator++(int) {
		auto const original = *this;
		operator++();
		return original;
	}
	constexpr Multidices& operator--() {
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
	constexpr Multidices operator--(int) {
		auto const original = *this;
		operator++();
		return original;
	}
	constexpr bool operator<(Multidices const &other) const {
		for (Dimension d = 0; d < D; d++) {
			if (idx_[d] < other.idx_[d]) {
				return true;
			}
		}
		return false;
	}
	constexpr bool operator<=(Multidices const &other) const {
		return !(*this > other);
	}
	constexpr bool operator==(Multidices const &other) const {
		for (Dimension d = 0; d < D; d++) {
			if (idx_[d] != other.idx_[d]) {
				return false;
			}
		}
		return true;
	}
	constexpr bool operator>=(Multidices const &other) const {
		return !(*this < other);
	}
	constexpr bool operator>(Multidices const &other) const {
		return other < *this;
	}
	constexpr bool operator!=(Multidices const &other) const {
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
	constexpr static Multidices ibegin() {
		std::array<Index, D> idx;
		idx.fill(0);
		return Multidices(idx);
	}
	constexpr static Multidices iend() {
		std::array<Index, D> idx;
		idx.fill(N);
		return Multidices(idx);
	}
	constexpr static Multidices irbegin() {
		std::array<Index, D> idx;
		idx.fill(N - 1);
		return Multidices(idx);
	}
	constexpr static Multidices irend() {
		std::array<Index, D> idx;
		idx.fill(~Index(0));
		return Multidices(idx);
	}
private:
	constexpr Multidices(std::array<Index, D> const &idx) :
			idx_(idx) {
	}
	std::array<Index, D> idx_;
};
