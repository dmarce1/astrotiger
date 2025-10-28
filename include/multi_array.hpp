/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include <array>
#include <numeric>
#include <functional>

#include "interval.hpp"
#include "util.hpp"
#include "vector.hpp"

namespace detail {

template<int J, int K>
constexpr auto swapIndices(auto A) {
	std::swap(A[J], A[K]);
	return A;
}

template<int J, int D>
constexpr auto deleteIndex(Interval<int, D> I) {
	constexpr auto B1 = I.literal().first;
	constexpr auto E1 = I.literal().end;
	std::array<int, D - 1> B2;
	std::array<int, D - 1> E2;
	for(size_t i = 0, j = 0; i < D; i++, j++) {
		if(i == J) {
			i++;
		}
		B2[j] = B1[i];
		E2[j] = E1[i];
	}
	return Interval<int, D - 1>(B2, E2);
}

template<int J, int D>
constexpr auto deleteIndex(std::array<int, D> A) {
	std::array<int, D - 1> B;
	for(size_t i = 0, j = 0; i < A.size(); i++, j++) {
		if(i == J) {
			i++;
		}
		B[j] = A[i];
	}
	return B;
}

}

template<typename T, int D, std::array<int, D> dims>
struct MultiArray { // @formatter:off
	DEFAULT_MEMBERS(MultiArray); // @formatter:on
	static constexpr int size() {
		return box_.size();
	}
	template<int J>
	auto slice(int k) {
		constexpr Interval<int, D> sliceBox(deleteIndex<J>(box_));
		constexpr auto map = [k, sliceBox]() {
			auto begin = box_.literal().first;
			auto end = box_.literal().second;
			std::array<int, sliceBox.size()> map;
			bool exit = false;
			begin[J] = k;
			end[J] = k + 1;
			int sliceIndex = 0;
			std::array<int, D> I = begin;
			while (!exit) {
				map[sliceIndex] = box_.flatten(I);
				int d = D - 1;
				while (++I[d] == end[d]) {
					if (d == 0) {
						exit = true;
						break;
					}
					I[d] = begin[d];
					d--;
				}
				sliceIndex++;
			}
			return map;
		}();
		MultiArray<T, D, sliceBox.span()> sliceBoxArray;
		for (size_t i = 0; i < sliceBox.size(); i++) {
			sliceBoxArray.data_[i] = data_[map[i]];
		}
		return sliceBoxArray;
	}
	template<int J, int K>
	auto transpose() {
		MultiArray const &Ain = *this;
		using namespace detail;
		constexpr auto map = []() {
			constexpr auto begin = box_.literal().first;
			constexpr auto end = box_.literal().second;
			std::array<int, size()> map;
			bool exit = false;
			int i = 0;
			std::array<int, D> I = begin;
			while (!exit) {
				map[i] = box_.flatten(swapIndices<J, K>(I));
				int d = D - 1;
				while (++I[d] == end[d]) {
					if (d == 0) {
						exit = true;
						break;
					}
					I[d] = begin[d];
					d--;
				}
				i++;
			}
			return map;
		}();
		MultiArray<T, D, swapIndices<J, K>(dims)> Aout;
		for (size_t i = 0; i < size(); i++) {
			Aout.data_[i] = Ain.data_[map[i]];
		}
		return Aout;
	}
	template<std::pair<std::array<int, D>, std::array<int, D>> literalBox>
	auto subarray() const {
		constexpr Interval<int, D> subBox(literalBox);
		constexpr auto map = [subBox]() {
			constexpr auto begin = literalBox.first;
			constexpr auto end = literalBox.second;
			std::array<int, subBox.size()> map;
			bool exit = false;
			int subIndex = 0;
			std::array<int, D> I = begin;
			while (!exit) {
				map[subIndex] = box_.flatten(I);
				int d = D - 1;
				while (++I[d] == end[d]) {
					if (d == 0) {
						exit = true;
						break;
					}
					I[d] = begin[d];
					d--;
				}
				subIndex++;
			}
			return map;
		}();
		MultiArray<T, D, subBox.span()> subArray;
		for (size_t i = 0; i < subBox.size(); i++) {
			subArray.data_[i] = data_[map[i]];
		}
		return subArray;
	}
	template<typename, int D2, std::array<int, D2>>
	friend struct MultiArray;
private:
	static constexpr auto box_ = Interval<int, D>(dims);
	std::array<T, size()> data_;
};

