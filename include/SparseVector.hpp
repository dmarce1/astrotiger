/*
 * SparseVector.hpp
 *
 *  Created on: Jan 22, 2026
 *      Author: dmarce1
 */

#ifndef INCLUDE_SPARSEVECTOR_HPP_
#define INCLUDE_SPARSEVECTOR_HPP_

#include <algorithm>
#include <utility>
#include <vector>

template<typename T>
struct SparseVector;

template<typename T>
struct SparseMatrix;

template<typename T>
struct SparseVector {
	friend class SparseMatrix<T> ;
	constexpr SparseVector() = default;
	constexpr SparseVector(int n) :
			N_(n) {
	}
	constexpr SparseVector(SparseVector const&) = default;
	constexpr SparseVector(SparseVector&&) = default;
	constexpr SparseVector& operator=(SparseVector const&) = default;
	constexpr SparseVector& operator=(SparseVector&&) = default;
	struct reference {
		constexpr reference(SparseVector<T> &ref, int index) :
				ref_(ref), index_(index) {
		}
		constexpr reference& operator=(T value) {
			if (value != zero) {
				ref_.set(index_, value);
			} else {
				ref_.remove(index_);
			}
			return *this;
		}
		constexpr reference& operator*=(T value) {
			if (value != one) {
				if (value == zero) {
					ref_.remove(index_);
				} else {
					auto *ptr = ref_.getPointer(index_);
					if (ptr != nullptr) {
						*ptr *= value;
					}
				}
			}
			return *this;
		}
		constexpr reference& operator/=(T value) {
			*this *= one / value;
			return *this;
		}
		constexpr reference& operator+=(T value) {
			if (value != zero) {
				auto *ptr = ref_.getPointer(index_);
				if (ptr == nullptr) {
					ref_.set(index_, value);
				} else {
					auto const sum = *ptr + value;
					if (sum == zero) {
						ref_.remove(index_);
					} else {
						ref_.set(index_, sum);
					}
				}
			}
			return *this;
		}
		constexpr reference& operator-=(T value) {
			*this += -value;
			return *this;
		}
		constexpr operator T() const {
			return ref_.get(index_).second;
		}
	private:
		SparseVector<T> &ref_ { };
		int index_;
	};
	constexpr SparseVector operator+() const {
		return *this;
	}
	constexpr SparseVector operator-() const {
		return -one * *this;
	}
	constexpr SparseVector operator+(SparseVector const &other) const {
		auto result = *this;
		result += other;
		return result;
	}
	constexpr SparseVector operator-(SparseVector const &other) const {
		auto result = *this;
		result -= other;
		return result;
	}
	constexpr SparseVector operator*(T scale) const {
		auto result = *this;
		result *= scale;
		return result;
	}
	constexpr T operator*(SparseVector const &other) const {
		T sum = zero;
		for (auto i = V_.begin(); i != V_.end(); i++) {
			auto const *ptr = other.getPointer(i->first);
			if (ptr == nullptr) {
				continue;
			}
			sum += i->second * *ptr;
		}
		return sum;
	}
	constexpr SparseVector operator/(T scale) const {
		return *this * (one / scale);
	}
	constexpr SparseVector& operator+=(SparseVector const &other) {
		for (auto const &i : other.V_) {
			auto [exists, value] = get(i.first);
			value += i.second;
			if (value == zero) {
				if (exists) {
					remove(i.first);
				}
			} else {
				set(i.first, value);
			}
		}
		return *this;
	}
	constexpr SparseVector& operator-=(SparseVector const &other) {
		*this += -other;
		return *this;
	}
	constexpr SparseVector& operator*=(T scale) {
		if (scale == zero) {
			V_.clear();
			return *this;
		}
		if (scale == one) {
			return *this;
		}
		for (auto &i : V_) {
			i.second *= scale;
		}
		return *this;
	}
	constexpr SparseVector& operator/=(T scale) {
		*this *= (one / scale);
		return *this;
	}
	constexpr T operator[](int i) const {
		return get(i).second;
	}
	constexpr reference operator[](int i) {
		return reference(*this, i);
	}
	constexpr int density() const {
		return V_.size();
	}
	constexpr int sparsity() const {
		return N_ - density();
	}
	constexpr auto begin() const {
		return V_.begin();
	}
	constexpr auto end() const {
		return V_.end();
	}
	friend constexpr SparseVector operator*(T scale, SparseVector const &v) {
		return v * scale;
	}
	static constexpr SparseVector unit(int n, int N) {
		SparseVector U(N);
		U[n] = one;
		return U;
	}
private:
	NUMERICAL_CONSTANTS (T);
	constexpr static auto cmp = [](auto const &a, auto const &b) {
		return a.first < b.first;
	};
	constexpr auto get(int i) const {
		auto it = std::lower_bound(V_.begin(), V_.end(), std::pair(i, zero), cmp);
		auto rc = it != V_.end();
		if (rc) {
			rc = it->first == i;
		}
		return std::pair(rc, rc ? it->second : zero);
	}
	constexpr auto* getPointer(int i) {
		auto it = std::lower_bound(V_.begin(), V_.end(), std::pair(i, zero), cmp);
		auto rc = it != V_.end();
		if (rc) {
			rc = it->first == i;
		}
		return rc ? &it->second : nullptr;
	}
	constexpr auto const* getPointer(int i) const {
		auto it = std::lower_bound(V_.begin(), V_.end(), std::pair(i, zero), cmp);
		auto rc = it != V_.end();
		if (rc) {
			rc = it->first == i;
		}
		return rc ? &it->second : nullptr;
	}
	constexpr void set(int i, T v) {
		if (v != zero) {
			auto it = std::lower_bound(V_.begin(), V_.end(), std::pair(i, zero), cmp);
			auto rc = it != V_.end();
			if (rc) {
				rc = it->first == i;
			}
			if (rc) {
				it->second = v;
			} else {
				V_.insert(it, { i, v });
			}
		} else {
			remove(i);
		}
	}
	constexpr bool remove(int i) {
		auto it = std::lower_bound(V_.begin(), V_.end(), std::pair(i, zero), cmp);
		if (it == V_.end()) {
			return false;
		}
		if (it->first == i) {
			V_.erase(it);
			return true;
		}
		return false;
	}
	std::vector<std::pair<int, T>> V_ { };
	int N_ { };
};

#endif /* INCLUDE_SPARSEVECTOR_HPP_ */
