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

template <typename T>
struct SparseVector;

template <typename T>
struct SparseMatrix;

template <typename T>
struct SparseVector {
	friend class SparseMatrix<T>;
	constexpr SparseVector() = default;
	constexpr SparseVector(int n) :
		N_(n) {
	}
	constexpr SparseVector(SparseVector const &) = default;
	constexpr SparseVector(SparseVector &&) = default;
	constexpr SparseVector &operator=(SparseVector const &) = default;
	constexpr SparseVector &operator=(SparseVector &&) = default;
	struct reference {
		constexpr reference(SparseVector<T> &ref, int index) :
			ref_(ref), index_(index) {
		}
		constexpr reference &operator=(T value) {
			if (!isZero(value)) {
				ref_.set(index_, value);
			} else {
				ref_.remove(index_);
			}
			return *this;
		}
		constexpr reference &operator*=(T value) {
			if (value != one) {
				if (isZero(value)) {
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
		constexpr reference &operator/=(T value) {
			*this *= one / value;
			return *this;
		}
		constexpr reference &operator+=(T value) {
			if (!isZero(value)) {
				auto *ptr = ref_.getPointer(index_);
				if (ptr == nullptr) {
					ref_.set(index_, value);
				} else {
					auto const sum = *ptr + value;
					if (isZero(sum)) {
						ref_.remove(index_);
					} else {
						ref_.set(index_, sum);
					}
				}
			}
			return *this;
		}
		constexpr reference &operator-=(T value) {
			*this += -value;
			return *this;
		}
		constexpr operator T() const {
			return ref_.get(index_).second;
		}

	private:
		SparseVector<T> &ref_{};
		int index_;
	};
	template <typename U>
	constexpr operator SparseVector<U>() const {
		SparseVector<U> vec(N_);
		for (auto i = begin(); i != end(); i++) {
			vec[i->first] = static_cast<U>(i->second);
		}
		return vec;
	}
	friend constexpr auto normalize(SparseVector const &v) {
		using std::sqrt;
		auto const v2 = v * v;
		if (isZero(v2)) {
			return v;
		} else {
			return v / sqrt(v2);
		}
	}
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
	constexpr SparseVector &operator+=(SparseVector const &other) {
		for (auto const &i : other.V_) {
			auto [exists, value] = get(i.first);
			value += i.second;
			if (isZero(value)) {
				if (exists) {
					remove(i.first);
				}
			} else {
				set(i.first, value);
			}
		}
		return *this;
	}
	constexpr SparseVector &operator-=(SparseVector const &other) {
		*this += -other;
		return *this;
	}
	constexpr SparseVector &operator*=(T scale) {
		if (isZero(scale)) {
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
	constexpr SparseVector &operator/=(T scale) {
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
	constexpr int size() const {
		return N_;
	}
	constexpr auto begin() {
		return V_.begin();
	}
	constexpr auto end() {
		return V_.end();
	}
	constexpr auto begin() const {
		return V_.cbegin();
	}
	constexpr auto end() const {
		return V_.cend();
	}
	friend constexpr SparseVector operator*(T scale, SparseVector const &v) {
		return v * scale;
	}
	static constexpr SparseVector unit(int n, int N) {
		SparseVector U(N);
		U[n] = one;
		return U;
	}
	friend constexpr std::ostream &operator<<(std::ostream &os, SparseVector<T> const &A) {
		NUMERICAL_CONSTANTS(T);
		auto const N = A.size();
		SparseVector<int> widths(N);
		int maxWidth = 3;
		for (int n = 0; n < N; n++) {
			if (!isZero(A[n])) {
				std::ostringstream ss;
				ss << A[n];
				widths[n] = ss.str().size();
				maxWidth = std::max(maxWidth, (int)ss.str().size());
			}
		}
		for (int m = 0; m < N; m++) {
			os << '|';
			if (!isZero(A[m])) {
				int const padding = maxWidth - widths[m];
				os << std::string(padding / 2, ' ') << A[m] << std::string((padding + 1) / 2, ' ');
			} else {
				os << std::string(maxWidth, ' ');
			}
		}
		os << "|\n";
		return os;
	}

private:
	NUMERICAL_CONSTANTS(T);
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
	constexpr auto *getPointer(int i) {
		auto it = std::lower_bound(V_.begin(), V_.end(), std::pair(i, zero), cmp);
		auto rc = it != V_.end();
		if (rc) {
			rc = it->first == i;
		}
		return rc ? &it->second : nullptr;
	}
	constexpr auto const *getPointer(int i) const {
		auto it = std::lower_bound(V_.begin(), V_.end(), std::pair(i, zero), cmp);
		auto rc = it != V_.end();
		if (rc) {
			rc = it->first == i;
		}
		return rc ? &it->second : nullptr;
	}
	constexpr void set(int i, T v) {
		if (!isZero(v)) {
			auto it = std::lower_bound(V_.begin(), V_.end(), std::pair(i, zero), cmp);
			auto rc = it != V_.end();
			if (rc) {
				rc = it->first == i;
			}
			if (rc) {
				it->second = v;
			} else {
				V_.insert(it, {i, v});
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
	static constexpr bool isZero(T const &value) {
		if constexpr (std::is_floating_point<T>::value) {
			constexpr T ε = std::sqrt(std::numeric_limits<T>::epsilon());
			return ((-ε < value) && (value < +ε));
		} else {
			return value == zero;
		}
	}
	std::vector<std::pair<int, T>> V_{};
	int N_{};
};

template <typename T>
constexpr bool areLinearlyIndependent(SparseVector<T> const &u, SparseVector<T> const &v) {
	NUMERICAL_CONSTANTS(T);
	if ((u.density() == 0) || (v.density() == 0) || (u.density() != v.density())) {
		return true;
	}
	auto i = u.begin();
	auto j = v.begin();
	if (i->first != j->first) {
		return true;
	}
	auto const scale = j->second / i->second;
	while (1) {
		i++;
		j++;
		if (i == u.end()) {
			break;
		}
		if (i->first != j->first) {
			return true;
		}
		if (j->second != scale * i->second) {
			return true;
		}
	};
	return false;
}

#endif /* INCLUDE_SPARSEVECTOR_HPP_ */
