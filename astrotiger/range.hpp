#pragma once

#include <astrotiger/defs.hpp>
#include <astrotiger/vect.hpp>

#include <vector>
#include <string>
#include <utility>

template<class T>
struct range {
	vect<T> min;
	vect<T> max;

public:
	range() = default;
	range(const vect<T> &pt) {
		min = pt;
		for (int dim = 0; dim < NDIM; dim++) {
			max[dim] = pt[dim] + 1;
		}
	}
	bool empty() const {
		return volume() == 0;
	}
	bool operator==(const range<T> &other) const {
		const bool e = empty();
		if (e == other.empty()) {
			if (!e) {
				for (int dim = 0; dim < NDIM; dim++) {
					if (max[dim] != other.max[dim]) {
						return false;
					}
					if (min[dim] != other.min[dim]) {
						return false;
					}
				}
				return true;
			} else {
				return true;
			}
		} else {
			return false;
		}
	}
	bool operator!=(const range<T> &other) const {
		return !(*this == other);
	}
	range<T> half() const {
		range<T> rc;
		for (int dim = 0; dim < NDIM; dim++) {
			rc.min[dim] = (min[dim] - std::abs(min[dim] % 2)) / 2;
			rc.max[dim] = (max[dim] + std::abs(max[dim] % 2)) / 2;
		}
		return rc;
	}
	std::string to_string() const {
		std::string str;
		for (int dim = 0; dim < NDIM; dim++) {
			str += "(";
			str += std::to_string(min[dim]);
			str += ",";
			str += std::to_string(max[dim]);
			str += ")";
		}
		return str;
	}
	range<T> union_(const range<T> &other) const {
		range<T> u;
		if (!other.empty()) {
			if (!empty()) {
				for (int dim = 0; dim < NDIM; dim++) {
					u.max[dim] = std::max(max[dim], other.max[dim]);
					u.min[dim] = std::min(min[dim], other.min[dim]);
				}
			} else {
				u = other;
			}
		} else {
			u = *this;
		}
		return u;
	}
	std::vector<range<T>> subtract(const range<T> &sub) const {
		std::vector<range<T>> ranges;
		ranges.reserve(2 * NDIM);
		range<T> mid = *this;
//		printf( "Subtracting (%i,%i),(%i,%i) from (%i,%i),(%i,%i)\n", sub.min[0],sub.max[0],sub.min[1],sub.max[1],(*this).min[0],(*this).max[0],(*this).min[1],(*this).max[1]);
		for (int dim = 0; dim < NDIM; dim++) {
			auto lo = mid;
			auto hi = mid;
			lo.max[dim] = std::min(sub.min[dim], lo.max[dim]);
			hi.min[dim] = std::max(sub.max[dim], hi.min[dim]);
			mid.max[dim] = std::min(hi.min[dim], mid.max[dim]);
			mid.min[dim] = std::max(lo.max[dim], mid.min[dim]);
			if (!lo.empty()) {
				ranges.push_back(lo);
			}
			if (!hi.empty()) {
				ranges.push_back(hi);
			}
			if (mid.empty()) {
				break;
			}
		}
//		for( int i = 0; i < ranges.size(); i++) {
//			printf( "++++ %i %i %i %i\n", ranges[i].min[0],ranges[i].max[0],ranges[i].min[1],ranges[i].max[1]);
//		}
		return ranges;
	}
	range<T> double_() const {
		range<T> rc;
		for (int dim = 0; dim < NDIM; dim++) {
			rc.min[dim] = 2 * min[dim];
			rc.max[dim] = 2 * max[dim];
		}
		return rc;
	}
	void set_null() {
		for (int dim = 0; dim < NDIM; dim++) {
			min[dim] = max[dim] = 0;
		}
	}
	bool contains(const range<T> &other) const {
		for (int dim = 0; dim < NDIM; dim++) {
			if (other.min[dim] < min[dim] || other.max[dim] > max[dim]) {
				return false;
			}
		}
		return true;
	}
	std::pair<range<T>, range<T>> split() const {
		int max_dim;
		int max_span = 0;
		for (int dim = 0; dim < NDIM; dim++) {
			const auto span = max[dim] - min[dim];
			if (span > max_span) {
				max_span = span;
				max_dim = dim;
			}
		}
		std::pair<range<T>, range<T>> rc;
		rc.first = rc.second = *this;
		const auto mid = (min[max_dim] + max[max_dim]) / 2;
		rc.first.max[max_dim] = rc.second.min[max_dim] = mid;
		return rc;
	}
	range<T> intersection(const range<T> &other) const {
		range<T> i;
		for (int dim = 0; dim < NDIM; dim++) {
			i.max[dim] = std::min(max[dim], other.max[dim]);
			i.min[dim] = std::max(min[dim], other.min[dim]);
		}
		return i;
	}
	range<T> shift(vect<T> d) const {
		range<T> s;
		for (int dim = 0; dim < NDIM; dim++) {
			s.min[dim] = min[dim] + d[dim];
			s.max[dim] = max[dim] + d[dim];
		}
		return s;
	}
	range pad(int i) const {
		range<T> padded;
		for (int dim = 0; dim < NDIM; dim++) {
			padded.min[dim] = min[dim] - i;
			padded.max[dim] = max[dim] + i;
		}
		return padded;
	}
	vect<T> dims() const {
		vect<T> d;
		for (int dim = 0; dim < NDIM; dim++) {
			d[dim] = max[dim] - min[dim];
		}
		return d;
	}
	T volume() const {
		T vol = T(1);
		for (int dim = 0; dim < NDIM; dim++) {
			vol *= std::max(max[dim] - min[dim], T(0));
		}
		return vol;
	}
	bool contains(const vect<T> &pt) const {
		for (int dim = 0; dim < NDIM; dim++) {
			if (pt[dim] < min[dim] || pt[dim] >= max[dim]) {
				return false;
			}
		}
		return true;
	}
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & min;
		arc & max;
	}
};

