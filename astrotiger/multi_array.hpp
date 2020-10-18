/*
 * multi_array.hpp
 *
 *  Created on: Oct 14, 2020
 *      Author: dmarce1
 */

#ifndef ASTROTIGER_MULTI_ARRAY_HPP_
#define ASTROTIGER_MULTI_ARRAY_HPP_
#include <astrotiger/range.hpp>
#include <astrotiger/vect.hpp>
#include <cassert>
#include <cstdint>
#include <vector>

using index_type = std::int64_t;

using multi_index = vect<index_type>;
using multi_range = range<index_type>;

class multi_iterator {
	const multi_range box;
	multi_index index;
public:
	multi_iterator(const multi_range &box_) :
			box(box_) {
		index = box.min;
	}
	operator multi_index() const {
		return index;
	}
	bool end() const {
		return index == box.max;
	}
	index_type operator[](int i) const {
		return index[i];
	}
	multi_iterator& operator++(int) {
		int dim = NDIM - 1;
		while (++index[dim] == box.max[dim]) {
			index[dim] = box.min[dim];
			dim--;
			if (dim == -1) {
				index = box.max;
				break;
			}
		}
		return *this;
	}
};

template<class T>
class multi_array {
	multi_index dims;
	multi_range box;
	std::size_t size;
	std::vector<T> data;
	inline index_type index(const multi_index &I) const {
		index_type i = I[0] - box.min[0];
		for (int dim = 1; dim < NDIM; dim++) {
			i = i * dims[dim] + (I[dim] - box.min[dim]);
		}
		assert(i >= 0 && i < data.size());
		return i;
	}
public:
	inline multi_array() {
		dims = 0.0;
		size = 0;
	}
	inline multi_array(const multi_range &box_) {
		resize(box_);
	}
	multi_array& operator=(const multi_array&) = default;
	multi_array& operator=(multi_array&&) = default;
	multi_array(const multi_array&) = default;
	multi_array(multi_array&&) = default;
	inline void resize(const multi_range &box_) {
		box = box_;
		dims = box.dims();
		size = dims[0];
		for (int i = 1; i < NDIM; i++) {
			size *= dims[i];
		}
		data.resize(size);
	}
	inline T operator[](const multi_index &I) const {
		return data[index(I)];
	}
	inline T& operator[](const multi_index &I) {
		return data[index(I)];
	}

	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & box;
		arc & dims;
		arc & size;
		arc & data;
	}
};

#endif /* ASTROTIGER_MULTI_ARRAY_HPP_ */
