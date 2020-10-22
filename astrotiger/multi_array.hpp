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
	multi_index index_;
public:
	multi_iterator(const multi_range &box_) :
			box(box_) {
		index_ = box.min;
	}
	operator multi_index() const {
		return index_;
	}
	bool end() const {
		return index_ == box.max;
	}
	index_type operator[](int i) const {
		return index_[i];
	}
	multi_index index() const {
		return index_;
	}
	multi_iterator& operator++(int) {
		int dim = NDIM - 1;
		while (++index_[dim] == box.max[dim]) {
			index_[dim] = box.min[dim];
			dim--;
			if (dim == -1) {
				index_ = box.max;
				break;
			}
		}
		return *this;
	}
};

template<class T>
class multi_array {
	multi_index dims;
	multi_index stride;
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
	multi_range box;
	inline multi_array() {
		dims = 0.0;
		size = 0;
		stride = 0;
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
		stride[NDIM - 1] = 1;
		for (int dim = NDIM - 1; dim > 0; dim--) {
			stride[dim - 1] = dims[dim] * stride[dim];
		}
		data.resize(size);
	}
	inline T operator[](const multi_index &I) const {
		return data[index(I)];
	}
	inline T& operator[](const multi_index &I) {
		return data[index(I)];
	}
	T smooth_gradient(int dim, const multi_index &I) const {
		const auto i = index(I);
		const auto ip = i + stride[dim];
		const auto im = i - stride[dim];
		assert(ip < data.size());
		assert(ip >= 0);
		assert(im < data.size());
		assert(im >= 0);
		assert(i < data.size());
		assert(i >= 0);
		return data[ip] - data[im];
	}
	T minmod_gradient(int dim, const multi_index &I) const {
		const auto i = index(I);
		const auto ip = i + stride[dim];
		const auto im = i - stride[dim];
		assert(ip < data.size());
		assert(ip >= 0);
		assert(im < data.size());
		assert(im >= 0);
		assert(i < data.size());
		assert(i >= 0);
		const auto a = data[ip] - data[i];
		const auto b = data[i] - data[im];
		return (std::copysign(0.5, a) + std::copysign(0.5, b)) * std::min(std::abs(a), std::abs(b));
	}
	multi_array restrict_(const multi_range &res_box) const {
		multi_array res(res_box);
		for (multi_iterator i(res_box); !i.end(); i++) {
			const multi_range cbox = multi_range(i.index()).double_();
			res[i] = 0.0;
			for (multi_iterator j(cbox); !j.end(); j++) {
				res[i] += (*this)[j];
			}
			res[i] /= std::pow(2, NDIM);
		}
		return res;
	}
	multi_array prolong(const multi_range &pro_box) const {
		multi_range grad_box;
		for (int dim = 0; dim < NDIM; dim++) {
			grad_box.min[dim] = (pro_box.min[dim] - (pro_box.min[dim] % 2)) / 2;
			grad_box.max[dim] = (pro_box.max[dim] + (pro_box.max[dim] % 2)) / 2;
		}
		multi_array grad[NDIM];
		for (int dim = 0; dim < NDIM; dim++) {
			grad[dim].resize(grad_box);
			for (multi_iterator i(grad_box); !i.end(); i++) {
				grad[dim][i] = minmod_gradient(dim, i);
			}
		}
		multi_array pro(pro_box);
		for (multi_iterator i(pro_box); !i.end(); i++) {
			vect<double> c0;
			for (int dim = 0; dim < NDIM; dim++) {
				c0[dim] = 0.25 * (2 * (i[dim] % 2) - 1);
			}
			multi_index j = i.index();
			for (int dim = 0; dim < NDIM; dim++) {
				j[dim] = (j[dim] - j[dim] % 2) / 2;
			}
			double sum = (*this)[j];
			for (int dim = 0; dim < NDIM; dim++) {
				sum += c0[dim] * grad[dim][j];
			}
			pro[i] = sum;
		}
		return std::move(pro);
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
