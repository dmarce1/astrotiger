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
#include <cmath>

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
	std::vector<T> a;
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
	inline index_type index(const multi_index &I) const {
		index_type i = I[0] - box.min[0];
		for (int dim = 1; dim < NDIM; dim++) {
			i = i * dims[dim] + (I[dim] - box.min[dim]);
		}
		assert(i >= 0 && i < a.size());
		return i;
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
		a.resize(size);
	}
	inline multi_index get_strides() const {
		return stride;
	}
	inline T operator[](const multi_index &I) const {
		return a[index(I)];
	}
	inline T& operator[](const multi_index &I) {
		return a[index(I)];
	}
	T smooth_gradient(int dim, const multi_index &I) const {
		const auto i = index(I);
		const auto ip = i + stride[dim];
		const auto im = i - stride[dim];
		assert(ip < a.size());
		assert(ip >= 0);
		assert(im < a.size());
		assert(im >= 0);
		assert(i < a.size());
		assert(i >= 0);
		return 0.5 * (a[ip] - a[im]);
	}
	T minmod_gradient(int dim, const multi_index &I) const {
		const auto i = index(I);
		const auto ip = i + stride[dim];
		const auto im = i - stride[dim];
		assert(ip < a.size());
		assert(ip >= 0);
		assert(im < a.size());
		assert(im >= 0);
		assert(i < a.size());
		assert(i >= 0);
		const auto c = a[ip] - a[i];
		const auto b = a[i] - a[im];
		return (std::copysign(0.5, c) + std::copysign(0.5, b)) * std::min(std::abs(c), std::abs(b));
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
	const T* data() const {
		return a.data();
	}
	T* data() {
		return a.data();
	}
	multi_array prolong(const multi_range &pro_box, bool smooth = false) const {
		multi_range grad_box;
		for (int dim = 0; dim < NDIM; dim++) {
			grad_box.min[dim] = (pro_box.min[dim] - (std::abs(pro_box.min[dim]) % 2)) / 2;
			grad_box.max[dim] = (pro_box.max[dim] + (std::abs(pro_box.max[dim]) % 2)) / 2;
		}
		multi_array grad[NDIM];
		for (int dim = 0; dim < NDIM; dim++) {
			grad[dim].resize(grad_box);
			for (multi_iterator i(grad_box); !i.end(); i++) {
				if (smooth) {
					grad[dim][i] = smooth_gradient(dim, i);
				} else {
					grad[dim][i] = minmod_gradient(dim, i);
				}
			}
		}
		multi_array pro(pro_box);
		for (multi_iterator i(pro_box); !i.end(); i++) {
			vect<double> c0;
			for (int dim = 0; dim < NDIM; dim++) {
				c0[dim] = 0.25 * (2 * (std::abs(i[dim] % 2)) - 1);
			}
			multi_index j = i.index();
			for (int dim = 0; dim < NDIM; dim++) {
				j[dim] = (j[dim] - std::abs(j[dim] % 2)) / 2;
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
		arc & a;
	}
};


struct boundary {
	std::vector<multi_range> boxes;
	std::vector<std::vector<double>> data;

	template<class A>
	void serialize(A&& arc, unsigned) {
		arc & boxes;
		arc & data;
	}
};



#endif /* ASTROTIGER_MULTI_ARRAY_HPP_ */
