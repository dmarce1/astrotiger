///******************************************************************************
// Copyright (C) 2024  Dominic C. Marcello
// *******************************************************************************/

#pragma once

#include <array>
#include <numeric>
#include <vector>

template<typename >
struct Array;

struct Dimensions;

struct Dimensions: public std::vector<size_t> {
	using base_type = std::vector<size_t>;
	Dimensions(auto ...args) :
			base_type(std::initializer_list < size_t > ( {size_t(args)...})) {
	}
};

template<typename T>
struct Array {
	Array() :
			ptr_(nullptr), rank_(0), size_(0) {
	}
	Array(Array const &other) {
		dims_ = other.dims_;
		size_ = other.size_;
		ptr_ = (T*) malloc(sizeof(T) * size_);
		memcpy(ptr_, other.ptr_, sizeof(T) * size_);
	}
	Array(Dimensions const &dims) :
			dims_(dims), rank_(dims.size()) {
		constexpr std::multiplies<size_t> mul;
		size_ = std::accumulate(dims_.begin(), dims_.end(), 1, mul);
		ptr_ = (T*) malloc(sizeof(T) * size_);
	}
	~Array() {
		if (ptr_) {
			free(ptr_);
		}
	}
	struct pointers {
		bool operator++() {
			size_t d = indices_.size() - 1;
			while (++indices_[d] == sizes_[d]) {
				indices_[d] -= sizes_[d];
				vPtr_ -= sizes_[d] * strides_[d];
				if (d == 0) {
					return false;
				}
				d--;
			}
			vPtr_ += strides_[d];
			aPtr_++;
			return true;
		}
		T* viewPointer() {
			return vPtr_;
		}
		T* arrayPointer() {
			return aPtr_;
		}
		size_t chunkSize() const {
			return chunkSize_;
		}
		pointers(std::vector<size_t> const &sizes, std::vector<size_t> const &strides, T *aptr, T *vptr) :
				sizes_(sizes), strides_(strides), vPtr_(vptr), aPtr_(aptr), indices_(sizes.size(), 0) {
			chunkSize_ = 1;
			int d = sizes_.size() - 1;
			while (strides_[d] == chunkSize_) {
				indices_.pop_back();
				chunkSize_ *= sizes_[d];
				if (d == 0) {
					break;
				}
				d--;
			}
		}
	private:
		std::vector<size_t> const &sizes_;
		std::vector<size_t> const &strides_;
		std::vector<size_t> indices_;
		size_t chunkSize_;
		T *vPtr_;
		T *aPtr_;
	};
	struct reference {
		reference(T *ptr, std::vector<size_t> dims) {
			dims_ = std::move(dims);
			ptr_ = ptr;
		}
		reference operator[](size_t index) {
			ptr_ += index * strides_[idx_];
			return *this;
		}
		reference operator()(size_t begin, size_t end) {
			ptr_ += begin * strides_[idx_];
			sizes_[idx_] = end - begin;
			return *this;
		}
		reference operator()(size_t begin, size_t end, size_t stride) {
			ptr_ += begin * strides_[idx_];
			strides_[idx_] *= stride;
			sizes_[idx_] = (end - begin) / stride;
			return *this;
		}
		reference operator=(Array<T> const &other) const {
			pointers ptrs(sizes_, strides_, other.ptr_, ptr_);
			size_t const chunkSize = ptrs.chunkSize();
			if (chunkSize == 1) {
				do {
					*ptrs.viewPointer() = *ptrs.arrayPointer();
				} while (++ptrs);
			} else {
				do {
					memcpy(ptrs.viewPointer(), ptrs.arrayPointer(), chunkSize * sizeof(T));
				} while (++ptrs);
			}
			return *this;
		}
		operator Array<T>() const {
			Array<T> arr(sizes_);
			pointers ptrs(sizes_, strides_, arr.ptr_, ptr_);
			size_t const chunkSize = ptrs.chunkSize();
			if (chunkSize == 1) {
				do {
					*ptrs.arrayPointer() = *ptrs.viewPointer();
				} while (++ptrs);
			} else {
				do {
					memcpy(ptrs.arrayPointer(), ptrs.viewPointer(), chunkSize * sizeof(T));
				} while (++ptrs);
			}
			return arr;
		}
		reference operator=(T const &value) const {
			*ptr_ = value;
		}
		operator T() const {
			return *ptr_;
		}
	private:
		std::vector<size_t> dims_;
		std::vector<size_t> sizes_;
		std::vector<size_t> strides_;
		T *ptr_;
	};
	reference operator[](size_t i) {
		reference ref(ptr_, dims_);
		return ref[i];
	}
	reference operator()(size_t b, size_t e) {
		reference ref(ptr_, dims_);
		return ref(b, e);
	}
	size_t size() const {
		return size_;
	}
private:
	std::vector<size_t> dims_;
	size_t rank_;
	size_t size_;
	T *ptr_;
};

