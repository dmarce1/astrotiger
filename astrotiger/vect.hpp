/*
 * vect.hpp
 *
 *      Author: dmarce1
 */

#ifndef VECT_HPP_
#define VECT_HPP_

#include <astrotiger/defs.hpp>
#include <cmath>

template<class T>
class vect {
	T v[NDIM];
public:
	vect() {
	}
	vect(T a);
	T& operator[](int i);
	T operator[](int i) const;
	vect operator-() const;
	vect operator-(const vect &other) const;
	vect operator+(const vect &other) const;
	vect operator*(T r) const;
	vect operator/(T r) const;
	vect operator-=(const vect &other);
	vect operator+=(const vect &other);
	vect operator*=(T r);
	vect operator/=(T r);
	bool operator<(const vect &other) const;
	bool operator<=(const vect &other) const;
	bool operator>(const vect &other) const;
	bool operator>=(const vect &other) const;
	bool operator==(const vect &other) const;
	bool operator!=(const vect &other) const;
	T dot(const vect &other) const;
	template<class U>
	operator vect<U>() const;
	template<class Arc>
	void serialize(Arc &&a, unsigned) {
		for (int i = 0; i < NDIM; i++) {
			a & v[i];
		}
	}
	inline T sum() const {
		T sum = v[0];
		for( int i = 1; i < NDIM; i++) {
			sum += v[i];
		}
		return sum;
	}

};

template<class T>
template<class U>
inline vect<T>::operator vect<U>() const {
	vect<U> a;
	for (int i = 0; i < NDIM; i++) {
		a[i] = (U) v[i];
	}
	return a;
}

template<class T>
bool inline vect<T>::operator<(const vect &other) const {
	for (int n = 0; n < NDIM; n++) {
		if ((*this)[n] < other[n]) {
			return true;
		} else if ((*this)[n] > other[n]) {
			return false;
		}
	}
	return false;
}

template<class T>
bool inline vect<T>::operator<=(const vect &other) const {
	return *this < other || *this == other;
}

template<class T>
bool inline vect<T>::operator>(const vect &other) const {
	return !(*this <= other);
}

template<class T>
bool inline vect<T>::operator>=(const vect &other) const {
	return !(*this < other);
}

template<class T>
inline vect<T>::vect(T a) {

	for (int i = 0; i < NDIM; i++) {
		v[i] = a;
	}
}

template<class T>
inline bool vect<T>::operator==(const vect<T> &other) const {

	for (int dim = 0; dim < NDIM; dim++) {
		if ((*this)[dim] != other[dim]) {
			return false;
		}
	}
	return true;
}

template<class T>
inline bool vect<T>::operator!=(const vect<T> &other) const {
	return !((*this) == other);
}

template<class T>
inline T& vect<T>::operator[](int i) {
	return v[i];
}

template<class T>
inline T vect<T>::operator[](int i) const {
	return v[i];
}

template<class T>
inline vect<T> vect<T>::operator-() const {
	vect<T> result;

	for (int dim = 0; dim < NDIM; dim++) {
		result[dim] = -v[dim];
	}
	return result;
}

template<class T>
inline vect<T> vect<T>::operator-(const vect<T> &other) const {
	vect<T> result;

	for (int dim = 0; dim < NDIM; dim++) {
		result[dim] = v[dim] - other[dim];
	}
	return result;
}

template<class T>
inline vect<T> vect<T>::operator-=(const vect<T> &other) {

	for (int dim = 0; dim < NDIM; dim++) {
		v[dim] -= other[dim];
	}
	return *this;
}

template<class T>
inline vect<T> vect<T>::operator+=(const vect<T> &other) {

	for (int dim = 0; dim < NDIM; dim++) {
		v[dim] += other[dim];
	}
	return *this;
}

template<class T>
inline vect<T> vect<T>::operator*=(T r) {

	for (int dim = 0; dim < NDIM; dim++) {
		v[dim] *= r;
	}
	return *this;
}

template<class T>
inline vect<T> vect<T>::operator/=(T r) {

	for (int dim = 0; dim < NDIM; dim++) {
		v[dim] /= r;
	}
	return *this;
}

template<class T>
inline vect<T> vect<T>::operator+(const vect<T> &other) const {
	vect<T> result;

	for (int dim = 0; dim < NDIM; dim++) {
		result[dim] = v[dim] + other[dim];
	}
	return result;
}

template<class T>
inline vect<T> vect<T>::operator*(T r) const {
	vect<T> result;

	for (int dim = 0; dim < NDIM; dim++) {
		result[dim] = v[dim] * r;
	}
	return result;
}

template<class T>
inline vect<T> vect<T>::operator/(T r) const {
	vect<T> result;

	for (int dim = 0; dim < NDIM; dim++) {
		result[dim] = v[dim] / r;
	}
	return result;
}

template<class T>
inline T vect<T>::dot(const vect<T> &other) const {
	T result = v[0] * other[0];

	for (int dim = 1; dim < NDIM; dim++) {
		result += v[dim] * other[dim];
	}
	return result;
}

template<class T>
inline T abs(const vect<T> &v) {
	return sqrt(v.dot(v));
}

template<class T>
inline vect<T> abs(const vect<T> &a, const vect<T> &b) {
	vect<T> c;

	for (int i = 0; i < NDIM; i++) {
		c[i] = abs(a[i] - b[i]);
	}
	return c;
}

template<class T>
inline vect<T> max(const vect<T> &a, const vect<T> &b) {
	vect<T> c;

	for (int i = 0; i < NDIM; i++) {
		c[i] = max(a[i], b[i]);
	}
	return c;
}

template<class T>
inline vect<T> min(const vect<T> &a, const vect<T> &b) {
	vect<T> c;

	for (int i = 0; i < NDIM; i++) {
		c[i] = min(a[i], b[i]);
	}
	return c;
}

#endif /* VECT_HPP_ */
