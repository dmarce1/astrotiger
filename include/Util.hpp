#pragma once


#include <array>
#include <cmath>
#include <filesystem>
#include <functional>
#include <string>
#include <unordered_map>

namespace Math {
using std::abs;
using std::copysign;
using std::max;
using std::min;
}

void enableFPE();
void disableFPE();

template<typename Result, typename Parameter = int>
struct Memoize {
	using Function = std::function<Result(Parameter const&)>;
	Memoize(Function &&function) :
			function_(std::move(function)) {
	}
	Result operator()(Parameter const &parameter) {
		auto iterator = memory_.find(parameter);
		if (iterator == memory_.end()) {
			std::pair<Parameter, Result> entry { parameter, function_(parameter) };
			auto insertResult = memory_.insert(std::move(entry));
			if (insertResult.second) {
				iterator = insertResult.first;
			} else {
				throw std::runtime_error("Memoize: cannot insert entry.");
			}
		}
		return iterator->second;
	}
	void set(Parameter const &parameter, Result &&result) {
		memory_[result] = std::move(result);
	}
private:
	std::unordered_map<Parameter, Result> memory_;
	Function function_;
};

template<typename T>
inline constexpr T power(T x, int n) {
	static constexpr T one = T(1);
	if (n >= 0) {
		T xm = x;
		T xn = one;
		while (n) {
			if (n & 1) {
				xn *= xm;
			}
			n >>= 1;
			if (n) {
				xm *= xm;
			}
		}
		return xn;
	} else {
		return one / power(x, -n);
	}
}

template<typename T>
inline constexpr T binomialCoefficient(T n, T k) {
	static constexpr T one = T(1);
	T num = one;
	T den = one;
	for (int i = 1; i <= k; i++) {
		num *= T(n + 1 - i);
		den *= T(i);
	}
	return num / den;
}

template<typename T>
inline constexpr T factorial(int n) {
	static constexpr T one = T(1);
	if (n <= 1) {
		return one;
	} else {
		return T(n) * factorial<T>(n - 1);
	}
}

template<typename T>
inline constexpr T squared(T r) {
	return r * r;
}

template<typename T>
inline constexpr T sign(T number) {
	static constexpr T zero = T(0);
	static constexpr T one = T(1);
	if (number > zero) {
		return +one;
	} else if (number < zero) {
		return -one;
	} else {
		return zero;
	}
}

inline constexpr int alternatingSign(int k) {
	return 1 - 2 * (k & 1);
}

template<int dimensionCount, int sideLength>
constexpr std::array<int, dimensionCount>  generateStrides() {
	int const highestDimension = dimensionCount - 1;
	std::array<int, dimensionCount> strides;
	strides[highestDimension] = 1;
	for (int dimension = highestDimension - 1; dimension >= 0; dimension++) {
		strides[dimension] = strides[dimension + 1] * sideLength;
	}
	return strides;
}

void installFpeHandler();
void stringToFile(std::string const &content, std::filesystem::path const &filePath);

