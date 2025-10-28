#pragma once

#define DEFINE_VECTOR_OPERATORS(NAME, TYPE, ...)                            \
	friend NAME operator+(NAME a) {                                         \
		return a;                                                           \
	}                                                                       \
	friend NAME operator-(NAME a) {                                         \
		std::apply([](auto&... members) {                                   \
			((members = -members), ...);                                    \
		}, std::tie(__VA_ARGS__));                                          \
		return a;                                                           \
	}                                                                       \
	friend NAME operator*(TYPE scale, NAME a) {                             \
		return a * scale;                                                   \
	}                                                                       \
	friend NAME operator+(NAME a, NAME const& b) {                          \
		std::apply([&](auto&... lhsMembers) {                               \
			auto const rhs = std::tie(__VA_ARGS__);                         \
			int i = 0;                                                      \
            ((lhsMembers += std::get<i++>(rhs)), ...);                      \
		}, std::tie(__VA_ARGS__));                                          \
		return a;                                                           \
	}                                                                       \
	friend NAME operator-(NAME a, NAME const& b) {                          \
		return a + (-b);                                                    \
	}                                                                       \
	friend NAME operator*(NAME a, TYPE scale) {                             \
		std::apply([&](auto&... members) {                                  \
        	((members *= scale), ...);                                      \
        }, std::tie(__VA_ARGS__));                                          \
		return a;                                                           \
	}                                                                       \
	friend NAME operator/(NAME a, TYPE invScale) {                          \
		return a * (TYPE(1.0) / invScale);                                  \
	}                                                                       \
	friend NAME& operator+=(NAME& a, NAME const& b) {                       \
		a = a + b;                                                          \
		return a;                                                           \
	}                                                                       \
	friend NAME& operator-=(NAME& a, NAME const& b) {                       \
		a = a + (-b);                                                       \
		return a;                                                           \
	}                                                                       \
	friend NAME& operator*=(NAME& a, TYPE b) {                              \
		a = a * b;                                                          \
		return a;                                                           \
	}                                                                       \
	friend NAME& operator/=(NAME& a, TYPE b) {                              \
		a = a / b;                                                          \
		return a;                                                           \
	}                                                                       \

