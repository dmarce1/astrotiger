/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
*******************************************************************************/

#pragma once

#include <array>
#include <ostream>

template<typename T, std::size_t N>
std::ostream &operator<<(std::ostream &os, std::array<T, N> const &a) {
    os << '{';
    for (std::size_t i = 0; i < N; i++) {
        if (i != 0) {
            os << ", ";
        }
        os << a[i];
    }
    os << '}';
    return os;
}


