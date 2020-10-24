#pragma once

#include <astrotiger/multi_array.hpp>
#include <astrotiger/geometry.hpp>


class geometry {
	multi_range box;
	std::vector<multi_index> directions;
	std::vector<multi_array<index_type>> multi_to_one;
	std::vector<std::vector<multi_index>> one_to_multi;
public:
	geometry(const multi_range&);
};
