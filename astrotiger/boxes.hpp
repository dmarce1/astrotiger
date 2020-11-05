#pragma once

#include <astrotiger/multi_array.hpp>
#include <memory>

std::shared_ptr<std::vector<int>> get_red_indices(const multi_range &box);
std::shared_ptr<std::vector<int>> get_black_indices(const multi_range &box);
