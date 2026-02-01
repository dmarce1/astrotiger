/******************************************************************************
 Copyright (C) 202#include "IntegerPartition.hpp"
 4  Dominic C. Marcello
 *******************************************************************************/

#include "IntegerPartition.hpp"
#include "TensorSymmetry.hpp"
#include "Young.hpp"
#include <hpx/init.hpp>

#include <fstream>
#include <string>
#include <vector>

struct TensorDefinition {
	int order{};
	int dimension{};
	std::vector<std::vector<int>> symmetries{};
	std::vector<std::vector<int>> antisymmetries{};
	std::vector<std::vector<int>> tracefree{};
	friend std::ostream &operator<<(std::ostream &os, TensorDefinition const &t) {
		auto const printList = [&](auto const &list) {
			for (auto const &v : list) {
				os << "(";
				for (unsigned i = 0; i < v.size(); i++) {
					os << v[i];
					if (i + 1 != v.size()) {
						os << ',';
					}
				}
				os << ") ";
			}
		};
		os << "Order = " << t.order << " ";
		os << "Dimension = " << t.dimension << " ";
		if (t.symmetries.size()) {
			os << "Symmetries = ";
			printList(t.symmetries);
		}
		if (t.antisymmetries.size()) {
			os << "Antisymmetries = ";
			printList(t.antisymmetries);
		}
		if (t.tracefree.size()) {
			os << "Tracefree = ";
			printList(t.tracefree);
		}
		os << std::endl;
		return os;
	}
};

std::vector<TensorDefinition> readTensorDefinitionFile(std::string const &filename) {
	std::vector<TensorDefinition> definitions;
	std::ifstream file(filename);
	if (!file) {
		throw std::runtime_error("Failed to open file: " + filename);
	}
	std::string line;
	while (std::getline(file, line)) {
		auto const readSpaces = [&](std::string::iterator &it) {
			while (isspace(*it) && (it != line.end())) {
				it++;
			}
		};
		auto const readNum = [&](std::string::iterator &it) {
			readSpaces(it);
			int i = 0;
			while (isdigit(*it) && (it != line.end())) {
				i = 10 * i + (*it - '0');
				it++;
			}
			return i;
		};
		auto const readChar = [&](std::string::iterator &it) {
			readSpaces(it);
			return std::tolower(*it);
		};
		auto const readNumList = [&](std::string::iterator &it) {
			std::vector<int> list;
			while (1) {
				list.push_back(readNum(it));
				if (readChar(it) == ',') {
					it++;
				} else {
					break;
				}
			}
			return list;
		};
		int k = 0;
		TensorDefinition tensor{};
		for (auto it = line.begin(); it != line.end(); it++) {
			switch (readChar(it)) {
			case '#':
				it = line.end();
				break;
			case '[':
				tensor.symmetries.push_back(readNumList(++it));
				break;
			case '(':
				tensor.antisymmetries.push_back(readNumList(++it));
				break;
			case '|':
				tensor.tracefree.push_back(readNumList(++it));
				break;
			default:
				if (isdigit(*it) && k < 2) {
					(k++ ? tensor.order : tensor.dimension) = readNum(it);
				}
				break;
			}
		}
		if (tensor.order) {
			definitions.push_back(std::move(tensor));
		}
	}
	return definitions;
}

// x x x x
// x x x
// x x
// x

// c r Λₐ[c] Λₛ[Λₐ[c]]
// 1 0   0      4
// 1 1   1      3
// 1 2   2      2
// 1 3   3      1
// 1 4   4
// auto view = idx | std::views::transform([&](std::size_t i) {
//    return data[i];
//});

#include "Permutation.hpp"
#include <algorithm>
#include <functional>
#include <numeric>
#include <ranges>
#include <utility>

int hpx_main(int argc, char *argv[]) {
	constexpr IntegerPartition<3, 3> Λ;
	constexpr int O = Λ.size();
	constexpr int D = 4;
	constexpr auto ele = symmetrizer<Λ, D>();
	for (auto e : ele) {
		std::cout << e << std::endl;
	}
	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
	hpx::init_params init_params;
	std::string sizeString = "hpx.stacks.small_size=" + std::to_string(size_t(1) << size_t(22)) + "\n";
	init_params.cfg.push_back("hpx.commandline.allow_unknown=1");
	init_params.cfg.push_back(sizeString);
	return hpx::init(argc, argv, init_params);
}
