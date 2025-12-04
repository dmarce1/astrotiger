/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#include "Definitions.hpp"

#include "Box.hpp"
#include "Child.hpp"
#include "Definitions.hpp"
#include "Face.hpp"
#include "FixedPrecision.hpp"
#include "Gas.hpp"
#include "Integer.hpp"
#include "IO.hpp"
#include "Leaf.hpp"
#include "Math.hpp"
#include "MortonKey.hpp"
#include "Multidices.hpp"
#include "Octree.hpp"
#include "Options.hpp"
#include "Permutation.hpp"
#include "Point.hpp"
#include "Tensor.hpp"
#include <hpx/hpx_init.hpp>
#include <functional>

double test(int, int) {
	return 0.0;
}

unsigned findKey(std::vector<int> data) {
	for (unsigned m = rand();; m = rand()) {
		unsigned a = rand();
		std::vector<bool> touched(data.size(), false);
		bool flag = true;
		for (auto d : data) {
			if (touched[(d * m + a) % data.size()]) {
				flag = false;
				break;
			}
			touched[(d * m + a) % data.size()] = true;
		}
		if (flag) {
			printf("%i %i \n", m, a);
			for (auto d : data) {
				printf("%i %i\n", d, (d * m + a) % data.size());
			}
			return m;
		}
	}
}

int hpx_main(int argc, char *argv[]) {
	Options options(argc, argv);
	using namespace Tensors;
	std::vector<int> data;
	data.push_back('i');
	data.push_back('j');
	data.push_back('l');
	data.push_back('m');
//	findKey(data);
	constexpr FreeIndex<'i'> i;
	constexpr FreeIndex<'j'> j;
	constexpr FreeIndex<'k'> k;
	printf("%i\n", CommonIndices<0, 1, 2, 3>::template type<3, 0, 5, 5>::value);
//	printf("%i %i\n", 0, std::get<0>(tup)());
//	printf("%i %i\n", 1, std::get<1>(tup)());
//	printf("%i %i\n", 2, std::get<2>(tup)());
//	printf("%i %i\n", 3, std::get<3>(tup)());
//	printf("%i %i\n", 4, std::get<4>(tup)());
//	printf("%i %i\n", 5, std::get<5>(tup)());
//	for(int i = 0; i < A.size(); i++) {
//		printf( "%i %i\n", i, A[i]);
//	}
	Tensor<double, 2, 3> A;
	A(0, 0) = A(0, 1) = A(0, 2) = 1.0;
	A(1, 0) = A(1, 1) = A(1, 2) = 1.0;
	A(2, 0) = A(2, 1) = A(2, 2) = 0.0;

	double B = A(j, j);
	printf("%e\n", B);
	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
	hpx::init_params init_params;
	std::string sizeString = "hpx.stacks.small_size=" + std::to_string(size_t(1) << size_t(18)) + "\n";
	init_params.cfg.push_back("hpx.commandline.allow_unknown=1");
	init_params.cfg.push_back(sizeString);
	return hpx::init(argc, argv, init_params);
}
