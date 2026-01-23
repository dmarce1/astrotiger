/******************************************************************************
 Copyright (C) 202#include "IntegerPartition.hpp"
 4  Dominic C. Marcello
 *******************************************************************************/
#include "LinearCombination.hpp"
#include "Matrix.hpp"
#include "Permutation.hpp"
#include "SparseVector.hpp"
#include "Tensor.hpp"

#include <hpx/hpx_init.hpp>
#include <unordered_map>

int hpx_main(int argc, char *argv[]) {
	using namespace Tensors;
	constexpr auto Λ = IntegerPartition<3> { };
	constexpr int R = Λ.size();
	constexpr int D = 3;
	constexpr int N = ipow(D, R);

	constexpr auto p = createTensorSymmetry<Λ, D>();
	auto const A = SparseMatrix<Rational>(std::get<0>(p), N, std::get<1>(p)) / Rational(factorial(R));
	auto const iA = Rational(factorial(R)) *  SparseMatrix<Rational>(N, std::get<0>(p), std::get<2>(p));
	std::cout << A;
	std::cout << iA;
	rankReduce(iA);
	std::cout << iA;
	//Tensor<double, 3, 3> T;
//	std::cout << pseudoinverse(rankReduce<T, N / (D * D), N / (D * D), B.literal()>());
//	constexpr auto A = genTransform<sym>();
//	std::cout << A << std::endl;

//	constexpr SpechtModule<IntegerPartition<2>{}> sm2{};
//	constexpr SpechtModule<IntegerPartition<1, 1>{}> sm11{};
//	constexpr auto eT = sm2.template operator()<0>();
//	constexpr auto eT2 = directProduct<R, R, eT.size(), eT.size(), eT, eT>();
//	constexpr auto A = toMatrix<D>(eT);
//	constexpr auto B = rankReduce<T, N, N, A.literal()>();
//	std::cout << B << std::endl;
//	std::cout << factorial(R) * pseudoinverse(B) << std::endl;
//	constexpr int N = ipow(D, R);
//	constexpr SpechtModule<2, 2> sm22 { };
//	constexpr auto eT = sm22.template operator()<0>();
//	constexpr auto A = toMatrix < D > (eT);
//	constexpr auto B = rankReduce<T, N, N, A.literal()>();
//	std::cout << factorial(R) * pseudoinverse(B) << std::endl;
	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
	hpx::init_params init_params;
	std::string sizeString = "hpx.stacks.small_size=" + std::to_string(size_t(1) << size_t(22)) + "\n";
	init_params.cfg.push_back("hpx.commandline.allow_unknown=1");
	init_params.cfg.push_back(sizeString);
	return hpx::init(argc, argv, init_params);
}
