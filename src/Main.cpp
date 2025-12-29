/******************************************************************************
 Copyright (C) 202#include "IntegerPartition.hpp"
4  Dominic C. Marcello
 *******************************************************************************/
#include "SparseMatrix.hpp"
#include "YoungTableau.hpp"
#include <bit>

#include "Definitions.hpp"
#include "SymmetricGroup.hpp"

#include "Approximate.hpp"
#include "Box.hpp"
#include "Child.hpp"
#include "Definitions.hpp"
#include "Face.hpp"
#include "FixedPrecision.hpp"
#include "Gas.hpp"
#include "IO.hpp"
#include "Integer.hpp"
#include "Leaf.hpp"
#include "Math.hpp"
#include "Matrix.hpp"
#include "MortonKey.hpp"
#include "Multidices.hpp"
#include "Octree.hpp"
#include "Options.hpp"
#include "Point.hpp"
#include "Rational.hpp"
#include "Tensor.hpp"
#include "Vector.hpp"
#include <functional>
#include <hpx/hpx_init.hpp>

int hpx_main(int argc, char *argv[]) {
	//	using TupleType = std::tuple<MatrixElement<0, 0, 1.0>, MatrixElement<1, 1, 1.0>, MatrixElement<2, 2, 1.0>>;
	//	printf("%i\n", FindMatrixElement<TupleType, 0, 0>::value);
	//	printf("%i\n", FindMatrixElement<TupleType, 0, 1>::value);
	//	printf("%i\n", FindMatrixElement<TupleType, 0, 2>::value);
	//	printf("%i\n", FindMatrixElement<TupleType, 1, 0>::value);
	//	printf("%i\n", FindMatrixElement<TupleType, 1, 1>::value);
	//	printf("%i\n", FindMatrixElement<TupleType, 1, 2>::value);
	//	printf("%i\n", FindMatrixElement<TupleType, 2, 0>::value);
	//	printf("%i\n", FindMatrixElement<TupleType, 2, 1>::value);
	//	printf("%i\n", FindMatrixElement<TupleType, 2, 2>::value);
	//	printf("%i \n", YoungTableau<1, 1, 1>::rowPermutationCount());
	//	printf("%i \n", YoungTableau<1, 1, 1>::columnPermutationCount());
	//	printf("%i \n", YoungTableau<1, 1, 1>::permutationCount());
	//	Matrix<double, 20, 20> A(std::array<std::array<double, 20>, 20>{
	//		{{{21.4, 0.3, -1.1, 0.7, -0.4, 1.2, 0.0, -0.9, 0.6, -0.8, 0.5, 1.1, -0.2, 0.4, -0.6, 0.9, -1.0, 0.2, -0.3, 0.8}},
	//		 {{-0.5, 22.1, 0.6, -1.2, 0.4, -0.7, 1.0, 0.3, -0.9, 0.5, -0.4, 0.8, 1.1, -0.6, 0.2, -0.3, 0.7, -0.8, 0.9, -0.1}},
	//		 {{0.7, -0.4, 23.0, 0.5, -1.3, 0.6, -0.8, 1.2, 0.1, -0.9, 0.4, -0.5, 0.3, 0.8, -0.6, 1.0, -0.7, 0.9, -0.2, 0.4}},
	//		 {{-0.6, 0.8, -0.3, 24.2, 0.5, -1.1, 0.9, -0.4, 0.7, 0.2, -0.8, 0.6, -0.9, 1.0, 0.4, -0.5, 0.3, -0.7, 0.1, 0.8}},
	//		 {{0.4, -0.7, 0.6, -0.5, 25.3, 0.9, -1.0, 0.8, -0.2, 0.7, 0.1, -0.4, 1.2, -0.6, 0.5, -0.3, 0.9, -0.8, 0.4, -0.1}},
	//		 {{-0.9, 0.5, -0.6, 0.8, -0.4, 26.1, 0.7, -1.2, 0.9, -0.3, 0.6, -0.8, 0.4, 0.5, -0.1, 1.0, -0.7, 0.2, -0.4, 0.9}},
	//		 {{0.3, -0.9, 0.4, -0.7, 0.6, -0.5, 27.2, 0.8, -1.1, 0.9, -0.2, 0.7, -0.6, 0.4, 1.0, -0.8, 0.5, -0.3, 0.9, -0.1}},
	//		 {{-0.8, 0.4, -0.9, 0.6, -0.7, 0.3, 0.5, 28.0, 0.8, -1.2, 0.4, -0.6, 0.9, -0.3, 0.7, 0.1, -0.5, 1.0, -0.2, 0.6}},
	//		 {{0.6, -0.3, 0.8, -0.4, 0.2, -0.9, 0.7, 0.5, 29.1, 0.6, -1.0, 0.4, -0.8, 0.9, -0.5, 0.7, 0.3, -0.6, 1.0, -0.2}},
	//		 {{-0.7, 0.5, -0.2, 0.3, -0.8, 0.4, -0.9, 0.7, 0.6, 30.4, 0.8, -1.1, 0.5, -0.4, 0.9, -0.6, 0.7, 0.2, -0.3, 1.0}},
	//		 {{0.5, -0.4, 0.3, -0.8, 0.1, -0.6, 0.2, -0.4, -0.9, 0.7, 31.3, 0.8, -1.0, 0.6, -0.5, 0.9, -0.3, 0.7, 0.4, -0.2}},
	//		 {{-0.4, 0.6, -0.5, 0.7, -0.3, 0.8, -0.6, 0.5, 0.4, -0.9, 32.2, 0.9, 0.7, -1.1, 0.6, -0.5, 0.8, -0.2, 0.4, 0.3}},
	//		 {{0.2, -0.7, 0.3, -0.9, 1.0, -0.4, 0.6, -0.8, 0.5, 0.4, -0.6, 33.5, 0.9, 0.7, -1.2, 0.8, -0.5, 0.6, -0.3, 0.4}},
	//		 {{-0.3, 0.8, -0.6, 1.0, -0.5, 0.4, -0.7, 0.3, 0.9, -0.4, 0.6, 0.8, 34.1, -1.0, 0.7, -0.6, 0.5, 0.9, -0.2, 0.4}},
	//		 {{0.6, -0.5, 0.8, -0.4, 0.7, -0.1, 0.9, 0.6, -0.5, 0.8, -0.6, 0.4, 0.7, 35.0, -1.1, 0.9, -0.8, 0.5, 0.6, -0.3}},
	//		 {{-0.7, 0.3, -0.6, 0.5, -0.4, 1.0, -0.8, 0.1, 0.7, -0.6, 0.9, -0.5, 0.8, 0.9, 36.2, -1.0, 0.6, -0.4, 0.5, 0.3}},
	//		 {{0.8, -0.6, 0.7, -0.3, 0.9, -0.7, 0.5, -0.4, 0.3, 0.7, -0.3, 0.8, -0.5, 0.6, 0.7, 37.1, -1.2, 0.9, -0.4, 0.5}},
	//		 {{-1.0, 0.7, -0.9, 0.3, -0.8, 0.2, -0.3, 1.0, -0.6, 0.2, 0.7, -0.2, 0.9, -0.4, -0.6, 0.9, 38.0, -1.1, 0.8, -0.3}},
	//		 {{0.4, -0.8, 0.9, -0.7, 0.4, -0.4, 0.9, -0.2, 1.0, -0.3, 0.4, 0.6, -0.3, 0.6, 0.5, -0.4, 0.8, 39.2, -1.0, 0.7}},
	//		 {{-0.6, 0.9, -0.2, 0.8, -0.1, 0.9, -0.1, 0.6, -0.2, 1.0, -0.2, 0.3, 0.4, -0.3, 0.3, 0.5, -0.3, 0.7, 0.8, 40.5}}}});
	//	std::cout << A;
	//	auto I = decltype(A)::identity();
	//	auto const [_, iA] = gaussianElimination<GEType::full>(A, I);
	//	std::cout << iA;
	// std::cout << A * iA;
	constexpr auto yS = YoungSymmetrizer<3, 2, 1>();
	std::cout << yS << std::endl;
	//	constexpr int N = 10;
	//	constexpr Matrix<Rational, 4> A1 = { { 1, 2, 3, 4}, { -2, 3, 0,5 }, { 1, -1, 1, 6 }, {2, 4, 6, 8} };
	//	constexpr Matrix<double, 4> A2 =  { { 1, 2, 3, 4}, { -2, 3, 0,5 }, { 1, -1, 1, 6 }, {2, 4, 6, 8} };
	//	std::cout << A1 << std::endl;
	//	std::cout << A2 << std::endl;
	//	constexpr auto A1lit = (MatrixLiteral<Rational, 4, 4> ) A1;
	//	constexpr auto A2lit = (MatrixLiteral<double, 4, 4> ) A2;
	//	std::cout << rankReduce<Rational, 4, 4, A1lit>() << std::endl;
	//	std::cout << rankReduce<double, 4, 4, A2lit>() << std::endl;
	//	std::cout <<  pseudoinverse(rankReduce<Rational, 4, 4, A1lit>()) << std::endl;
	//	std::cout <<  pseudoinverse(rankReduce<double, 4, 4, A2lit>()) << std::endl;
	//	constexpr Permutation<N> P = {9, 2, 3, 7, 4, 5, 0, 8, 1, 6};
	//	std::cout << P << std::endl;
	//	auto const cycles = P.cycles();
	//	for (auto const &cycle : cycles) {
	//		std::cout << cycle;
	//	}
	// x x x
	// x x
	// x x
	//	std::cout << std::endl;
	//	static constexpr YoungSymmetrizer<2, 2> Y{};
	//	std::cout << Y.Rt << std::endl;
	//	std::cout << Y << std::endl;
	//	static auto constexpr A = Y.template genMatrix<4>();
	//	//	std::cout << A;
	//	static auto constexpr Alit = (MatrixLiteral<double, A.rowCount(), A.columnCount()>)A;
	//	static auto constexpr rA = rankReduce<double, A.rowCount(), A.columnCount(), Alit>();
	//	std::cout << rA << std::endl;
	//	std::cout << std::endl;
	//	std::cout << pseudoinverse(rA) << std::endl;
	//	//	std::cout << std::get<1>(Y.Sc) << std::endl;
	//	//	using namespace Tensors;
	//	//	constexpr int D = 3;
	//	constexpr int R = 3;
	//	TensorSymmetry<R, D> sym;
	//	constexpr Cycle<2> cyc1 = {0, 1};
	//	constexpr Cycle<2> cyc2 = {0, 2};
	//	constexpr Cycle<2> cyc3 = {1, 2};
	//	constexpr Cycle<3> cyc4 = {0, 1, 2};
	//	constexpr Cycle<3> cyc5 = {0, 2, 1};
	//	constexpr TensorSymmetry<R, D> sym0{};
	//	constexpr TensorSymmetry<R, D> sym1(cyc1);
	//	constexpr TensorSymmetry<R, D> sym2(cyc2);
	//	constexpr TensorSymmetry<R, D> sym3(cyc3);
	//	constexpr TensorSymmetry<R, D> sym4(cyc4);
	//	constexpr TensorSymmetry<R, D> sym5(cyc5);
	//	constexpr auto sym = sym0 + sym1 + sym2 + sym3 + sym4 + sym5;
	//	std::cout << sym;
	//	std::cout << rankReduce<Rational<int>, ipow(D, R), ipow(D, R), MatrixLiteral<Rational<int>, ipow(D, R)>(sym)>();
	//	std::cout << pseudoinverse(rankReduce<Rational<int>, ipow(D, R), ipow(D, R), MatrixLiteral<Rational<int>, ipow(D, R)>(sym)>());
	//	constexpr Matrix<double, D * D> A = {{2, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 1, 0, 0},
	//										 {0, 1, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 2, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 1, 0},
	//										 {0, 0, 1, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 2}};
	////	constexpr Matrix<double, D * D> A = {{2, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0},
	////										 {0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 2, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0, 0},
	////										 {0, 0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 2}};
	//	constexpr auto I = Matrix<double, D * D>::identity();
	//	constexpr auto Alit = MatrixLiteral<double, D * D>(A);
	//	constexpr auto R = rankReduce<double, D * D, D * D, Alit>();
	//	auto iR = pseudoinverse(R);
	//
	//	std::cout << A << std::endl;
	//	std::cout << I << std::endl;
	//	std::cout << R << std::endl;
	//	std::cout << iR << std::endl;
	//	std::cout << std::to_string(rank<double, D * D, D * D, Alit>()) << std::endl;
	//	YoungTableau<int, 4, 2, 2> tab;
	//	auto const otab = tab;
	//	for (int n = 0; n < 21; n++) {
	//		std::cout << std::to_string(n + 1) << ".\n";
	//		std::cout << tab << std::endl;
	//		tab++;
	//		if (tab == otab) {
	//			break;
	//		}
	//	}
	//	std::cout << std::to_string(tab.count()) << "\n";
	//	FreeIndex<'i'> i;
	//	FreeIndex<'j'> j;
	//	FreeIndex<'k'> k;
	//	FreeIndex<'l'> l;
	//	Tensor<double, 2, 3> B;
	//	B(i, i, 3);
	//	//	Options options(argc, argv);
	//	//	using namespace Tensors;
	//	std::vector<int> data;
	//	data.push_back('i');
	//	data.push_back('j');
	//	data.push_back('l');
	//	data.push_back('m');
	////	findKey(data);
	////	constexpr FreeIndex<'i'> i;
	////	constexpr FreeIndex<'j'> j;
	////	constexpr FreeIndex<'k'> k;
	////	auto tup = CommonIndices<0, 1, 2, 3>::template type<3, 0, 5, 6>::value;
	////	printf("%i %i\n", 0, std::get<0>(tup)());
	////	printf("%i %i\n", 1, std::get<1>(tup)());
	//////	printf("%i %i\n", 2, std::get<2>(tup)());
	//////	printf("%i %i\n", 3, std::get<3>(tup)());
	//////	printf("%i %i\n", 4, std::get<4>(tup)());
	//////	printf("%i %i\n", 5, std::get<5>(tup)());
	//////	for(int i = 0; i < A.size(); i++) {
	//////		printf( "%i %i\n", i, A[i]);
	//////	}
	////	Tensor<double, 2, 3> A;
	////	A(0, 0) = A(0, 1) = A(0, 2) = 1.0;
	////	A(1, 0) = A(1, 1) = A(1, 2) = 1.0;
	////	A(2, 0) = A(2, 1) = A(2, 2) = 0.0;
	////
	////	double B = A(j, j);
	////	printf("%e\n", B);
	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
	hpx::init_params init_params;
	std::string sizeString = "hpx.stacks.small_size=" + std::to_string(size_t(1) << size_t(18)) + "\n";
	init_params.cfg.push_back("hpx.commandline.allow_unknown=1");
	init_params.cfg.push_back(sizeString);
	return hpx::init(argc, argv, init_params);
}
