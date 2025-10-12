#include "numbers.hpp"
#include "units.hpp"
#include "BellPolynomial.hpp"
#include "MultiPrecision.hpp"
#include "DoubleReal.hpp"
#include "Integrate.hpp"
#include "FermiDirac.hpp"
#include "Quadrature.hpp"
#include "AutoDiff.hpp"
#include "Polynomial.hpp"
#include "Constants.hpp"
#include <algorithm>
#include <memory>
#include <numeric>
#include <type_traits>
#include <vector>
#include <hpx/hpx_init.hpp>
#include "Octogrid.hpp"
#include "HyperSubgrid.hpp"
#include "dgTransforms.hpp"
#include "Options.hpp"
#include "MultiIndex.hpp"
#include "EulerState.hpp"
#include "Real.hpp"
#include "RungeKutta.hpp"
#include "BiconjugateGradient.hpp"
#include "Polynomial.hpp"
#include <cassert>
void testRadiation();

#include <tuple>
#include <utility>
#include <cstddef>
#include <unordered_map>
#include <ranges>
#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <array>
#include <cstddef>
#include <functional>

#include <array>
#include <cmath>
#include <algorithm>
#include <span>


constexpr int P = 3;
constexpr int D = 2;
constexpr int N = 64;
using T = double;
using RK = RungeKutta<T, P>::type;

using SubgridType = HyperSubgrid<T, D, N, P, RK, EulerStateHLLC>;
using OctogridServerType = OctogridServer<SubgridType>;
using OctogridType = hpx::components::component<OctogridServerType>;

HPX_REGISTER_COMPONENT_MODULE();
HPX_REGISTER_COMPONENT (OctogridType);

using GetFaceChildrenAction = typename OctogridServerType::refineAction;
using RefineAction = typename OctogridServerType::getFaceChildrenAction;
using SetAction = typename OctogridServerType::setAction;

HPX_REGISTER_ACTION(GetFaceChildrenAction, getFaceChildrenAction);
HPX_REGISTER_ACTION(RefineAction, refineAction);
HPX_REGISTER_ACTION(SetAction, setAction);

struct FDTest {
	double k;
	double η;
	double β;
	double FD;
};



void polytest();
void radiation_test();

int hpx_main(int argc, char *argv[]) {
	radiation_test();
//	for (Real ρ = 1e10; ρ < Real(1e20); ρ *= Real(10)) {
//		printf( " ρ          T          n          ε          p\n");
//		for (Real T = 10; T < Real(1e13); T *= Real(10)) {
//			printf("%e %e %e %e %e\n", (double) eos.ρ, (double) eos.T, (double) eos.n, (double) eos.e, (double) eos.p);
//		}
//		printf("\n");
//	}

//	RadiationSource<double, 3> test;
//	auto tmp = bellPolynomial<6>(6, 3);
//	std::cout << tmp << "\n";

//	std::cout << "/**************************************/\n";
//	fluxExpression<3>(0, { 0, 0, 0, 0, 0 });
//	std::cout << "/**************************************/\n";
//	fluxExpression<3>(0, { 2, 0, 0, 0, 0 });
//	std::cout << "/**************************************/\n";
//	fluxExpression<3>(0, { 0, 2, 0, 0, 0 });
//	std::cout << "/**************************************/\n";
//	fluxExpression<3>(0, { 0, 0, 2, 0, 0 });
//	std::cout << "/**************************************/\n";
//	fluxExpression<3>(0, { 0, 0, 0, 2, 0 });
//	std::cout << "/**************************************/\n";
//	fluxExpression<3>(0, { 0, 0, 0, 0, 2 });
//	std::cout << "/**************************************/\n";
//	constexpr int N = 3;
//	Auto<double> x(10.0, 1.0);
//	Auto<double> dfxdx(10.0, 1.0);
//	auto f = [](auto x) {
//		return x * x * x;
//	};
//	auto g = derivative2(f, 10.0);
//	printf( "%e \n", derivative3(f, 10.0).derivative().derivative().value().value());

//	SymbolicVariable<0> rho("rho");
//	SymbolicVariable<1> S("S");
//	SymbolicVariable<2> E("E");
//	SymbolicConstant gamma(5./3.);
//	SymbolicOne one;
//	SymbolicConstant half(0.5);Constant
//	auto const p = (gamma - one) * (E - half * S * S / rho);
//	auto const U = std::make_tuple(S, S * S / rho + p, (E + p) * S / rho);
//	auto J = jacobian(U);
//	std::cout << jacobianToString(J) << "\n";
//	std::cout << "S = " << static_cast<std::string>(S) << "\n";
//	processOptions(argc, argv);
//	printf("\nPrologue complete\n");
//	RK const rk;
//	SubgridType grid;
//	grid.initialize(initSodShockTube<T, D>);
//	grid.applyLimiter();
//	grid.output("X", 0, Real(0.0));
//	T t = T(0);
//	T tmax = T(.125);
//	T dt;
//	int iter = 0;
//	while (t < tmax) {
//		grid.output("X", iter, t);
//		std::cout << "i = " << std::to_string(iter);
//		std::cout << "  t = " << std::to_string(t);
//		dt = grid.beginStep();
//		std::cout << "  dt = " << dt << std::endl;
//		for (int s = 0; s < rk.stageCount(); s++) {
//			grid.subStep(dt, s);
//		}
//		grid.endStep();
//		iter++;
//		t += dt;
//	}
	printf("\nStopping\n");
	return hpx::local::finalize();
}

thread_local FpeThreadInit _fpeThreadInitGuard;

int main(int argc, char *argv[]) {
	installFpeHandler();
	_fpeThreadInitGuard.touch();
	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };
	cfg.push_back("hpx.stacks.small_size=1048576");
	hpx::init_params init_params;
	init_params.cfg = std::move(cfg);
	auto rc = hpx::init(argc, argv, init_params);
	return rc;
}
