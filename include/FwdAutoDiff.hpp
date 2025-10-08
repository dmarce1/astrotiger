///******************************************************************************
// Cartesian Multivariate Forward AutoDiff (Jets)
// — Heterogeneous tuple storage (homogeneous placeholder)
// Author: Dominic C. Marcello (2025)
// *******************************************************************************/
//#ifndef INCLUDE_FWDAUTODIFF_JET_HPP_
//#define INCLUDE_FWDAUTODIFF_JET_HPP_
//
//#include <tuple>
//#include <array>
//#include <utility>
//#include <type_traits>
//#include <functional>
//#include <cmath>
//#include <ostream>
//
//// ============================================================
//// Forward declaration and traits
//// ============================================================
//template<typename, int, class ...>
//struct FwdAutoDiff;
//
//template<typename T>
//struct IsFwdAutoDiff: std::false_type {
//};
//
//template<class Tdep, int Order, class ... Ts>
//struct IsFwdAutoDiff<FwdAutoDiff<Tdep, Order, Ts...>> : std::true_type {
//};
//
//template<typename T>
//inline constexpr bool IsFwdAutoDiff_v = IsFwdAutoDiff<T>::value;
//
//// ============================================================
//// Compile-time Cartesian index iterator (0..Order inclusive)
//// ============================================================
//template<int Dim, int Order, typename Func>
//constexpr void forEachIndex(std::array<int, Dim> &idx, Func &&f, int depth = 0) {
//	if constexpr (Dim == 0) {
//		f();
//	} else {
//		if (depth == Dim - 1) {
//			for (idx[depth] = 0; idx[depth] <= Order; ++idx[depth])
//				f(idx);
//		} else {
//			for (idx[depth] = 0; idx[depth] <= Order; ++idx[depth])
//				forEachIndex<Dim, Order>(idx, std::forward < Func > (f), depth + 1);
//		}
//	}
//}
//
//template<int Dim, int Order, typename Func>
//constexpr void forEachIndex(Func &&f) {
//	std::array<int, Dim> idx { };
//	forEachIndex<Dim, Order>(idx, std::forward < Func > (f), 0);
//}
//
//template<int Dim>
//constexpr std::size_t flatten(std::array<int, Dim> const &i, int order) {
//	std::size_t idx = 0, stride = 1;
//	for (int k = 0; k < Dim; ++k) {
//		idx += static_cast<std::size_t>(i[k]) * stride;
//		stride *= static_cast<std::size_t>(order + 1);
//	}
//	return idx;
//}
//
//// Small helpers for product/quotient rules
//template<int Dim>
//constexpr auto addIndex(std::array<int, Dim> a, std::array<int, Dim> const &b) {
//	for (int k = 0; k < Dim; ++k)
//		a[k] += b[k];
//	return a;
//}
//template<int Dim>
//constexpr auto subIndex(std::array<int, Dim> a, std::array<int, Dim> const &b) {
//	for (int k = 0; k < Dim; ++k)
//		a[k] -= b[k];
//	return a;
//}
//template<int Dim>
//constexpr bool leqIndices(std::array<int, Dim> const &a, std::array<int, Dim> const &b) {
//	for (int k = 0; k < Dim; ++k)
//		if (a[k] > b[k])
//			return false;
//	return true;
//}
//template<int Dim>
//constexpr bool insideBox(std::array<int, Dim> const &i) {
//	for (int k = 0; k < Dim; ++k)
//		if (i[k] < 0)
//			return false;
//	return true;
//}
//template<int Dim>
//constexpr bool isZeroIndex(std::array<int, Dim> const &i) {
//	for (int k = 0; k < Dim; ++k)
//		if (i[k] != 0)
//			return false;
//	return true;
//}
//template<int Dim>
//constexpr int sumIndex(std::array<int, Dim> const &i) {
//	int s = 0;
//	for (int k = 0; k < Dim; ++k)
//		s += i[k];
//	return s;
//}
//
//// ============================================================
//// FwdAutoDiff<Tdep, Order, Ts...> — Jet implementation
//// ============================================================
//template<typename Tdep, int Order, typename ... Ts>
//struct FwdAutoDiff {
//	static constexpr int dim = sizeof...(Ts);
//	static constexpr int order = Order;
//
//	static constexpr std::size_t storageSize = [] {
//		std::size_t n = 1;
//		for (int i = 0; i < dim; ++i)
//			n *= (Order + 1);
//		return n;
//	}();
//
//	using Storage = std::tuple<std::array<Tdep, storageSize>>;
//
//	Storage C_;
//
//	// ---------------------------------------------------------
//	// Construction
//	// ---------------------------------------------------------
//	constexpr FwdAutoDiff() {
//		std::get < 0 > (C_).fill(Tdep(0));
//	}
//	explicit constexpr FwdAutoDiff(Tdep const &val) {
//		std::get < 0 > (C_).fill(Tdep(0));
//		std::get < 0 > (C_)[0] = val;
//	}
//
//	constexpr auto& raw() {
//		return std::get < 0 > (C_);
//	}
//	constexpr auto const& raw() const {
//		return std::get < 0 > (C_);
//	}
//
//	// ---------------------------------------------------------
//	// Element access
//	// ---------------------------------------------------------
//	constexpr Tdep& operator[](auto const &idx) {
//		if constexpr (dim == 1 && std::is_integral_v<std::decay_t<decltype(idx)>>)
//			return raw()[static_cast<std::size_t>(idx)];
//		else
//			return raw()[flatten<dim>(idx, order)];
//	}
//	constexpr Tdep const& operator[](auto const &idx) const {
//		if constexpr (dim == 1 && std::is_integral_v<std::decay_t<decltype(idx)>>)
//			return raw()[static_cast<std::size_t>(idx)];
//		else
//			return raw()[flatten<dim>(idx, order)];
//	}
//
//	explicit constexpr operator Tdep() const {
//		return raw()[0];
//	}
//
//	// ---------------------------------------------------------
//	// Independent variable construction
//	// ---------------------------------------------------------
//	static constexpr FwdAutoDiff independentVariable(Tdep const &x0, int axis = 0) {
//		FwdAutoDiff a;
//		a.raw()[0] = x0;
//		if (axis >= 0 && axis < dim) {
//			std::array<int, dim> e { };
//			e[axis] = 1;
//			a[e] = Tdep(1);
//		}
//		return a;
//	}
//
//	// ---------------------------------------------------------
//	// Compose helper (Taylor composition)
//	// ---------------------------------------------------------
//	template<int dorder, typename Function>
//	static constexpr void composeHelper(Function const &dkdF, Tdep const &g0, FwdAutoDiff const &delta, FwdAutoDiff &h, FwdAutoDiff &term) {
//		if constexpr (dorder > Order)
//			return;
//		else {
//			auto const a_n = dkdF(g0, dorder) / factorial<double>(dorder);
//			for (std::size_t i = 0; i < storageSize; ++i)
//				h.raw()[i] += a_n * term.raw()[i];
//			for (std::size_t i = 0; i < storageSize; ++i)
//				term.raw()[i] *= delta.raw()[i];
//			composeHelper<dorder + 1>(dkdF, g0, delta, h, term);
//		}
//	}
//
//	template<typename Function>
//	friend constexpr auto compose(Function const &dkdF, FwdAutoDiff const &g) {
//		FwdAutoDiff h;
//		Tdep const g0 = g.raw()[0];
//		FwdAutoDiff delta = g;
//		delta.raw()[0] = Tdep(0);
//		FwdAutoDiff term(Tdep(1));
//		composeHelper<0>(dkdF, g0, delta, h, term);
//		return h;
//	}
//
//	// ---------------------------------------------------------
//	// Elementary example: exp
//	// ---------------------------------------------------------
//	friend constexpr auto exp(FwdAutoDiff const &x) {
//		return compose([](Tdep const &g0, int) {
//			using std::exp;
//			return exp(g0);
//		}, x);
//	}
//
//	// ---------------------------------------------------------
//	// Output
//	// ---------------------------------------------------------
//	friend std::ostream& operator<<(std::ostream &os, FwdAutoDiff const &A) {
//		forEachIndex<dim, order>([&](auto const &idx) {
//			os << "[";
//			for (int k = 0; k < dim; ++k)
//				os << idx[k] << (k + 1 < dim ? "," : "");
//			os << "]=" << A[idx] << " ";
//		});
//		return os;
//	}
//
//private:
//	template<typename T>
//	static constexpr T factorial(int n) {
//		T v { 1 };
//		for (int k = 1; k <= n; ++k)
//			v *= T(k);
//		return v;
//	}
//};
//
//// ============================================================
//// Jet × Jet (heterogeneous units)
//// ============================================================
//template<class A, int O, class ...TA, class B, class ...TB, class R = decltype(std::declval<A>() * std::declval<B>())>
//constexpr auto operator*(FwdAutoDiff<A, O, TA...> const &X, FwdAutoDiff<B, O, TB...> const &Y) {
//	FwdAutoDiff<R, O, TA...> Z;
//	forEachIndex<Z.dim, Z.order>([&](auto const &i) {
//		forEachIndex<Z.dim, Z.order>([&](auto const &j) {
//			auto s = addIndex<Z.dim>(i, j);
//			if (insideBox<Z.dim>(s))
//				Z[s] += X[i] * Y[j];
//		});
//	});
//	return Z;
//}
//
//template<class A, int O, class ...TA, class B, class ...TB, class R = decltype(std::declval<A>() / std::declval<B>())>
//constexpr auto operator/(FwdAutoDiff<A, O, TA...> const &X, FwdAutoDiff<B, O, TB...> const &Y) {
//	FwdAutoDiff<R, O, TA...> Z;
//	auto const y0 = Y.raw()[0];
//	Z.raw()[0] = X.raw()[0] / y0;
//	forEachIndex<Z.dim, Z.order>([&](auto const &I) {
//		if (sumIndex<Z.dim>(I) == 0)
//			return;
//		auto acc = X[I];
//		forEachIndex<Z.dim, Z.order>([&](auto const &K) {
//			if (isZeroIndex<Z.dim>(K) || !leqIndices<Z.dim>(K, I))
//				return;
//			auto ImK = subIndex<Z.dim>(I, K);
//			acc -= Y[K] * Z[ImK];
//		});
//		Z[I] = acc / y0;
//	});
//	return Z;
//}
//
//// ============================================================
//// Jet × scalar  (true scalar, not a jet)
//// ============================================================
//template<class A,int O,class...TA,class S>
//requires(!IsFwdAutoDiff_v<S>)
//constexpr auto operator*(FwdAutoDiff<A,O,TA...> const& X,S const& s) {
//	using R = decltype(std::declval<A>()*std::declval<S>());
//	FwdAutoDiff<R,O,TA...> Z;
//	for(std::size_t i=0;i<Z.storageSize;++i) Z.raw()[i]=X.raw()[i]*s;
//	return Z;
//}
//template<class S,class A,int O,class...TA>
//requires(!IsFwdAutoDiff_v<S>)
//constexpr auto operator*(S const& s,FwdAutoDiff<A,O,TA...> const& X) {return X*s;}
//
//template<class A,int O,class...TA,class S>
//requires(!IsFwdAutoDiff_v<S>)
//constexpr auto operator/(FwdAutoDiff<A,O,TA...> const& X,S const& s) {
//	using R = decltype(std::declval<A>()/std::declval<S>());
//	FwdAutoDiff<R,O,TA...> Z;
//	for(std::size_t i=0;i<Z.storageSize;++i) Z.raw()[i]=X.raw()[i]/s;
//	return Z;
//}
//template<class S,class A,int O,class...TA>
//requires(!IsFwdAutoDiff_v<S>)
//constexpr auto operator/(S const& s,FwdAutoDiff<A,O,TA...> const& Y) {
//	FwdAutoDiff<S,O,TA...> X(s);
//	return X/Y;
//}
//
//// ============================================================
//// Constant and oneLike helpers
//// ============================================================
//template<class Tdep, int O, class ...Ts>
//constexpr auto makeConstantLike(FwdAutoDiff<Tdep, O, Ts...> const&, Tdep const &v) {
//	FwdAutoDiff<Tdep, O, Ts...> z;
//	z.raw().fill(Tdep { });
//	z.raw()[0] = v;
//	return z;
//}
//template<class AnyJet>
//constexpr auto oneLike(AnyJet const &ref) {
//	using Dep = decltype(ref.raw()[0]/ref.raw()[0]);
//	// dimensionless
//	return makeConstantLike(ref, Dep(1));
//}
//
//// ========= Jet ± Jet =========
//template<class A, int O, class ...TA, class B, class ...TB, class R = decltype(std::declval<A>() + std::declval<B>())>
//constexpr auto operator+(FwdAutoDiff<A, O, TA...> const &X, FwdAutoDiff<B, O, TB...> const &Y) {
//	static_assert(O==O,"orders must match");
//	FwdAutoDiff<R, O, TA...> Z;
//	for (std::size_t i = 0; i < Z.storageSize; ++i)
//		Z.raw()[i] = X.raw()[i] + Y.raw()[i];
//	return Z;
//}
//
//template<class A, int O, class ...TA, class B, class ...TB, class R = decltype(std::declval<A>() - std::declval<B>())>
//constexpr auto operator-(FwdAutoDiff<A, O, TA...> const &X, FwdAutoDiff<B, O, TB...> const &Y) {
//	static_assert(O==O,"orders must match");
//	FwdAutoDiff<R, O, TA...> Z;
//	for (std::size_t i = 0; i < Z.storageSize; ++i)
//		Z.raw()[i] = X.raw()[i] - Y.raw()[i];
//	return Z;
//}
//
//// ========= Jet ± scalar =========
//template<class A,int O,class...TA,class S>
//requires(!IsFwdAutoDiff_v<S>)
//constexpr auto operator+(FwdAutoDiff<A,O,TA...> const& X,S const& s)
//{
//	using R = decltype(std::declval<A>()+std::declval<S>());
//	FwdAutoDiff<R,O,TA...> Z;
//	Z.raw() = X.raw();
//	Z.raw()[0] += s;
//	return Z;
//}
//
//template<class S,class A,int O,class...TA>
//requires(!IsFwdAutoDiff_v<S>)
//constexpr auto operator+(S const& s,FwdAutoDiff<A,O,TA...> const& X) {return X+s;}
//
//template<class A,int O,class...TA,class S>
//requires(!IsFwdAutoDiff_v<S>)
//constexpr auto operator-(FwdAutoDiff<A,O,TA...> const& X,S const& s)
//{
//	using R = decltype(std::declval<A>()-std::declval<S>());
//	FwdAutoDiff<R,O,TA...> Z;
//	Z.raw() = X.raw();
//	Z.raw()[0] -= s;
//	return Z;
//}
//
//template<class S,class A,int O,class...TA>
//requires(!IsFwdAutoDiff_v<S>)
//constexpr auto operator-(S const& s,FwdAutoDiff<A,O,TA...> const& X)
//{
//	using R = decltype(std::declval<S>()-std::declval<A>());
//	FwdAutoDiff<R,O,TA...> Z;
//	for(std::size_t i=0;i<Z.storageSize;++i) Z.raw()[i] = -X.raw()[i];
//	Z.raw()[0] += s;
//	return Z;
//}
//
//
//// ============================================================
//// Elementary functions: sqrt, log, pow<Rational>
//// ============================================================
//
//// --- sqrt(x) ------------------------------------------------
//template<class A, int O, class ... Ts>
//constexpr auto sqrt(FwdAutoDiff<A, O, Ts...> const &x) {
//	using std::sqrt;
//	using R = decltype(sqrt(std::declval<A>()));
//	FwdAutoDiff<R, O, Ts...> y;
//	R const x0 = x.raw()[0];
//	R const y0 = sqrt(x0);
//	y.raw()[0] = y0;
//
//	if constexpr (O >= 1) {
//		auto inv2y = R(0.5) / y0;
//		for (std::size_t i = 1; i < y.storageSize; ++i)
//			y.raw()[i] = inv2y * x.raw()[i];
//	}
//	return y;
//}
//
//// --- log(x) -------------------------------------------------
//template<class A, int O, class ... Ts>
//constexpr auto log(FwdAutoDiff<A, O, Ts...> const &x) {
//	using std::log;
//	using R = decltype(log(std::declval<A>()));
//	FwdAutoDiff<R, O, Ts...> y;
//	R const x0 = x.raw()[0];
//	y.raw()[0] = log(x0);
//
//	if constexpr (O >= 1) {
//		auto invx = R(1) / x0;
//		for (std::size_t i = 1; i < y.storageSize; ++i)
//			y.raw()[i] = invx * x.raw()[i];
//	}
//	return y;
//}
//
//template<Rational Power, class A, int O, class... Ts>
//constexpr auto pow(FwdAutoDiff<A, O, Ts...> const& x)
//{
//    using std::pow;
//    constexpr double p = static_cast<double>(Power);
//
//    // Result type (unit-aware)
//    using R = decltype(pow<Power>(std::declval<A>()));
//
//    FwdAutoDiff<R, O, Ts...> h;
//    A const x0 = x.raw()[0];
//    R const xPow = pow<Power>(x0); // correct physical unit
//
//    h.raw().fill(R(0));
//    h.raw()[0] = xPow; // base term
//
//    if constexpr (O >= 1) {
//        // Make delta/x0 dimensionless
//        FwdAutoDiff<A, O, Ts...> delta = x;
//        delta.raw()[0] = A(0);
//
//        using Dimless = Quantity<Unit<0, 0, 0, 0, 0>, double>;
//        FwdAutoDiff<Dimless, O, Ts...> eps;
//        for (std::size_t i = 0; i < eps.storageSize; ++i) {
//            eps.raw()[i] = Dimless(delta.raw()[i] / x0); // dimensionless
//        }
//
//        // Build (1 + eps)^p
//        FwdAutoDiff<Dimless, O, Ts...> term(Dimless(1));
//        FwdAutoDiff<Dimless, O, Ts...> sum(Dimless(1));
//
//        auto fact = [](int n) constexpr {
//            double v = 1.0;
//            for (int k = 1; k <= n; ++k)
//                v *= k;
//            return v;
//        };
//        auto factPow = [](double p, int n) constexpr {
//            double r = 1.0;
//            for (int k = 0; k < n; ++k)
//                r *= (p - k);
//            return r;
//        };
//
//        for (int n = 1; n <= O; ++n) {
//            double const coeff = factPow(p, n) / fact(n);
//            for (std::size_t i = 0; i < sum.storageSize; ++i)
//                sum.raw()[i] += Dimless(coeff) * term.raw()[i];
//            for (std::size_t i = 0; i < term.storageSize; ++i)
//                term.raw()[i] *= eps.raw()[i];
//        }
//
//        // Now multiply back the unit-aware xPow
//        for (std::size_t i = 0; i < h.storageSize; ++i)
//            h.raw()[i] = xPow * sum.raw()[i];
//    }
//
//    return h;
//}
//
//#endif /* INCLUDE_FWDAUTODIFF_JET_HPP_ */
