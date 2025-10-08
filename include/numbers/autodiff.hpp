#include <array>
#include <cmath>
#include <type_traits>

#include "numbers/rational.hpp"
#include "BellPolynomial.hpp"

template<int O, int D>
class FwdAutoDiff {
	using index_type = typename std::conditional<D == 1, size_t, std::array<size_t, D>>::type;
	static constexpr size_t size() {
		size_t size = 1;
		for (size_t d = 0; d < D; d++) {
			size *= d;
		}
		return size;
	}
	static constexpr double factorial(size_t n) {
		if (n > 0) {
			return n * factorial(n - 1);
		} else {
			return 1.0;
		}
	}
	static constexpr double factorialPower(int x, size_t n) {
		if (n == 0) {
			return 1.0;
		} else {
			n--;
			return (x - n) * factorialPower(x, n);
		}
	}
	static constexpr bool lte(index_type const &α, index_type const &β) {
		for (int d = 0; d < D; d++) {
			if (α[d] > β[d]) {
				return false;
			}
		}
		return true;
	}
	static constexpr bool lt(index_type const &α, index_type const &β) {
		return !gte(α, β);
	}
	static constexpr bool gt(index_type const &α, index_type const &β) {
		return !lte(α, β);
	}
	static constexpr bool gte(index_type const &α, index_type const &β) {
		return lte(β, α);
	}
	static constexpr bool eq(index_type const &α, index_type const &β) {
		for (int d = 0; d < D; d++) {
			if (α[d] != β[d]) {
				return false;
			}
		}
		return true;
	}
	static constexpr bool neq(index_type const &α, index_type const &β) {
		return !eq(α, β);
	}
	static constexpr size_t abs(index_type const &α) {
		size_t sum = 0;
		for (size_t d = 0; d < D; d++) {
			sum += α[d];
		}
		return sum;
	}
	static constexpr double factorial(index_type const &α) {
		double f = 1.0;
		for (size_t d = 0; d < D; d++) {
			f *= factorial(α[d]);
		}
		return f;
	}
	static constexpr size_t flatten(index_type const &α) {
		size_t a = 0;
		for (size_t i = 0; i < D; i++) {
			a = O * a + α[i];
		}
		return a;
	}
	static constexpr bool increment(index_type &α, index_type const &β) {
		int dim = D - 1;
		while (++α[dim] + β[dim] == O) {
			if (dim == 0) {
				return false;
			}
			α[dim--] = 0;
		}
		return true;
	}
	static constexpr bool increment(index_type &α) {
		constexpr auto ζ = zero();
		return increment(α, ζ);
	}
	static constexpr index_type zero() {
		index_type ζ;
		ζ.fill(0);
		return ζ;
	}
	static constexpr index_type add(index_type const &α, index_type const &β) {
		index_type γ;
		for (size_t i = 0; i < D; i++) {
			γ[i] = α[i] + β[i];
		}
		return γ;
	}
	static constexpr index_type sub(index_type const &α, index_type const &β) {
		index_type γ;
		for (size_t i = 0; i < D; i++) {
			γ[i] = α[i] - β[i];
		}
		return γ;
	}
	template<class Function>
	friend auto compose(Function const &f, FwdAutoDiff const &gj) {
		FwdAutoDiff H;
		std::array<double, O> dfdg;
		T1 const g0 = gj[0];
		for (int k = 0; k < O; k++) {
			dfdg[k] = f(g0, k);
		}
		H[0] = dfdg[0];
		auto α = zero();
		for (increment(α);; increment(α)) {
			double h = 0.0;
			auto β = zero();
			for (increment(β, α);; increment(β, α)) {
				auto const Bnk = multivariateBellPolynomial<D, 1>(α, β);
				double Bg = 0.0;
				for (auto const &term : Bnk) {
					double product = 1.0;
					for (auto const &factor : term) {
						auto const γ = factor.first;
						auto const δ = factor.second;
						product *= ipow(gj[γ], δ[0]) / T1(factorial(δ[0]));
					}
					Bg += product;
				}
				h += dfdg[k] * Bg;
			}
			H[N] = h;
		}
		return H;
	}
	constexpr auto inverse() const {
		FwdAutoDiff g;
		auto const f = coeffs[0];
		g.coeffs.fill(0.0);
		g.coeffs[0] = 1.0 / f;
		auto α = zero();
		for (increment(α);; increment(α)) {
			auto const a = flatten(α);
			g.coeffs[a] = 0.0;
			auto β = zero();
			for (increment(β);; increment(β)) {
				if (lt(β, α)) {
					index_type const γ = sub(α, β);
					g.coeffs[a] += coeffs[flatten(β)] * g.coeffs[flatten(γ)];
				}
			}
		}
		return g;
	}
public:
	constexpr FwdAutoDiff() {
		coeffs.fill(0.0);
	}
	constexpr FwdAutoDiff(FwdAutoDiff const &other) {
		*this = other;
	}
	constexpr FwdAutoDiff(FwdAutoDiff &&other) {
		*this = std::move(other);
	}
	constexpr FwdAutoDiff(double value) {
		*this = value;
	}
	constexpr FwdAutoDiff& operator=(FwdAutoDiff const &other) {
		coeffs = other.coeffs;
		return *this;
	}
	constexpr FwdAutoDiff& operator=(FwdAutoDiff &&other) {
		coeffs = std::move(other.coeffs);
		return *this;
	}
	constexpr FwdAutoDiff& operator=(double value) {
		coeffs.fill(0.0);
		coeffs[0] = value;
		return *this;
	}
	static constexpr FwdAutoDiff independent(double value, size_t dim) {
		FwdAutoDiff diff;
		index_type indices;
		indices.fill(0);
		indices[dim] = 1;
		diff.coeffs.fill(0.0);
		diff.coeffs[0] = value;
		diff.coeffs[flatten(indices)] = 1.0;
	}
	auto operator+() const {
		return *this;
	}
	auto operator-() const {
		auto A = *this;
		for (int i = 0; i < size(); i++) {
			A.coeffs[i] = -this->coeffs[i];
		}
		return A;
	}
	auto operator+(FwdAutoDiff A) const {
		for (int i = 0; i < size(); i++) {
			A.coeffs[i] += this->coeffs[i];
		}
		return A;
	}
	auto operator-(FwdAutoDiff A) const {
		return *this + (-A);
	}
	auto operator*(FwdAutoDiff const &B) const {
		auto const &A = *this;
		FwdAutoDiff C;
		C.coeffs.fill(0.0);
		for (auto α = zero();; increment(α)) {
			for (auto β = zero();; increment(β, α)) {
				auto const γ = add(α, β);
				C.coeffs[flatten(γ)] += A.coeffs[flatten(α)] * B.coeffs[flatten(β)];
			}
		}
		return C;
	}
	auto operator/(FwdAutoDiff const &B) const {
		return operator*(inverse(B));
	}
	auto operator+(double b) const {
		auto A = *this;
		for (int i = 0; i < size(); i++) {
			A.coeffs[i] += b;
		}
		return A;
	}
	auto operator-(double b) const {
		auto A = *this;
		return A + (-b);
	}
	auto operator*(double b) const {
		auto A = *this;
		for (int i = 0; i < size(); i++) {
			A.coeffs[i] *= b;
		}
		return A;
	}
	auto operator/(double b) const {
		return operator*(1.0 / b);
	}
	friend auto operator+(double b, FwdAutoDiff A) {
		return A + b;
	}
	friend auto operator-(double b, FwdAutoDiff A) {
		return -A + b;
	}
	friend auto operator*(double b, FwdAutoDiff A) {
		return A * b;
	}
	friend auto operator/(double b, FwdAutoDiff A) {
		return b * inverse(A);
	}
	friend constexpr auto exp(FwdAutoDiff const &f) {
		using std::exp;
		FwdAutoDiff g { };
		g.coeffs.fill(0.0);
		g.coeffs[0] = exp(f.coeffs[0]);
		auto α = zero();
		for (increment(α);; increment(α)) {
			auto const a = flatten(α);
			g.coeffs[a] = 0.0;
			auto β = zero();
			for (increment(β, α);; increment(β, α)) {
				auto const fβ = f.coeffs[flatten(β)];
				auto const gαβ = g.coeffs[flatten(sub(α, β))];
				g.coeffs[a] += fβ * gαβ / factorial(β);
			}
		}
		return g;
	}

private:
	std::array<double, size()> coeffs;
};
