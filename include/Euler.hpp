#pragma once

#include "Definitions.hpp"
#include "SquareMatrix.hpp"

#include <array>
#include <cmath>
#include <vector>

#define EULERS_CONSTRUCTION : \
	rho(*std::launder(reinterpret_cast<T*>(BaseClass::data() + 0))), \
	eg (*std::launder(reinterpret_cast<T*>(BaseClass::data() + 1))), \
	S  (*std::launder(reinterpret_cast<std::array<T, Parameters::dimensionCount>*>(BaseClass::data() + 2)))

template<typename T>
struct EulerState: public std::array<T, 2 + Parameters::dimensionCount> {
	static constexpr int fieldCount_ = 2 + Parameters::dimensionCount;
	static constexpr Real gamma_ = Real(5) / Real(3);
	using BaseClass = std::array<T, fieldCount_>;
	using Eigensystem = std::pair<std::array<T, fieldCount_>, SquareMatrix<T, fieldCount_>>;
	static int fieldCount() noexcept {
		return fieldCount_;
	}
	EulerState() EULERS_CONSTRUCTION {
	}
	EulerState(BaseClass const &other) EULERS_CONSTRUCTION {
		((BaseClass&)(*this)).operator=(other);
	}
	EulerState(EulerState const &other) EULERS_CONSTRUCTION {
		*this = other;
	}
	EulerState(T const &other) EULERS_CONSTRUCTION {
		for( auto& u : *this) {
			u = other;
		}
	}
	EulerState(EulerState &&other) EULERS_CONSTRUCTION {
		*this = std::move(other);
	}
	EulerState& operator=(EulerState const &other) {
		static_cast<BaseClass&>(*this) = static_cast<BaseClass const&>(other);
		return *this;
	}
	EulerState& operator=(EulerState &&other) {
		static_cast<BaseClass&>(*this) = std::move(static_cast<BaseClass&&>(other));
		return *this;
	}
	std::array<T, fieldCount_> eigenvalues(int dim) const {
		std::array<T, fieldCount_> lambda;
		using std::sqrt;
		T const irho = T(1) / rho;
		T const v = S[dim] * irho;
		T const ek = T(0.5) * irho * dot(S, S);
		T const ei = std::max(eg - ek, T(0));
		T const p = (gamma - T(1)) * ei;
		T const a = sqrt(gamma * p * irho);
		std::fill(lambda.begin(), lambda.end(), v);
		lambda.front() -= a;
		lambda.back() += a;
		return lambda;

	}
	Eigensystem eigenSystem(int dim) const {
		Eigensystem rc;
		using namespace Parameters;
		using std::sqrt;
		auto& lambda = rc.first;
		auto& R = rc.second;
		std::fill(R.begin(), R.end(), T(0));
		T const invρ = T(1) / rho;
		auto const v = S * invρ;
		T ke = T(0);
		for(int d = 0; d < dimensionCount; d++)
		ke += T(0.5) * v[d] * S[d];
		T const ei = std::max(T(0), eg - ke);
		T const p = (gamma - T(1)) * ei;
		T const c = sqrt(gamma * p * invρ);
		T const h = (eg + p) * invρ;
		T const u_n = v[dim];
		for(int i = 0; i < fieldCount_; i++) {
			lambda[i] = u_n;
		}
		lambda[0] = u_n - c;
		lambda[fieldCount_ - 1] = u_n + c;
		int col = 0;
		R(0, col) = T(1);
		for(int i = 0; i < dimensionCount; ++i) {
			R(1 + i, col) = (i == dim ? u_n - c : v[i]);
		}
		R(dimensionCount + 1, col) = h - u_n * c;
		col++;
		for(int d = 0; d < dimensionCount; ++d) {
			if (d == dim) continue;
			R(0, col) = T(0);
			R(1 + d, col) = T(1);
			R(dimensionCount + 1, col) = T(0);
			col++;
		}
		R(0, col) = T(1);
		for(int i = 0; i < dimensionCount; i++) {
			R(1 + i, col) = v[i];
		}
		R(dimensionCount + 1, col) = T(0);
		col++;
		R(0, col) = T(1);
		for(int i = 0; i < dimensionCount; ++i) {
			R(1 + i, col) = (i == dim ? u_n + c : v[i]);
		}
		R(dimensionCount + 1, col) = h + u_n * c;
		return rc;
	}
	EulerState flux(int d) const noexcept {
		EulerState F;
		T const irho = T(1) / rho;
		auto const v = S * irho;
		T const ek = T(0.5) * dot(v, S);
		T const ei = std::max(T(0), eg - ek);
		T const p = (gamma - T(1)) * ei;
		T const u = v[d];
		F.rho = S[d];
		F.eg = u * (eg + p);
		F.S = u * S;
		F.S[d] += p;
		return F;
	}
	friend EulerState solveRiemannProblem(const EulerState &uL, const EulerState &uR, int dim) noexcept {
		using namespace Parameters;
		using std::sqrt;
		using State = EulerState;
		constexpr int iN2 = 2;
		constexpr int iN3 = 3;
		constexpr int iL = 0;
		constexpr int iR = 1;
		constexpr int iStar = 2;
		std::array<State, iN3> u = {uL, uR};
		std::array<State, iN2> f;
		std::array<T, iN2> irho, v, ek, ei, a;
		std::array<T, iN3> s, p;
		for(int i = 0; i < iN2; i++) {
			irho[i] = T(1) / u[i].rho;
			v[i] = irho[i] * u[i].S[dim];
			ek[i] = T(0);
			for(int d = 0; d < dimensionCount; d++) {
				ek[i] += T(0.5) * irho[i] * sqr(u[i].S[d]);
			}
			ei[i] = std::max(T(0), u[i].eg - ek[i]);
			p[i] = (gamma - T(1)) * ei[i];
			a[i] = sqrt(gamma * p[i] * irho[i]);
		}
		s[iL] = std::min(v[iL] - a[iL], v[iR] - a[iR]);
		s[iR] = std::max(v[iL] + a[iL], v[iR] + a[iR]);
		T const num = p[iR] - p[iL] + u[iL].rho * v[iL] * (s[iL] - v[iL]) - u[iR].rho * v[iR] * (s[iR] - v[iR]);
		T const den = u[iL].rho * (s[iL] - v[iL]) - u[iR].rho * (s[iR] - v[iR]);
		s[iStar] = num / den;
		for(int i = 0; i < iN2; i++) {
			f[i] = u[i].flux(dim);
		}
		if (T(0) < s[iL]) {
			return f[iL];
		} else if (T(0) > s[iR]) {
			return f[iR];
		} else {
			int const i = ((s[iStar] > T(0)) ? iL : iR);
			u[iStar].rho = u[i].rho * (s[i] - v[i]) / (s[i] - s[iStar]);
			for(int dir = 0; dir < dimensionCount; dir++) {
				if(dim != dir ) {
					u[iStar].S[dir] = u[iStar].rho * irho[i] * u[i].S[dir];
				} else {
					u[iStar].S[dir] = u[iStar].rho * s[iStar];
				}
			}
			p[iStar] = p[i] + u[i].rho * (s[i] - v[i]) * (s[iStar] - v[i]);
			u[iStar].eg = ((s[i] - v[i]) * u[i].eg - p[i] * v[i] + p[iStar] * s[iStar]) / (s[i] - s[iStar]);
			return f[i] + s[i] * (u[iStar] - u[i]);
		}
	}
	T const& getDensity() const {
		return rho;
	}
	T const& getEnergy() const {
		return eg;
	}
	std::array<T, Parameters::dimensionCount> const& getMomentum() const {
		return S;
	}
	T const& getMomentum(int d) const {
		return S[d];
	}
	void setDensity(T const& value) {
		rho = value;
	}
	void setEnergy(T const& value) {
		eg = value;
	}
	void setMomentum(std::array<T, Parameters::dimensionCount> const& value) {
		S = value;
	}
	void setMomentum(int d, T const& value) {
		S[d] = value;
	}
	static std::vector<std::string> getFieldNames() {
		using namespace Parameters;
		static std::vector<std::string> const fieldNames = []() {
			std::vector<std::string> names;
			names.push_back("rho");
			names.push_back("eg");
			for(int d = 0; d < dimensionCount; d++) {
				std::string const name = std::string("s_") + std::string(1, 'x' + d);
				names.push_back(name);
			}
			return names;
		}();
		return fieldNames;
	}
private:
	T& rho;
	T& eg;
	std::array<T, Parameters::dimensionCount>& S;
};

