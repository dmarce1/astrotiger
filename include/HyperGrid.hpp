#pragma once
#include "ContainerArithmetic.hpp"
#include "Hdf5.hpp"
#include "Matrix.hpp"
#include "dgTransforms.hpp"

#include <numeric>
#include <valarray>

template<typename T>
std::valarray<T> copysign(std::valarray<T> x, std::valarray<T> y) {
	int const N = x.size();
	assert(int(y.size()) == N);
	std::valarray<T> z(N);
	for (int i = 0; i < N; i++) {
		z[i] = std::copysign(x[i], y[i]);
	}
	return z;
}

template<typename T>
std::valarray<T> copysign(std::valarray<T> x, T const &y) {
	int const N = x.size();
	std::valarray<T> z(N);
	for (int i = 0; i < N; i++) {
		z[i] = std::copysign(x[i], y);
	}
	return z;
}

template<typename T>
std::valarray<T> copysign(T const &x, std::valarray<T> y) {
	int const N = y.size();
	std::valarray<T> z(N);
	for (int i = 0; i < N; i++) {
		z[i] = std::copysign(x, y[i]);
	}
	return z;
}

template<typename T>
std::valarray<T> min(std::valarray<T> x, std::valarray<T> y) {
	int const N = x.size();
	assert(int(y.size()) == N);
	std::valarray<T> z(N);
	for (int i = 0; i < N; i++) {
		z[i] = std::min(x[i], y[i]);
	}
	return z;
}

template<typename T>
std::valarray<T> min(std::valarray<T> x, T const &y) {
	int const N = x.size();
	std::valarray<T> z(N);
	for (int i = 0; i < N; i++) {
		z[i] = std::min(x[i], y);
	}
	return z;
}

template<typename T>
std::valarray<T> min(T const &x, std::valarray<T> y) {
	int const N = y.size();
	std::valarray<T> z(N);
	for (int i = 0; i < N; i++) {
		z[i] = std::min(x, y[i]);
	}
	return z;
}

template<typename T>
std::valarray<T> minmod(std::valarray<T> a, std::valarray<T> b) {
	std::valarray<T> const sgn = copysign(T(0.5), a) + copysign(T(0.5), b);
	a = std::abs(a);
	b = std::abs(b);
	std::valarray<T> const mag = min(a, b);
	return sgn * mag;
}

template<template<typename, int> typename S, typename T, int D, int N, int O, typename RK>
struct HyperGrid {
	static constexpr int nB = 2;
	static constexpr int nF = S<T, D>::fieldCount();
	static constexpr int nM = binco(O + D - 1, D);
	static constexpr int nN = ipow(O, D);
	static constexpr int nV = ipow(N + 2 * nB, D);
	static constexpr int nRK = RK::stageCount();
	static constexpr RK rk { };
	static constexpr Range<int, D> exteriorBox_ { repeat<D>(-nB), repeat<D>(N + nB) };
	static constexpr Range<int, D> interiorBox_ { repeat<D>(0), repeat<D>(N) };
	using ScalarState = S<T, D>;
	using InteriorIndex = MultiIndex<exteriorBox_, interiorBox_>;
	using ModeIndex = TriIndex<O, D>;
	HyperGrid(T const &xNint = T(1)) :
			Um_(nF * nM * nV), dx_(xNint / T(N)), dxInv_(T(N) / xNint) {
	}
	void initialize(std::function<ScalarState(std::array<T, D> const&)> const &initialState) {
		std::valarray<T> nodeState(nF * nN);
		std::valarray<T> thisModeState(nF * nM);
		std::valarray<size_t> sizes( { nF, nM });
		std::valarray<size_t> strides( { nM * nV, nV });
		for (auto iI = InteriorIndex::begin(); iI != InteriorIndex::end(); iI++) {
			int const ii = iI;
			for (int ni = 0; ni < nN; ni++) {
				auto const q = getQuadraturePoint<T, D, O>(ni);
				std::array<T, D> x;
				for (int d = 0; d < D; d++) {
					x[d] = (T(2 * iI[d] + 1) + q[d]) * dx_;
				}
				auto const thisState = initialState(x);
				for (int fi = 0; fi < nF; fi++) {
					nodeState[nN * fi + ni] = thisState[fi];
				}
			}
			for (int fi = 0; fi < nF; fi++) {
				dgAnalyze<T, D, O>(std::begin(nodeState) + nN * fi, std::begin(thisModeState) + nM * fi);
				dgMassInverse<T, D, O>(std::begin(thisModeState) + nM * fi, std::begin(thisModeState) + nM * fi);
			}
			Um_[std::gslice(ii, sizes, strides)] = thisModeState;
		}
	}
	void enforceBoundaryConditions() {
		std::valarray<size_t> sizes0(N + 2 * nB, D);
		std::valarray<size_t> strides0(D);
		int stride = 1;
		for (int d = 0; d < D; d++) {
			strides0[D - d - 1] = stride;
			stride *= N + 2 * nB;
		}
		for (int d = 0; d < D; d++) {
			auto sizes = sizes0;
			auto strides = strides0;
			sizes[d] = nB;
			int const from1 = strides0[d] * nB;
			int const to1 = from1 + N * strides0[d];
			int const to2 = 0;
			int const from2 = N * strides0[d];
			for (int fi = 0; fi < nF; fi++) {
				for (int mi = 0; mi < nM; mi++) {
					int const offset = (nM * fi + mi) * nV;
					std::valarray<T> const tmp1 = Um_[std::gslice(offset + from1, sizes, strides)];
					std::valarray<T> const tmp2 = Um_[std::gslice(offset + from2, sizes, strides)];
					Um_[std::gslice(offset + to1, sizes, strides)] = tmp1;
					Um_[std::gslice(offset + to2, sizes, strides)] = tmp2;
				}
			}
		}
//		for (auto iI = InteriorIndex::begin(); iI != InteriorIndex::end(); iI++) {
//			int const i = iI;
//			for (int mi = 0; mi < 1; mi++) {
//				for (int fi = 0; fi < 2; fi++) {
//					printf("(%i, %i) fi=%i mi=%i  i=%i %e \n", iI[0], iI[1], fi, mi, i, Um_[nV * (fi * nM + mi) + i]);
//				}
//			}
//		}
	}
	void applyLimiter() {
		std::valarray<size_t> sizes0(N + 2, D);
		std::valarray<size_t> strides0(D);
		int stride = 1;
		for (int d = 0; d < D; d++) {
			strides0[d] = stride;
			stride *= N + 2 * nB;
		}
		constexpr int nMlo = binco(O + D - 2, D);
		int const corner = strides0.sum();
		for (int fi = 0; fi < nF; fi++) {
			std::array<std::valarray<T>, D> Ulim;
			Ulim.fill(std::valarray<T>(nV * nMlo));
			for (int d = 0; d < D; d++) {
				for (int mi = 0; mi < nMlo; mi++) {
					auto const tri = flatToTriangular<D, O>(mi);
					int const p = tri[d];
					int const j0 = mi * nV + strides0[d];
					int const km = (nM * fi + mi) * nV + corner;
					int const k0 = km + strides0[d];
					int const kp = k0 + strides0[d];
					std::valarray<T> const u0 = Um_[std::gslice(k0, sizes0, strides0)];
					std::valarray<T> const up = Um_[std::gslice(kp, sizes0, strides0)];
					std::valarray<T> const um = Um_[std::gslice(km, sizes0, strides0)];
					T const clim = T(1) / T(2 * p + 1);
					Ulim[d][std::gslice(j0, sizes0, strides0)] = clim * minmod(std::valarray<T>(up - u0), std::valarray<T>(u0 - um));
				}
			}
		}
	}
	T beginStep() {
		applyLimiter();
		std::valarray<size_t> sizes( { nM });
		std::valarray<size_t> strides( { nV });

		Um0_ = Um_;
		dUm_.fill(std::valarray<T>(T(0), Um_.size()));
		std::valarray<T> lambdaMax(T(0), D);
		std::valarray<T> nodeState(nF * nN);
		std::valarray<T> thisModeState(nF * nM);
		for (auto iI = InteriorIndex::begin(); iI != InteriorIndex::end(); iI++) {
			int const ii = iI;
			for (int fi = 0; fi < nF; fi++) {
				std::valarray<T> const thisModeState = Um_[std::gslice(ii + fi * nM * nV, sizes, strides)];
				dgSynthesize<T, D, O>(std::begin(thisModeState), std::begin(nodeState) + fi * nN);
			}
			std::valarray<T> const tmp = nodeState;
			nodeState[std::gslice(0, { nF, nN }, { 1, nF })] = tmp;
			for (int ni = 0; ni < nN; ni++) {
				ScalarState thisState;
				std::copy_n(std::begin(nodeState) + ni * nF, nF, thisState.begin());
				for (int d = 0; d < D; d++) {
					auto const eigenvalues = thisState.eigenvalues(d);
					for (auto lambda : eigenvalues) {
						lambdaMax[d] = std::max(lambdaMax[d], abs(lambda));
					}
				}
			}

		}
		T lambdaSum = T(0);
		for (int d = 0; d < D; d++) {
			lambdaSum += lambdaMax[d];
		}
		T const dt = (dx_ * rk.cfl()) / (T(2 * O - 1) * lambdaSum);
		return dt;
	}
private:
	std::valarray<T> Um_;
	std::valarray<T> Um0_;
	std::array<std::valarray<T>, nRK> dUm_;
	T dx_;
	T dxInv_;
};

//class HyperGrid {
//
//	static constexpr int dimensionCount = State::dimCount();
//	static constexpr int ghostWidth = 2;
//	static constexpr int fieldCount = State::fieldCount();
//	static constexpr int nRK = RungeKutta::stageCount();
//	static constexpr int cellsAcrossExterior = cellsAcrossInterior + 2 * ghostWidth;
//	static constexpr int exteriorVolume = ipow(cellsAcrossExterior, dimensionCount);
//	static constexpr RungeKutta butcherTable { };
//	static constexpr int modeVolume = binco(basisOrder + dimensionCount - 1, dimensionCount);
//	static constexpr int nodeVolume = ipow(basisOrder, dimensionCount);
//	static constexpr Range<int, dimensionCount> exteriorBox { repeat<dimensionCount>(-ghostWidth), repeat<dimensionCount>(cellsAcrossInterior + ghostWidth) };
//	static constexpr Range<int, dimensionCount> interiorBox { repeat<dimensionCount>(0), repeat<dimensionCount>(cellsAcrossInterior) };
//
//	using Type = State::value_type;
//	using InteriorIndex = MultiIndex<exteriorBox, interiorBox>;
//	using ModeIndex = TriIndex<basisOrder,dimensionCount>;
//
//	std::vector<std::array<std::array<std::array<Type, exteriorVolume>, modeVolume>, fieldCount>> stageDerivatives_;
//	std::vector<std::array<std::array<Type, exteriorVolume>, modeVolume>> nextState;
//	std::vector<std::array<std::array<Type, exteriorVolume>, modeVolume>> currentState;
//	const Type cellWidth;
//	const Type inverseCellWidth;
//
//	static constexpr auto interiorIndexMap(int i) {
//		static constexpr auto map = createMultiIndexMap<exteriorBox, interiorBox>();
//		return map[i];
//	}
//	static constexpr int stride(int d) {
//		return ipow(cellsAcrossExterior, dimensionCount - 1 - d);
//	}
//
//public:
//	HyperGrid(Type const &xNint = Type(1)) :
//			nextState(fieldCount), cellWidth(xNint / Type(cellsAcrossInterior)), inverseCellWidth(Type(cellsAcrossInterior) / xNint) {
//	}
//	void initialize(std::function<State(std::array<Type, dimensionCount> const&)> const &initialState) {
//		Type const halfCellWidth = Type(0.5) * cellWidth;
//		for (auto cellMultiIndex = InteriorIndex::begin(); cellMultiIndex != InteriorIndex::end(); cellMultiIndex++) {
//			int const cellFlatIndex = cellMultiIndex;
//			std::array<std::array<Type, nodeVolume>, fieldCount> nodeState;
//			for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
//				auto const quadraturePoint = getQuadraturePoint<Type, dimensionCount, basisOrder>(nodeIndex);
//				std::array<Type, dimensionCount> position;
//				for (int dimension = 0; dimension < dimensionCount; dimension++) {
//					position[dimension] = (Type(2 * cellMultiIndex[dimension] + 1) + quadraturePoint[dimension]) * halfCellWidth;
//				}
//				auto const thisState = initialState(position);
//				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//					nodeState[fieldIndex][nodeIndex] = thisState[fieldIndex];
//				}
//			}
//			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//				auto modeState = dgAnalyze<Type, dimensionCount, basisOrder>(nodeState[fieldIndex]);
//				modeState = dgMassInverse<Type, dimensionCount, basisOrder>(modeState);
//				for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
//					nextState[fieldIndex][modeIndex][cellFlatIndex] = modeState[modeIndex];
//				}
//			}
//		}
//	}
//	void output(const char *filenameBase, int timeStepNumber, Type const &time) {
//		std::string filename = std::string(filenameBase) + "." + std::to_string(timeStepNumber) + ".h5";
//		writeHdf5<Type, dimensionCount, cellsAcrossInterior, basisOrder, ghostWidth>(filename, cellWidth, nextState, State::getFieldNames());
//		writeList("X.visit", "!NBLOCKS 1\n", filename + ".xmf");
//	}
//	void enforceBoundaryConditions() {
//		using Index = MultiIndex<exteriorBox>;
//		for (auto ghostZoneIndex = Index::begin(); ghostZoneIndex != Index::end(); ghostZoneIndex++) {
//			Index interiorIndex;
//			bool isGhostZone = false;
//			for (int dimensionIndex = 0; dimensionIndex < dimensionCount; dimensionIndex++) {
//				if (ghostZoneIndex[dimensionIndex] < 0) {
//					isGhostZone = true;
//					interiorIndex[dimensionIndex] = ghostZoneIndex[dimensionIndex] + cellsAcrossInterior;
//				} else if (ghostZoneIndex[dimensionIndex] >= cellsAcrossInterior) {
//					isGhostZone = true;
//					interiorIndex[dimensionIndex] = ghostZoneIndex[dimensionIndex] - cellsAcrossInterior;
//				} else
//					interiorIndex[dimensionIndex] = ghostZoneIndex[dimensionIndex];
//			}
//			if (!isGhostZone) {
//				continue;
//			}
//			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//				for (int basisIndex = 0; basisIndex < modeVolume; basisIndex++) {
//					nextState[fieldIndex][basisIndex][ghostZoneIndex] = nextState[fieldIndex][basisIndex][interiorIndex];
//				}
//			}
//		}
//	}
//	Type beginStep() {
//		applyLimiter();
//		using namespace Math;
//		currentState = nextState;
//		stageDerivatives_.resize(nRK);
//		std::array<Type, dimensionCount> maximumEigenvalue;
//		maximumEigenvalue.fill(Type(0));
//		for (auto cellMultiIndex = InteriorIndex::begin(); cellMultiIndex != InteriorIndex::end(); cellMultiIndex++) {
//			int const cellFlatIndex = cellMultiIndex;
//			std::array<std::array<Type, modeVolume>, fieldCount> modeState;
//			std::array<std::array<Type, nodeVolume>, fieldCount> nodeState;
//			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//				for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
//					modeState[fieldIndex][modeIndex] = nextState[fieldIndex][modeIndex][cellFlatIndex];
//				}
//				nodeState[fieldIndex] = dgSynthesize<Type, dimensionCount, basisOrder>(modeState[fieldIndex]);
//			}
//			State thisState;
//			for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
//				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//					thisState[fieldIndex] = nodeState[fieldIndex][nodeIndex];
//				}
//			}
//			for (int dimension = 0; dimension < dimensionCount; dimension++) {
//				auto const eigenvalues = thisState.eigenvalues(dimension);
//				for (auto thisEigenvalue : eigenvalues) {
//					maximumEigenvalue[dimension] = max(maximumEigenvalue[dimension], abs(thisEigenvalue));
//				}
//			}
//		}
//		Type maximumEigenvalueSum = Type(0);
//		for (int dimension = 0; dimension < dimensionCount; dimension++) {
//			maximumEigenvalueSum += maximumEigenvalue[dimension];
//		}
//		Type const timeStepSize = (cellWidth * butcherTable.cfl()) / (Type(2 * basisOrder - 1) * maximumEigenvalueSum);
//		return timeStepSize;
//	}
//	void subStep(Type const &timeStepSize, int stageIndex) {
//		nextState = currentState;
//		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//			for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
//				for (int thisStage = 0; thisStage < stageIndex; thisStage++) {
//					for (auto cellMultiIndex = InteriorIndex::begin(); cellMultiIndex != InteriorIndex::end(); cellMultiIndex++) {
//						int const cellFlatIndex = cellMultiIndex;
//						nextState[fieldIndex][modeIndex][cellFlatIndex] +=
//								butcherTable.a(stageIndex, thisStage) * stageDerivatives_[thisStage][fieldIndex][modeIndex][cellFlatIndex];
//					}
//				}
//				stageDerivatives_[stageIndex][fieldIndex][modeIndex].fill(Type(0));
//			}
//		}
//		enforceBoundaryConditions();
//		applyLimiter();
//		computeDudt(timeStepSize, stageDerivatives_[stageIndex], std::make_integer_sequence<int, dimensionCount> { });
//	}
//	void endStep() {
//		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//			for (auto modeMultiIndex = ModeIndex::begin(); modeMultiIndex != ModeIndex::end(); modeMultiIndex++) {
//				int const modeFlatIndex = modeMultiIndex;
//				for (auto cellMultiIndex = InteriorIndex::begin(); cellMultiIndex != InteriorIndex::end(); cellMultiIndex++) {
//					int const cellFlatIndex = cellMultiIndex;
//					nextState[fieldIndex][modeFlatIndex][cellFlatIndex] = currentState[fieldIndex][modeFlatIndex][cellFlatIndex];
//					for (int stageIndex = 0; stageIndex < nRK; stageIndex++) {
//						nextState[fieldIndex][modeFlatIndex][cellFlatIndex] +=
//								butcherTable.b(stageIndex) * stageDerivatives_[stageIndex][fieldIndex][modeFlatIndex][cellFlatIndex];
//					}
//				}
//			}
//		}
//		stageDerivatives_ = { };
//		currentState = { };
//	}
//	void applyLimiter() {
//		constexpr auto limiterCoefficients = []() {
//			std::array<Type, basisOrder> limiterCoefficients;
//			for (int orderIndex = 0; orderIndex < basisOrder; orderIndex++) {
//				limiterCoefficients[orderIndex] = Type(1) / Type(2 * orderIndex + 1);
//			}
//			return limiterCoefficients;
//		}();
//
//		constexpr Range<int, dimensionCount> limiterBox { repeat<dimensionCount>(-1), repeat<dimensionCount>(cellsAcrossInterior + 1) };
//
//		using LimiterIndex = MultiIndex<exteriorBox, limiterBox>;
//
//		for (int polynomialDegree = basisOrder - 1; polynomialDegree > 0; polynomialDegree--) {
//			for (auto targetModeIndex = ModeIndex::begin(); targetModeIndex != ModeIndex::end(); targetModeIndex++) {
//				for (int dimension = 0; dimension < dimensionCount; dimension++) {
//
//					int totalDegree = 0;
//					for (int dimension = 0; dimension < dimensionCount; dimension++) {
//						totalDegree += targetModeIndex[dimension];
//						if (totalDegree > polynomialDegree) {
//							break;
//						}
//					}
//					if ((totalDegree != polynomialDegree) || (targetModeIndex[dimension] == 0)) {
//						continue;
//					}
//
//					auto lowerModeIndex = targetModeIndex.dec(dimension);
//					auto const lowerModeFlatIndex = lowerModeIndex;
//					auto const targetModeFlatIndex = targetModeIndex;
//					auto const strideOffset = stride(dimension);
//
//					for (auto cellMultiIndex = LimiterIndex::begin(); cellMultiIndex != LimiterIndex::end(); cellMultiIndex++) {
//						int const cellFlatIndex = cellMultiIndex;
//						State highOrderMode, referenceState, slopeDifference;
//
//						for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//							highOrderMode[fieldIndex] = nextState[fieldIndex][targetModeFlatIndex][cellFlatIndex];
//							referenceState[fieldIndex] = nextState[fieldIndex][0][cellFlatIndex];
//							auto const lowOrderMode = nextState[fieldIndex][lowerModeFlatIndex][cellFlatIndex];
//							auto const differenceRight = nextState[fieldIndex][lowerModeFlatIndex][cellFlatIndex + strideOffset] - lowOrderMode;
//							auto const differenceLeft = lowOrderMode - nextState[fieldIndex][lowerModeFlatIndex][cellFlatIndex - strideOffset];
//							slopeDifference[fieldIndex] = minmod(differenceRight, differenceLeft);
//						}
//
//						auto const [_, rightEigenvectors] = referenceState.eigenSystem(dimension);
//						auto const leftEigenvectors = matrixInverse(rightEigenvectors);
//						auto const limiterCoefficient = limiterCoefficients[lowerModeIndex[dimension]];
//
//						highOrderMode = leftEigenvectors * highOrderMode;
//						slopeDifference = leftEigenvectors * slopeDifference;
//
//						for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//							highOrderMode[fieldIndex] = minmod(highOrderMode[fieldIndex], limiterCoefficient * slopeDifference[fieldIndex]);
//						}
//
//						highOrderMode = rightEigenvectors * highOrderMode;
//
//						for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//							nextState[fieldIndex][targetModeFlatIndex][cellFlatIndex] = highOrderMode[fieldIndex];
//						}
//
//					}
//				}
//			}
//		}
//	}
//
//	template<int ... dimension>
//	void computeDudt(Type timeStepSize, std::array<std::array<std::array<Type, exteriorVolume>, modeVolume>, fieldCount> &stateDerivative,
//			std::integer_sequence<int, dimension...>) {
//		(computeDudtByDim<dimension>(timeStepSize, stateDerivative), ...);
//	}
//	template<int dimension>
//	void computeDudtByDim(Type timeStepSize, std::array<std::array<std::array<Type, exteriorVolume>, modeVolume>, fieldCount> &stateDerivative) {
//		constexpr Range<int, dimensionCount> interiorBoxPlusOne { repeat<dimensionCount>(0), repeat<dimensionCount>(cellsAcrossInterior) + unit<dimensionCount>(
//				dimension) };
//		using FluxIndexType = MultiIndex<exteriorBox, interiorBoxPlusOne>;
//		Type const lambda = Type(2) * timeStepSize * inverseCellWidth;
//		for (auto cellMultIndex = FluxIndexType::begin(); cellMultIndex != FluxIndexType::end(); cellMultIndex++) {
//			int const rightInterfaceIndex = cellMultIndex;
//			int const leftInterfaceIndex = rightInterfaceIndex - stride(dimension);
//
//			std::array<std::array<Type, triangleSize<dimensionCount, basisOrder>>, fieldCount> volumeModesLeft, volumeModesRight;
//			std::array<std::array<Type, squareSize<dimensionCount - 1, basisOrder>>, fieldCount> surfaceNodesLeft, surfaceNodesRight, surfaceNodesFlux;
//			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//				for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
//					volumeModesLeft[fieldIndex][modeIndex] = nextState[fieldIndex][modeIndex][leftInterfaceIndex];
//					volumeModesRight[fieldIndex][modeIndex] = nextState[fieldIndex][modeIndex][rightInterfaceIndex];
//				}
//				surfaceNodesLeft[fieldIndex] = dgSynthesize<Type, dimensionCount - 1, basisOrder>(
//						dgTrace<Type, dimensionCount, basisOrder>(2 * dimension + 1, volumeModesLeft[fieldIndex]));
//				surfaceNodesRight[fieldIndex] = dgSynthesize<Type, dimensionCount - 1, basisOrder>(
//						dgTrace<Type, dimensionCount, basisOrder>(2 * dimension + 0, volumeModesRight[fieldIndex]));
//			}
//
//			for (int nodeIndex = 0; nodeIndex<squareSize < dimensionCount - 1, basisOrder> ; nodeIndex++) {
//				State stateRight, stateLeft;
//				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//					stateLeft[fieldIndex] = surfaceNodesLeft[fieldIndex][nodeIndex];
//					stateRight[fieldIndex] = surfaceNodesRight[fieldIndex][nodeIndex];
//				}
//				auto const riemannFlux = solveRiemannProblem(stateLeft, stateRight, dimension);
//				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//					surfaceNodesFlux[fieldIndex][nodeIndex] = riemannFlux[fieldIndex];
//				}
//			}
//			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//				auto const modes = dgAnalyze<Type, dimensionCount - 1, basisOrder>(surfaceNodesFlux[fieldIndex]);
//				volumeModesLeft[fieldIndex] = dgTraceInverse<Type, dimensionCount, basisOrder>(2 * dimension + 1, modes);
//				volumeModesLeft[fieldIndex] = dgMassInverse<Type, dimensionCount, basisOrder>(volumeModesLeft[fieldIndex]);
//				volumeModesRight[fieldIndex] = dgTraceInverse<Type, dimensionCount, basisOrder>(2 * dimension + 0, modes);
//				volumeModesRight[fieldIndex] = dgMassInverse<Type, dimensionCount, basisOrder>(volumeModesRight[fieldIndex]);
//			}
//
//			if (cellMultIndex[dimension] > 0) {
//				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//					for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
//						stateDerivative[fieldIndex][modeIndex][leftInterfaceIndex] -= lambda * volumeModesLeft[fieldIndex][modeIndex];
//					}
//				}
//			}
//			if (cellMultIndex[dimension] < cellsAcrossInterior) {
//				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//					for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
//						stateDerivative[fieldIndex][modeIndex][rightInterfaceIndex] += lambda * volumeModesRight[fieldIndex][modeIndex];
//					}
//				}
//			}
//		}
//		for (auto cellMultiIndex = InteriorIndex::begin(); cellMultiIndex != InteriorIndex::end(); cellMultiIndex++) {
//			int const cellFlatIndex = cellMultiIndex;
//			std::array<State, nodeVolume> nodeStates;
//			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//				std::array<Type, modeVolume> modalState;
//				for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
//					modalState[modeIndex] = nextState[fieldIndex][modeIndex][cellFlatIndex];
//				}
//				auto thisValue = dgSynthesize<Type, dimensionCount, basisOrder>(modalState);
//				for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
//					nodeStates[nodeIndex][fieldIndex] = thisValue[nodeIndex];
//				}
//			}
//			for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
//				nodeStates[nodeIndex] = nodeStates[nodeIndex].flux(dimension);
//			}
//			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//				std::array<Type, nodeVolume> nodes;
//				for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
//					nodes[nodeIndex] = nodeStates[nodeIndex][fieldIndex];
//				}
//				auto flux = dgAnalyze<Type, dimensionCount, basisOrder>(nodes);
//				flux = dgMassInverse<Type, dimensionCount, basisOrder>(flux);
//				flux = dgStiffness<Type, dimensionCount, basisOrder>(dimension, flux);
//				flux = dgMassInverse<Type, dimensionCount, basisOrder>(flux);
//				for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
//					stateDerivative[fieldIndex][modeIndex][cellFlatIndex] += lambda * flux[modeIndex];
//				}
//			}
//		}
//	}
//
//}
//;
//
