#include "Definitions.hpp"
#include "Indent.hpp"
#include "Util.hpp"

#include <algorithm>
#include <iostream>
#include <functional>
#include <numeric>
#include <map>
#include <vector>
#include <utility>
constexpr auto tiny = std::numeric_limits<double>::epsilon();

template<typename T>
std::vector<int> sortWithPermutation(std::vector<T> &vec, std::function<bool(T const&, T const&)> const &less) {
	std::vector<std::pair<T, int>> paired;
	paired.reserve(vec.size());
	for (int i = 0; i < static_cast<int>(vec.size()); ++i) {
		paired.emplace_back(vec[i], i);
	}
	std::stable_sort(paired.begin(), paired.end(), [less](auto const &a, auto const &b) {
		return less(a.first, b.first);
	});
	std::vector<int> permutation(paired.size());
	for (int i = 0; i < int(permutation.size()); i++) {
		vec[i] = std::move(paired[i].first);
		permutation[i] = paired[i].second;
	}
	return permutation;
}

using Real = long double;

static Indent indent { };

struct Constants {
	Constants(std::string const &name = "constant") :
			name_(name) {
		setType("double");
	}
	std::string operator()(Real constant) {
		auto iterator = indexes_.find(constant);
		if (iterator == indexes_.end()) {
			indexes_.insert(std::make_pair(constant, nextIndex_));
			iterator = indexes_.find(constant);
			constants_.push_back(constant);
			nextIndex_++;
		}
		return name_ + "[" + std::to_string(iterator->second) + std::string("]");
	}
	std::string getCode() const {
		std::ostringstream code;
		std::ostringstream line;
		code << "std::array<Type, " << indexes_.size() << "> " << name_ << " {";
		indent++;
		bool first = true;
		code << std::string(indent);
		int counter = 0;
		for (int i = 0; i < int(constants_.size()); i++) {
			if (!first) {
				line << ", ";
			}
			if (counter++ % 3 == 0) {
				code << line.str();
				line = std::ostringstream { };
				line << "\n" << std::string(indent);
			}
			line << indent << "Type {" << std::setprecision(std::numeric_limits<double>::max_digits10 - 1) << std::scientific << constants_[i] << "}";
			first = false;
		}
		code << line.str();
		line.clear();
		indent--;
		return code.str();
	}
	void setType(std::string type) {
		type_ = type;
	}
private:
	std::unordered_map<double, int> indexes_;
	std::vector<double> constants_;
	std::string type_;
	std::string name_ = "C";
	int nextIndex_ = 0;
};

static Constants getConstant;
using Matrix = std::vector<std::vector<Real>>;

Matrix createMatrix(int N, int M = -1) {
	if (M < 0) {
		M = N;
	}
	return std::vector<std::vector<Real>>(N, std::vector<Real>(M, 0.0));
}

Matrix identityMatrix(int N) {
	auto I = createMatrix(N);
	for (int n = 0; n < N; n++) {
		I[n][n] = Real(1);
	}
	return I;
}

Matrix matrixTranspose(Matrix const &A) {
	int rows = A.size(), cols = A[0].size();
	Matrix AT(cols, std::vector<Real>(rows));
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			AT[j][i] = A[i][j];
		}
	}
	return AT;
}

Matrix matrixMultiply(Matrix const &A, Matrix const &B) {
	int rows = A.size(), cols = B[0].size(), inner = A[0].size();
	assert(A[0].size() == B.size());
	Matrix C(rows, std::vector<Real>(cols, 0.0));
	for (int i = 0; i < rows; ++i) {
		for (int k = 0; k < inner; ++k) {
			for (int j = 0; j < cols; ++j) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	return C;
}

Matrix matrixInverse(Matrix const &A) {
	int const N = A.size();
	for (auto const &row : A) {
		if (int(row.size()) != N) {
			throw std::runtime_error("Matrix is not square");
		}
	}
	Matrix I = identityMatrix(N);
	Matrix B = A;
	for (int i = 0; i < N; ++i) {
		//	std::cout << matrixToString(A) << "\n";
		int pivot = i;
		if (B[i][i] < std::sqrt(std::numeric_limits<double>::epsilon())) {
			int maxB = B[i][i];
			for (int row = i + 1; row < N; ++row) {
				if (std::abs(B[row][i]) > std::abs(maxB)) {
					maxB = B[row][i];
					pivot = row;
				}
			}
		}
		if (std::abs(B[pivot][i]) < std::sqrt(std::numeric_limits<Real>::epsilon())) {
			throw std::runtime_error("Matrix is singular or nearly singular");
		}
		if (pivot != i) {
			std::swap(B[i], B[pivot]);
			std::swap(I[i], I[pivot]);
		}
		Real diag = B[i][i];
		for (int j = 0; j < N; ++j) {
			B[i][j] /= diag;
			I[i][j] /= diag;
		}
		for (int k = 0; k < N; ++k) {
			if (k == i) {
				continue;
			}
			Real factor = B[k][i];
			for (int j = 0; j < N; ++j) {
				B[k][j] -= factor * B[i][j];
				I[k][j] -= factor * I[i][j];
			}
		}
	}
	return I;
}

Matrix kroneckerProduct(Matrix const &A, Matrix const &B) {
	int const NR = A.size();
	int const NC = A[0].size();
	int const MR = B.size();
	int const MC = B[0].size();
	int const LR = NR * MR;
	int const LC = NC * MC;
	auto C = createMatrix(LR, LC);
	for (int nr = 0; nr < NR; nr++) {
		for (int nc = 0; nc < NC; nc++) {
			for (int mr = 0; mr < MR; mr++) {
				for (int mc = 0; mc < MC; mc++) {
					C[nr * MR + mr][nc * MC + mc] = A[nr][nc] * B[mr][mc];
				}
			}
		}
	}
	return C;
}

std::string matrixToString(Matrix const &A) {
	char *ptr;
	std::string str;
	for (int n = 0; n < (int) A.size(); n++) {
		for (int m = 0; m < (int) A[n].size(); m++) {
			Real value = std::abs(A[n][m]) < tiny ? Real(0) : A[n][m];
			asprintf(&ptr, "%4.1f ", (double) value);
			str += ptr;
			free(ptr);
		}
		str += "\n";
	}
	str += "\n";
	return str;
}

enum class Quadrature : int {
	gaussLegendre, gaussLobatto
};

enum class TransformDirection : int {
	forward, backward
};

struct QuadraturePoint {
	Real position;
	Real weight;
};

int triangular2flat(std::vector<int> const &I) {
	int const D = I.size();
	int flat = 0;
	int sum = 0;
	for (int d1 = D - 1; d1 >= 0; d1--) {
		sum += I[d1];
		int num = 1;
		int den = 1;
		for (int d2 = 0; d2 < D - d1; d2++) {
			num *= sum + d2;
			den *= d2 + 1;
		}
		flat += num / den;
	}
	return flat;
}

std::vector<int> flat2triangular(int D, int flat) {
	std::vector<int> I(D);
	I[0] = std::floor(std::pow(factorial<int>(D) * flat, Real(1) / Real(D)));
	for (int dim = 0; dim < D; dim++) {
		int count = binco(I[dim] + D - dim - 1, D - dim);
		while (count > flat) {
			assert(I[dim]);
			I[dim]--;
			count = binco(I[dim] + D - dim - 1, D - dim);
		}
		flat -= count;
		if (dim + 1 < D) {
			I[dim + 1] = I[dim];
		}
	}
	for (int d = 0; d < D - 1; d++) {
		I[d] -= I[d + 1];
	}
	return I;
}

Real legendreP(int n, Real eta, int m = 0) {
	Real value;
	if (std::abs(eta) != Real(1)) {
		value = std::sqrt(Real(1) / ipow(Real(1) - eta * eta, m)) * std::assoc_legendre(n, m, eta);
	} else {
		value = Real(1);
		for (int k = 1; k < m; k++) {
			value *= Real(n + k) / Real(2 * k * (n - k));
		}
		value *= Real((eta < Real(0)) ? nonepow(n + m) : +1);
	}
	return value;
}

std::vector<QuadraturePoint> gaussQuadrature(int nodeCount, Quadrature quadrature) {
	constexpr Real pi = std::numbers::pi_v<Real>;
	std::vector<QuadraturePoint> results;
	switch (quadrature) {
	case Quadrature::gaussLegendre: {
		for (int pointIndex = 0; pointIndex < nodeCount; pointIndex++) {
			QuadraturePoint point;
			Real theta, rootEquationDerivative, rootEquation;
			Real newTheta = pi * (Real(1) - Real(2 * pointIndex + 1) / Real(2) / Real(nodeCount));
			do {
				theta = newTheta;
				point.position = std::cos(theta);
				rootEquation = legendreP(nodeCount, point.position);
				rootEquationDerivative = legendreP(nodeCount, point.position, 1) * std::sin(theta);
				newTheta = theta + rootEquation / rootEquationDerivative;
			} while (double(newTheta) != std::nextafter(double(theta), double(newTheta)));
			point.weight = Real(2) / (sqr(legendreP(nodeCount, point.position, 1)) * (Real(1) - sqr(point.position)));
			results.push_back(point);
		}
		break;
	}
	case Quadrature::gaussLobatto: {
		Real const baseWeight = Real(2) / Real(nodeCount * (nodeCount - 1));
		QuadraturePoint edgePoint = { -Real(1), baseWeight };
		results.push_back(edgePoint);
		for (int pointIndex = 1; pointIndex < nodeCount - 1; pointIndex++) {
			QuadraturePoint point;
			Real theta, rootEquationDerivative, rootEquation;
			Real newTheta = pi * (Real(1) - Real(2 * pointIndex + 1) / Real(2) / Real(nodeCount));
			do {
				theta = newTheta;
				point.position = std::cos(theta);
				rootEquation = legendreP(nodeCount - 1, point.position, 1);
				rootEquationDerivative = legendreP(nodeCount - 1, point.position, 2) * std::sin(theta);
				newTheta = theta + rootEquation / rootEquationDerivative;
			} while (double(newTheta) != std::nextafter(double(theta), double(newTheta)));
			point.weight = baseWeight / sqr(legendreP(nodeCount - 1, point.position));
			//	printf( "----> %e\n", double(point.position));
			results.push_back(point);
		}
		edgePoint.position = -edgePoint.position;
		results.push_back(edgePoint);
		break;
	}
	}
	return results;
}

Matrix transformMatrix1D(TransformDirection transformDirection, int modeCount, int nodeCount, Quadrature quadratureType, bool isDerivative = false) {
	Matrix transform;
	if (transformDirection == TransformDirection::forward) {
		transform = createMatrix(modeCount, nodeCount);
		auto const quadratureRules = gaussQuadrature(nodeCount, quadratureType);
		for (int modeIndex = 0; modeIndex < modeCount; modeIndex++) {
			for (int nodeIndex = 0; nodeIndex < nodeCount; nodeIndex++) {
				Real const position = quadratureRules[nodeIndex].position;
				Real const weight = quadratureRules[nodeIndex].weight;
				Real const inverseMass = Real(2 * modeIndex + 1) / Real(2);
				transform[modeIndex][nodeIndex] = weight * inverseMass * legendreP(modeIndex, position);
			}
		}
	} else {
		assert(!isDerivative);
		transform = createMatrix(nodeCount, modeCount);
		auto const quadratureRules = gaussQuadrature(nodeCount, quadratureType);
		for (int modeIndex = 0; modeIndex < modeCount; modeIndex++) {
			for (int nodeIndex = 0; nodeIndex < nodeCount; nodeIndex++) {
				Real const position = quadratureRules[nodeIndex].position;
				transform[nodeIndex][modeIndex] = legendreP(modeIndex, position);
			}
		}
	}
	return transform;
}

std::string matrixVectorProduct(std::vector<std::string> const &v, Matrix const &A, std::vector<std::string> const &x) {
	std::string code;
	int const M = x.size();
	int const N = A.size();
	assert((int ) v.size() == N);
	assert((int ) A[0].size() == M);

	for (int n = 0; n < N; n++) {
		if (v[n] == "") {
			continue;
		}
		bool first = true;
		char *ptr;
		for (int m = 0; m < M; m++) {
			Real const C = A[n][m];
			if (std::abs(C) < tiny) {
				continue;
			}
			std::string cons;
			if (std::abs(std::abs(C) - Real(1)) >= tiny) {
				cons = getConstant(std::abs(C));
			}
			if (first) {
				if (std::abs(C - Real(+1)) < tiny) {
					asprintf(&ptr, "%s = %s;\n", v[n].c_str(), x[m].c_str());
				} else if (std::abs(C - Real(-1)) < tiny) {
					asprintf(&ptr, "%s = -%s;\n", v[n].c_str(), x[m].c_str());
				} else if (C < Real(0)) {
					assert(cons != "");
					asprintf(&ptr, "%s = -%s * %s;\n", v[n].c_str(), cons.c_str(), x[m].c_str());
				} else/*if (C > Real(0))*/{
					assert(cons != "");
					asprintf(&ptr, "%s = %s * %s;\n", v[n].c_str(), cons.c_str(), x[m].c_str());
				}
				first = false;
			} else {
				if (std::abs(C - Real(+1)) < tiny) {
					asprintf(&ptr, "%s += %s;\n", v[n].c_str(), x[m].c_str());
				} else if (std::abs(C - Real(-1)) < tiny) {
					asprintf(&ptr, "%s -= %s;\n", v[n].c_str(), x[m].c_str());
				} else if (C < Real(0)) {
					assert(cons != "");
					asprintf(&ptr, "%s -= %s * %s;\n", v[n].c_str(), cons.c_str(), x[m].c_str());
				} else/*if (C > Real(0))*/{
					assert(cons != "");
					asprintf(&ptr, "%s += %s * %s;\n", v[n].c_str(), cons.c_str(), x[m].c_str());
				}
			}
			code += indent + ptr;
			free(ptr);
		}
		if (first) {
			asprintf(&ptr, "%s  = %s(0);\n", v[n].c_str(), (std::string("Type")).c_str());
			code += indent + ptr;
			free(ptr);
		}
	}
	return code;
}

Matrix permutationMatrix(int modeCount, int triIndexCount) {
	int N = binco(modeCount + triIndexCount - 1, triIndexCount);
	int M = ipow(modeCount, triIndexCount);
	int L = ipow(modeCount + 1, dimensionCount - triIndexCount);
	assert(N);
	assert(M);
	Matrix P = createMatrix(N, M);
	if (triIndexCount) {
		for (int n = 0; n < N; n++) {
			auto const i = flat2triangular(triIndexCount, n);
			int m = 0;
			for (int d = 0; d < triIndexCount; d++) {
				m = modeCount * m + i[d];
			}
			printf("%i %i %i %i\n", n, N, m, M);
			P[n][m] = Real(1);
		}
	} else {
		P[0][0] = 1;
	}
	P = kroneckerProduct(identityMatrix(L), P);
	if (!triIndexCount) {
		std::cout << matrixToString(P);
	}
	return P;
}

Matrix permutationToMatrix(std::vector<int> const &permutation) {
	auto P = createMatrix(permutation.size());
	int const N = P.size();
	for (int n = 0; n < N; n++) {
		P[n][permutation[n]] = Real(1);
	}
	return P;
}

Matrix transformMatrix(TransformDirection transformDirection, int transformDimension) {
	auto A = identityMatrix(1);
	Matrix B;
	std::vector<int> sizes;
	for (int thisDimension = dimensionCount - 1; thisDimension >= 0; thisDimension--) {
		int const nodeCount = modeCount + 1;
		if (thisDimension < transformDimension) {
			A = kroneckerProduct(identityMatrix(nodeCount), A);
		} else if (thisDimension > transformDimension) {
			A = kroneckerProduct(identityMatrix(modeCount), A);
		} else/*if( thisDimension == transformDimension*/{
			auto const transform = transformMatrix1D(transformDirection, modeCount, nodeCount, Quadrature::gaussLobatto);
			A = kroneckerProduct(transform, A);
		}
	}
	int const columnCount = A[0].size();
	int const rowCount = A.size();
	if (transformDirection == TransformDirection::forward) {
		auto const P1 = permutationMatrix(modeCount, dimensionCount - transformDimension - 1);
		auto const P2 = permutationMatrix(modeCount, dimensionCount - transformDimension);
		A = matrixMultiply(P2, matrixMultiply(A, matrixTranspose(P1)));
	} else {
		auto const P1 = permutationMatrix(modeCount, dimensionCount - transformDimension);
		auto const P2 = permutationMatrix(modeCount, dimensionCount - transformDimension - 1);
		A = matrixMultiply(P2, matrixMultiply(A, matrixTranspose(P1)));
	}
	return A;
}

std::vector<std::string> generateVariableNames(std::string const &name, int count, int offset = 0) {
	std::vector<std::string> names(count);
	for (int nameIndex = 0; nameIndex < count; nameIndex++) {
		names[nameIndex] = name + "[" + std::to_string(nameIndex + offset) + "]";
	}
	return names;
}

std::vector<std::string> generateVariableNames(char name, int count) {
	return generateVariableNames(std::string(1, name), count);
}

std::string generateArgumentDeclaration(char name, int count, std::string type) {
	return std::string("std::array<") + type + std::string(", ") + std::to_string(count) + std::string("> const& ") + std::string(1, name);
}

std::string generateVariableDeclaration(char name, int count, std::string type) {
	return indent + std::string("std::array<") + type + std::string(", ") + std::to_string(count) + std::string("> ") + std::string(1, name) + ";\n";
}

Matrix stiffnessMatrix(int modeCount, int direction) {
	auto S = createMatrix(modeCount);
	auto const rules = gaussQuadrature(modeCount, Quadrature::gaussLegendre);
	for (int n = 0; n < modeCount; n++) {
		for (int m = 0; m < modeCount; m++) {
			S[n][m] = Real(0);
			for (auto const &rule : rules) {
				Real const x = rule.position;
				Real const w = rule.weight;
				S[n][m] += w * legendreP(n, x, 1) * legendreP(m, x);
			}
		}
	}
	auto S3d = identityMatrix(1);
	for (int dimension = 0; dimension < dimensionCount; dimension++) {
		S3d = kroneckerProduct((direction == dimension) ? S : identityMatrix(modeCount), S3d);
	}
	std::cout << matrixToString(S3d) << "\n";
	return S3d;
}

Matrix inverseMassMatrix(int modeCount) {
	auto iM = createMatrix(modeCount);
	for (int m = 0; m < modeCount; m++) {
		iM[m][m] = Real(2 * m + 1) / Real(2);
	}
	auto iM3d = identityMatrix(1);
	for (int dimension = 0; dimension < dimensionCount; dimension++) {
		iM3d = kroneckerProduct(iM, iM3d);
	}
	return iM3d;
}

static auto genArray(int count) {
	return std::string("std::array<Type, ") + std::to_string(count) + ">";
}
;

std::string generateStiffnessMatrix() {
	std::string code;
	std::vector<std::string> inputs, outputs;
	std::array<Matrix, dimensionCount> S;
	for (int d = 0; d < dimensionCount; d++) {
		S[d] = stiffnessMatrix(modeCount, d);
		auto const permute = permutationMatrix(modeCount, dimensionCount);
		S[d] = matrixMultiply(S[d], matrixTranspose(permute));
		S[d] = matrixMultiply(permute, S[d]);
	}
	code += indent + "template<RealNumberType Type>\n";
	code +=
			indent + "TrialSpace<Type>::AnalysisVector TrialSpace<Type>::applyStiffnessOperator(int dimension, TrialSpace<Type>::AnalysisVector const& input) {\n";
	indent++;
	inputs = generateVariableNames("input", S[0][0].size());
	outputs = generateVariableNames("output", S[0].size());
	code += indent + "AnalysisVector output;\n";
	code += std::string(indent);
	for (int dim = 0; dim < dimensionCount; dim++) {
		code += "if( dimension == " + std::to_string(dim) + " ) {\n";
		indent++;
		code += matrixVectorProduct(outputs, S[dim], inputs);
		indent--;
		if (dim + 1 == dimensionCount) {
			code += indent + "}\n";
		} else {
			code += indent + "} else ";
		}
	}
	code += indent + "return output;\n";
	indent--;
	code += indent + "}\n\n";
	return code;
}

std::string generateInverseMassMatrix() {
	std::string code;
	std::vector<std::string> inputs, outputs;
	auto iM = inverseMassMatrix(modeCount);
	auto const permute = permutationMatrix(modeCount, dimensionCount);
	iM = matrixMultiply(iM, matrixTranspose(permute));
	iM = matrixMultiply(permute, iM);
	code += indent + "template<RealNumberType Type>\n";
	code += indent + "TrialSpace<Type>::AnalysisVector TrialSpace<Type>::applyInverseMassOperator(TrialSpace<Type>::AnalysisVector const& input) {\n";
	indent++;
	inputs = generateVariableNames("input", iM[0].size());
	outputs = generateVariableNames("output", iM.size());
	code += indent + "AnalysisVector output;\n";
	code += matrixVectorProduct(outputs, iM, inputs);
	code += indent + "return output;\n";
	indent--;
	code += indent + "}\n\n";
	return code;
}

std::string generateTransform(TransformDirection transformDirection) {
	std::string code, deferredCode;
	std::vector<Matrix> matrixFactors;
	std::vector<int> arraySizes;
	std::vector<std::string> inputs, outputs;
	inputs = decltype(inputs)();
	outputs = decltype(outputs)();
	bool const isForward = transformDirection == TransformDirection::forward;
	std::string functionName = isForward ? "analyze" : "synthesize";
	std::string outputSize = isForward ? std::to_string(modeCount3d) : std::to_string(nodeCount3d);
	std::string inputSize = isForward ? std::to_string(nodeCount3d) : std::to_string(modeCount3d);
	int inCount = isForward ? nodeCount3d : modeCount3d;
	int outCount = isForward ? modeCount3d : nodeCount3d;
	int bufferSize, currentSize, bufferOffset;
	outputs = generateVariableNames("input", inCount);
	code += "\n";
	if (isForward) {
		code += indent + "template<RealNumberType Type>\n";
		code += indent + "TrialSpace<Type>::AnalysisVector TrialSpace<Type>::analyze(TrialSpace<Type>::SynthesisVector const& input) {\n";
	} else {
		code += indent + "template<RealNumberType Type>\n";
		code += indent + "TrialSpace<Type>::SynthesisVector TrialSpace<Type>::synthesize(TrialSpace<Type>::AnalysisVector const& input) {\n";
	}
	arraySizes = std::vector<int>(1, inCount);
	auto const direction = isForward ? TransformDirection::forward : TransformDirection::backward;
	for (int transformDimension = 0; transformDimension < dimensionCount; transformDimension++) {
		matrixFactors.push_back(transformMatrix(direction, transformDimension));
	}
	if (isForward) {
		std::reverse(matrixFactors.begin(), matrixFactors.end());
	}
	for (int transformDimension = 0; transformDimension < dimensionCount; transformDimension++) {
		auto const sz = matrixFactors[transformDimension].size();
		arraySizes.push_back(sz);
	}
	bufferSize = *(std::max_element(arraySizes.begin() + 1, arraySizes.end() - 1));
	currentSize = arraySizes[0] + arraySizes[1];
	for (int i = 0; i < dimensionCount - 2; i++) {
		currentSize -= arraySizes[i];
		currentSize += arraySizes[i + 2];
		bufferSize = std::max(bufferSize, currentSize);
	}
	indent++;
	code += indent + genArray(bufferSize) + " buffer;\n";
	bufferOffset = 0;
	for (int transformDimension = 0; transformDimension < dimensionCount; transformDimension++) {
		auto const &transformMatrix = matrixFactors[transformDimension];
		inputs = std::move(outputs);
		if (transformDimension != (dimensionCount - 1)) {
			outputs = generateVariableNames("buffer", transformMatrix.size(), bufferOffset);
			if (bufferOffset == 0) {
				bufferOffset = transformMatrix.size();
			} else {
				bufferOffset = 0;
			}
		} else {
			outputs = generateVariableNames("output", outCount);
			code += indent + "std::array<Type, " + outputSize + "> output;\n";
		}
		deferredCode += matrixVectorProduct(outputs, transformMatrix, inputs);
	}
	code += deferredCode;
	code += indent + "return output;\n";
	indent--;
	code += indent + "}\n\n";
	return code;
}

int main(int, char*[]) {
//	constexpr int D = 3;
//	for (int flat = 0; flat < 51; flat++) {
//		std::array<int, D> i = flat2triangular<D>(flat);
//		int const flat2 = triangular2flat<D>(i);
//		printf("%i: ", flat);
//		for (int j = 0; j < D; j++) {
//			printf("%i ", i[j]);
//		}
//		printf("(%i)\n", flat2);
//		assert(flat2 == flat);
//	}
//	return 0;

	auto const analyzeCode = generateTransform(TransformDirection::forward);
	auto const synthesizeCode = generateTransform(TransformDirection::backward);
	auto const inverseMassCode = generateInverseMassMatrix();
	auto const stiffnessCode = generateStiffnessMatrix();
	std::ostringstream code;
	code << indent << "#pragma once\n";
	code << indent << "\n";
	code << indent << "#include <array>\n";
	code << indent << "#include <concepts>\n";
	code << indent << "#include <type_traits>\n";
	code << indent << "#include <experimental/simd>\n";
	code << indent << "\n";
	code << indent << "template<typename Type>\n";
	code << indent << "concept RealNumberType = (\n";
	indent++;
	code << indent << "std::floating_point<Type> ||\n";
	code << indent << "(std::experimental::is_simd<Type>::value && std::floating_point<typename Type::value_type>)\n";
	indent--;
	code << indent << ");\n";
	code << indent << "\n";
	std::string const constants = std::string(getConstant.getCode());

	code << "template<RealNumberType Type, bool = std::is_literal_type_v<Type>>\n"
			"struct TrialSpaceConstants;\n\n"
			"template<RealNumberType Type>\n"
			"struct TrialSpaceConstants<Type, false> {\n"
			"   inline static const " << constants << "\n"
			"   };\n"
			"};\n\n"
			"template<RealNumberType Type>\n"
			"struct TrialSpaceConstants<Type, true> {\n"
			"   inline static constexpr " << constants << "\n"
			"   };\n"
			"};\n\n"
			"";
	code << indent << "template<RealNumberType Type>\n";
	code << indent << "struct TrialSpace : public TrialSpaceConstants<Type> {\n\n";
	indent++;
	code << indent << "static constexpr int modalSize = " << modeCount3d << ";\n";
	code << indent << "static constexpr int nodalSize = " << nodeCount3d << ";\n";
//	code << modalIndicesCode;
	code << indent << "using TrialSpaceConstants<Type>::constant;\n";
	code << indent << "using AnalysisVector = std::array<Type, modalSize>;\n";
	code << indent << "using SynthesisVector = std::array<Type, nodalSize>;\n\n";
	code << indent << "static AnalysisVector analyze(SynthesisVector const&);\n";
	code << indent << "static SynthesisVector synthesize(AnalysisVector const&);\n";
	code << indent << "static AnalysisVector applyInverseMassOperator(AnalysisVector const&);\n";
	code << indent << "static AnalysisVector applyStiffnessOperator(int, AnalysisVector const&);\n\n";
	indent--;
	code << indent << "};\n";
	code << indent << "\n";
	code << analyzeCode;
	code << synthesizeCode;
	code << inverseMassCode;
	code << stiffnessCode;
	toFile(code.str(), "./generated_source/TrialSpace.hpp");
	return 0;
}
