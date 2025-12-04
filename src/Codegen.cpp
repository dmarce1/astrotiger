/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#include "Indent.hpp"

#include <algorithm>
#include <bitset>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>

constexpr int maxRank = 4;

Indent indent;

std::string codeFreeIndex() {
	std::string code;
	code += indent + "template<char C>\n";
	code += indent + "struct FreeIndex {\n";
	indent++;
	code += indent + "static constexpr char value = C;\n";
	indent--;
	code += indent + "};\n";
	code += indent + "\n\n";
	return code;
}

using PairList = std::vector<std::pair<size_t, size_t>>;

std::string codeTensor(int rank) {
	std::string code;
	code += indent + "template<typename T, int D>\n";
	code += indent + "struct Tensor" + std::to_string(rank) + " {\n";
	indent++;
	// T operator() const
	code += indent + "T operator()(";
	for (int r = 0; r < rank; r++) {
		code += (r ? std::string(", ") : std::string("")) + "int " + std::string(1, 'i' + r);
	}
	code += ") const {\n";
	indent++;
	code += indent + "return data_[flatten(";
	for (int r = 0; r < rank; r++) {
		code += (r ? std::string(", ") : std::string("")) + std::string(1, 'i' + r);
	}
	code += ")];\n";
	indent--;
	code += indent + "}\n";
	// T& operator()
	code += indent + "T& operator()(";
	for (int r = 0; r < rank; r++) {
		code += (r ? std::string(", ") : std::string("")) + "int " + std::string(1, 'i' + r);
	}
	code += ") {\n";
	indent++;
	code += indent + "return data_[flatten(";
	for (int r = 0; r < rank; r++) {
		code += (r ? std::string(", ") : std::string("")) + std::string(1, 'i' + r);
	}
	code += ")];\n";
	indent--;
	code += indent + "}\n";
	// auto operator()
	code += indent + "template<";
	for (int r = 0; r < rank; r++) {
		code += (r ? std::string(", ") : std::string("")) + "char " + std::string(1, 'I' + r);
	}
	code += ">\n";
	code += indent + "auto operator()(";
	for (int r = 0; r < rank; r++) {
		code += (r ? std::string(", ") : std::string("")) + "FreeIndex<" + std::string(1, 'I' + r) + ">";
	}
	code += ") {\n";
	indent++;
	code += indent + "return Tensor" + std::to_string(rank) + "Expression<Tensor" + std::to_string(rank) + "&, D, ";
	for (int r = 0; r < rank; r++) {
		code += (r ? std::string(", ") : std::string("")) + std::string(1, 'I' + r);
	}
	code += ">(*this);\n";
	indent--;
	code += indent + "}\n";
	// traces
	{
		std::vector<char> free(rank);
		for (int n = 1; n <= rank / 2; n++) {
			char c = 'I';
			for (int k = 0; k < rank; k++) {
				free[k] = c;
				if ((k & 1) || (k >= 2 * n)) {
					c++;
				}
			}
			std::set<PairList> pairs;
			do {
				std::vector<std::pair<size_t, size_t>> list;
				size_t count = 0;
				for (int j = 0; j < rank; j++) {
					for (int k = j + 1; k < rank; k++) {
						if (free[k] == free[j]) {
							list.push_back( { j, k });
							count++;
						}
					}
				}
				if (pairs.find(list) == pairs.end()) {
					pairs.insert(list);
					int const traceRank = rank - 2 * count;
					code += indent + "template<";
					for (size_t r = 0; r < rank - count; r++) {
						code += (r ? std::string(", ") : std::string("")) + "char " + std::string(1, 'I' + r);
					}
					code += ">\n";
					code += indent + "auto operator()(";
					for (int r = 0; r < rank; r++) {
						code += (r ? std::string(", ") : std::string("")) + "FreeIndex<" + free[r] + ">";
					}
					code += ") {\n";
					indent++;
					code += indent + "auto const lambda = [this](";
					for (int i = 0; i < traceRank; i++) {
						code += (i ? std::string(", ") : std::string("")) + "int " + std::string(1, 'i' + count + i);
					}
					code += ") {\n";
					indent++;
					code += indent + "T sum = T(0);\n";
					for (int i = 0; i < count; i++) {
						auto const idx = std::string(1, 'i' + i);
						code += indent + "for( int " + idx + "; " + idx + " < D; " + idx + "++ ) {\n";
						indent++;
					}
					code += indent + "sum += (*this)(";
					for (int r = 0; r < rank; r++) {
						code += (r ? std::string(", ") : std::string("")) + std::string(1, std::tolower(free[r]));
					}
					code += ");\n";
					for (int i = 0; i < count; i++) {
						indent--;
						code += indent + "}\n";
					}
					code += indent + "return sum;\n";
					indent--;
					code += indent + "};\n";
					code += indent + "return Tensor" + std::to_string(traceRank) + "Expression<decltype(lambda), D";
					for (int r = 0; r < traceRank; r++) {
						code += ", " + std::string(1, 'I' + count + r);
					}
					code += ">(*this);\n";
					indent--;
					code += indent + "}\n";
				}
			} while (std::next_permutation(free.begin(), free.end()));
		}
	}
	code += "private:\n";
	code += indent + "static constexpr int size = ";
	for (int r = 0; r < rank; r++) {
		code += (r ? std::string(" * ") : std::string("")) + "D";
	}
	code += ";\n";
	// flatten
	code += indent + "int flatten(";
	for (int r = 0; r < rank; r++) {
		code += (r ? std::string(", ") : std::string("")) + "int " + std::string(1, 'i' + r);
	}
	code += ") const {\n";
	indent++;
	code += indent + "return ";
	for (int r = 0; r < rank; r++) {
		if (r) {
			code += " + D * ";
			if (r + 1 < rank) {
				code += "(";
			}
		}
		code += std::string(1, 'i' + rank - 1 - r);
	}
	for (int r = 2; r < rank; r++) {
		code += ")";
	}
	code += ";\n";
	indent--;
	code += indent + "}\n";
	code += indent + "std::array<T, size> data_;\n";
	indent--;
	code += indent + "};\n";
	code += indent + "\n\n";
	return code;
}

auto codeTensorExpression(int rank) {
	std::string code;
	std::string impl;
	std::string tmpl;
	tmpl += "template<typename H, int D";
	for (int r = 0; r < rank; r++) {
		tmpl += ", char " + std::string(1, 'I' + r);
	}
	tmpl += ">\n";
	code += indent + tmpl;
	std::string const name = "TensorExpression" + std::to_string(rank);
	code += indent + "struct " + name + " {\n";
	indent++;
	code += indent + "" + name + "(H const& h) : handle_(h) {\n";
	code += indent + "}\n";
	std::string astr[] = { "=", "+=", "-=" };
	for (int ai = 0; ai < 3; ai++) {
		std::vector<char> indices(rank);
		std::iota(indices.begin(), indices.end(), 'I');
		do {
			code += indent + "template<typename F>\n";
			code += indent + "auto& operator" + astr[ai] + "(TensorExpression" + std::to_string(rank) + "<F, D";
			for (size_t r = 0; r < rank; r++) {
				code += ", " + std::string(1, indices[r]);
			}
			code += "> const& other) {\n";
			indent++;
			for (int r = 0; r < rank; r++) {
				std::string idx(1, 'i' + r);
				code += indent + "for(int " + idx + "; " + idx + " < D; " + idx + "++) {\n";
				indent++;
			}
			code += indent + "(*this)(";
			std::string comma = "";
			for (int r = 0; r < rank; r++) {
				std::string idx(1, 'i' + r);
				code += comma + idx;
				comma = ", ";
			}
			code += ") " + astr[ai] + " other(";
			comma = "";
			for (int r = 0; r < rank; r++) {
				std::string idx(1, tolower(indices[r]));
				code += comma + idx;
				comma = ", ";
			}
			code += ");\n";
			for (int r = 0; r < rank; r++) {
				indent--;
				code += indent + "}\n";
			}
			code += indent + "return *this;\n";
			indent--;
			code += indent + "}\n";
		} while (std::next_permutation(indices.begin(), indices.end()));
	}
	if (rank > 0) {
		std::vector<char> indices(rank);
		std::iota(indices.begin(), indices.end(), 'I');
		for (unsigned otherRank = 1; otherRank <= maxRank; otherRank++) {
			for (unsigned contract = 0; contract < unsigned(1 << rank); contract++) {
				std::bitset < 32 > bits(contract);
				if ((bits.count() > otherRank) || (otherRank + rank - 2 * bits.count() > maxRank)) {
					continue;
				}
				std::vector<char> otherIndices(otherRank);
				std::iota(otherIndices.begin(), otherIndices.end(), 'I' + rank);
				int next = 0;
				for (int r = 0; r < rank; r++) {
					if (bits[r]) {
						otherIndices[next++] = indices[r];
					}
				}
				for (size_t r = 0; r < otherRank; r++) {
					if (otherIndices[r] >= 'I' + rank) {
						otherIndices[r] -= bits.count();
					}
				}
				while (std::next_permutation(otherIndices.begin(), otherIndices.end())) {
				}
				std::set<PairList> pairs;
				do {
					std::vector<std::pair<size_t, size_t>> list;
					size_t count = 0;
					for (int j = 0; j < rank; j++) {
						for (int k = 0; k < otherRank; k++) {
							if (indices[j] == otherIndices[k]) {
								list.push_back( { j, k });
								count++;
							}
						}
					}
					if (pairs.find(list) == pairs.end()) {
						pairs.insert(list);\
						std::string pack;
						code += indent + "template<typename F";
						pack += "<F, D";
						for (size_t r = 0; r < otherRank; r++) {
							pack += ", " + std::string(1, otherIndices[r]);
							if (otherIndices[r] < 'I' + rank) {
								continue;
							}
							code += ", char " + std::string(1, otherIndices[r]);
						}
						pack += ">";
						code += ">\n";
						code += indent + "auto operator*(TensorExpression" + std::to_string(otherRank) + pack + " const& other);\n";
						indent--;
						impl += tmpl;
						impl += indent + "template<typename F";
						for (size_t r = 0; r < otherRank; r++) {
							if (otherIndices[r] < 'I' + rank) {
								continue;
							}
							impl += ", char " + std::string(1, otherIndices[r]);
						}
						impl += ">\n";
						impl += indent + "auto ";
						impl += name + "<H, D";
						for (int r = 0; r < rank; r++) {
							impl += ", " + std::string(1, 'I' + r);
						}
						impl += ">::operator*(TensorExpression" + std::to_string(otherRank) + pack + " const& other) {\n";
						indent++;
						std::vector<char> fixed;
						std::vector<char> free;
						for (int r = 0; r < rank; r++) {
							if (std::find(otherIndices.begin(), otherIndices.end(), indices[r]) == otherIndices.end()) {
								fixed.push_back(indices[r]);
							} else {
								free.push_back(indices[r]);
							}
						}
						for (int r = 0; r < otherRank; r++) {
							if (std::find(indices.begin(), indices.end(), otherIndices[r]) == indices.end()) {
								fixed.push_back(otherIndices[r]);
							}
						}
						impl += indent + "auto const lambda = [this, other](";
						for (size_t r = 0; r < fixed.size(); r++) {
							impl += (r ? std::string(", ") : std::string("")) + "int " + std::string(1, tolower(fixed[r]));
						}
						impl += ") {\n";
						indent++;
						impl += indent + "using type = std::remove_cvref_t<decltype((*this)(";
						for (int r = 0; r < rank; r++) {
							impl += (r ? std::string(", ") : std::string("")) + "0";
						}
						impl += ") * other(";
						for (size_t r = 0; r < otherRank; r++) {
							impl += (r ? std::string(", ") : std::string("")) + "0";
						}
						impl += "))>;\n";
						impl += indent + "auto sum = type(0);\n";
						for (size_t r = 0; r < free.size(); r++) {
							std::string idx(1, tolower(free[r]));
							impl += indent + "for(int " + idx + " = 0; " + idx + " < D; " + idx + "++ ) {\n";
							indent++;
						}
						impl += indent + "sum += (*this)(";
						for (int r = 0; r < rank; r++) {
							impl += (r ? ", " : "") + std::string(1, 'i' + r);
						}
						impl += ") * other(";
						for (size_t r = 0; r < otherRank; r++) {
							impl += (r ? ", " : "") + std::string(1, tolower(otherIndices[r]));
						}
						impl += ");\n";
						for (size_t r = 0; r < free.size(); r++) {
							indent--;
							impl += indent + "}\n";
						}
						impl += indent + "return sum;\n";
						indent--;
						impl += indent + "};\n";
						size_t const traceRank = rank + otherRank - 2 * bits.count();
						impl += indent + "return Tensor" + std::to_string(traceRank) + "Expression<decltype(lambda), D";
						for (size_t r = 0; r < fixed.size(); r++) {
							impl += ", " + std::string(1, fixed[r]);
						}
						impl += ">(lambda);\n";
						indent--;
						impl += "}\n\n";
						indent++;
					}
				} while (std::next_permutation(otherIndices.begin(), otherIndices.end()));
			}
		}
	}
	code += "private:\n";
	code += indent + "H handle_;\n";
	indent--;
	code += indent + "};\n";
	code += indent + "\n\n";
	return std::pair<std::string, std::string>(code, impl);
}

int main(int argc, char *argv[]) {
	std::string code;
	std::string impl;

	code += indent + "#include <array>\n";
	code += indent + "\n\n";
	code += codeFreeIndex();
	for (int rank = 0; rank <= maxRank; rank++) {
		code += indent + "template<typename, int";
		for (int r = 0; r < rank; r++) {
			code += ", char";
		}
		code += ">\n";
		code += "struct TensorExpression" + std::to_string(rank) + ";\n";
		code += "\n";
	}
	code += "\n";
	for (int r = 0; r <= maxRank; r++) {
		auto const [c, i] = codeTensorExpression(r);
		code += c;
		impl += i;
	}
	for (int r = 1; r <= maxRank; r++) {
		code += codeTensor(r);
	}
	std::cout << code;
	std::cout << impl;
	return 0;
}
