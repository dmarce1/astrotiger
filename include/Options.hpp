/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_OPTIONS_HPP_
#define INCLUDE_OPTIONS_HPP_

#include "Definitions.hpp"

#include <string>
#include <type_traits>

#define OPTIONS                                \
	OPTION(dim_count, int,         2,     ""); \
	OPTION(max_level, int,         2,     ""); \
	OPTION(max_time,  double,      1.0,   ""); \
	OPTION(problem,   std::string, "sod", ""); \
	OPTION(scale,     double,      1.0,   ""); \
	OPTION(verbose,   int,         1,     ""); \

class Options {
	bool status;
	std::string config;
#define OPTION(name, type, dflt, desc) type name = dflt;
	OPTIONS
#undef OPTION
	template<typename T>
	static std::string optionString(std::string const &name, T value, size_t);
public:
	template<typename Arc, unsigned>
	void serialize(Arc &&arc, unsigned) {
#define OPTION(name, type, dflt, desc) arc & name;
		OPTIONS
#undef OPTION
	}
	Options(int, char*[]);
	operator std::string() const;
};

#endif /* INCLUDE_OPTIONS_HPP_ */
