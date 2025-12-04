/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#include "Definitions.hpp"

#include <cmath>
#include <fstream>

#include <hpx/modules/runtime_distributed.hpp>

#include "Options.hpp"

template<typename T>
std::string Options::optionString(std::string const &name, T value, size_t nameLen) {
	std::string str = name;
	while (str.size() < nameLen) {
		str += " ";
	}
	str += " = ";
	if constexpr (std::is_same<T, std::string>::value) {
		str += value;
	} else if constexpr (std::is_same<T, double>::value) {
		str += std::to_string(value);
	} else if constexpr (std::is_same<T, int>::value) {
		str += std::to_string(value);
	} else {
		static_assert(false, "Unknown option type\n");
	}
	str += "\n";
	return str;
}

Options::operator std::string() const {
	using std::max;
	size_t maxNameLen = 0;
	std::string str;
#define OPTION(name, type, dflt, desc) maxNameLen = std::max(maxNameLen, strlen(#name));
	OPTIONS
#undef OPTION
#define OPTION(name, type, dflt, desc) str += optionString<type>(#name, name, maxNameLen);
	OPTIONS
#undef OPTION
	return str;
}

Options::Options(int argc, char *argv[]) :
		status(false) {
	using namespace hpx::program_options;
	options_description options("options");
#define OPTION(name, type, dflt, desc) \
    options.add_options() \
        (#name, value<type>()->default_value(dflt), desc);
	OPTIONS
	;
#undef OPTION
	options.add_options()("help", "show help message")("config", value<std::string>()->default_value(""), "Path to .ini configuration file");
	variables_map vm;
	store(parse_command_line(argc, argv, options), vm);
	notify(vm);
	if (vm.count("help")) {
		std::cout << options << "\n";
		return;
	}
	std::string cfg = vm["config"].as<std::string>();
	if (!cfg.empty()) {
		std::ifstream cfgStream(cfg);
		if (!cfgStream) {
			std::cerr << "Unable to read configuration file: " << cfg << "\n";
			return;
		}
		store(parse_config_file(cfgStream, options), vm);
		notify(vm);
	}
	status = true;
}
