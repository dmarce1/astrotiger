/*
 * Options.cpp
 *
 *  Created on: Jan 13, 2025
 *      Author: dmarce1
 */

#include "Options.hpp"

#include <fstream>

#include <hpx/hpx.hpp>

#define SHOW( opt ) show(#opt, opts.opt)

static void show(const char*, bool);
static void show(const char*, int);
static void show(const char*, double);
static void show(const char*, std::string);
static void setOptions(const Options &opts);
const Options& getOptions();

static Options globalOptions;

HPX_PLAIN_ACTION (setOptions);

const Options& getOptions() {
	return globalOptions;
}

bool processOptions(int argc, char *argv[]) {
	namespace po = hpx::program_options;
	Options opts;
	bool success;
	po::options_description commandOptions("options");
	commandOptions.add_options()                                                                            //
	("help", "produce help message")                                                                        //
	("configFile", po::value < std::string > (&(opts.configFile))->default_value(""), "configuration file") //
	("fluidGamma", po::value<double>(&(opts.fluidGamma))->default_value(5.0 / 3.0), "fluid gamma")          //
	("dualEnergyPressureSwitch", po::value<double>(&(opts.dualEnergyPressureSwitch))->default_value(0.001),
			"dual energy pressure switch")                                                                  //
	("dualEnergyUpdateSwitch", po::value<double>(&(opts.dualEnergyUpdateSwitch))->default_value(0.1),
			"dual energy update switch")                                                                    //
			;
	printf("Processing options\n");
	po::variables_map variableMap;
	po::store(po::parse_command_line(argc, argv, commandOptions), variableMap);
	po::notify(variableMap);
	if (variableMap.count("help")) {
		std::cout << commandOptions << "\n";
		success = false;
	} else {
		if (!opts.configFile.empty()) {
			std::ifstream configFile { variableMap["configFile"].as<std::string>() };
			if (configFile) {
				po::store(po::parse_config_file(configFile, commandOptions), variableMap);
				success = true;
			} else {
				printf("Configuration file %s not found!\n", opts.configFile.c_str());
				return false;
			}
		} else {
			success = true;
		}
	}

	if (success) {
		po::notify(variableMap);
	}
	SHOW(configFile);
	SHOW(dualEnergyPressureSwitch);
	SHOW(dualEnergyUpdateSwitch);
	SHOW(fluidGamma);
	setOptions(opts);
	return success;
}

void setOptions(const Options &opts) {
	if (hpx::get_locality_id() == 0) {
		std::vector<hpx::future<void>> futs;
		const auto localities = hpx::find_all_localities();
		int const cnt = localities.size();
		for (int i = 1; i < cnt; i++) {
			futs.push_back(hpx::async < setOptions_action > (localities[i], opts));
		}
		hpx::wait_all(futs.begin(), futs.end());
	}
	globalOptions = opts;
}

static void show(const char *name, bool opt) {
	printf("%-20s: %c\n", name, opt ? 'T' : 'F');
}

static void show(const char *name, int opt) {
	printf("%-20s: %i\n", name, opt);
}

static void show(const char *name, double opt) {
	printf("%-20s: %e\n", name, opt);
}

static void show(const char *name, std::string opt) {
	printf("%-20s: %s\n", name, opt.c_str());
}
