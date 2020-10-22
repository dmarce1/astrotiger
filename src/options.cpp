#define OPTIONS_CPP
#include <astrotiger/defs.hpp>
#include <astrotiger/hpx.hpp>
#include <astrotiger/options.hpp>
#include <boost/program_options.hpp>

#include <fstream>
#include <iostream>

HPX_PLAIN_ACTION(options::set, set_options_action);

options opts;

void options::set(options o) {
	opts = o;
}

bool options::process_options(int argc, char *argv[]) {
	namespace po = boost::program_options;

	po::options_description command_opts("options");

	command_opts.add_options() //
	("help", "produce help message") //
	("config_file", po::value < std::string > (&config_file)->default_value(""), "configuration file") //
	("problem", po::value < std::string > (&problem)->default_value("blast"), "Problem 1 - sod 2 - blast\n") //
	("max_box", po::value<int>(&max_box)->default_value(32), "maximum (box volume)^(1/3)") //
	("window", po::value<int>(&window)->default_value(1), "refinement window size") //
	("min_box", po::value<int>(&min_box)->default_value(8), "minimum (box volume)^(1/3)") //
	("max_level", po::value<int>(&max_level)->default_value(3), "maximum refinement level") //
	("tmax", po::value<double>(&tmax)->default_value(0.25), "maximum simulation time") //
	("cfl", po::value<double>(&cfl)->default_value(0.2), "cfl factor") //
	("efficiency", po::value<double>(&efficiency)->default_value(0.5), "refinement efficiency") //
	("gamma", po::value<double>(&gamma)->default_value(5.0 / 3.0), "fluid gamma") //
	("refine_slope", po::value<double>(&refine_slope)->default_value(0.1), "refinement slope criteria") //
	("ngroup", po::value<int>(&ngroup)->default_value(1), "number of frequency groups") //
			;

	boost::program_options::variables_map vm;
	po::store(po::parse_command_line(argc, argv, command_opts), vm);
	po::notify(vm);
	if (vm.count("help")) {
		std::cout << command_opts << "\n";
		return false;
	}
	if (!config_file.empty()) {
		std::ifstream cfg_fs { vm["config_file"].as<std::string>() };
		if (cfg_fs) {
			po::store(po::parse_config_file(cfg_fs, command_opts), vm);
		} else {
			printf("Configuration file %s not found!\n", config_file.c_str());
			return false;
		}
	}
	po::notify(vm);

	nhydro = 2 + NDIM;
	hbw = 2;
	max_bw = 2;
	window = 1;
	nrk = 2;

	const auto loc = hpx::find_all_localities();
	const auto sz = loc.size();
	std::vector<hpx::future<void>> futs;
	set(*this);
	for (int i = 1; i < sz; i++) {
		futs.push_back(hpx::async<set_options_action>(loc[i], *this));
	}
	hpx::wait_all(futs.begin(), futs.end());
#define SHOW( opt ) std::cout << std::string( #opt ) << " = " << std::to_string(opt) << '\n';
#define SHOW_STR( opt ) std::cout << std::string( #opt ) << " = " << opt << '\n';
	return true;
}
