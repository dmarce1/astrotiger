#define OPTIONS_CPP
#include <astrotiger/defs.hpp>
#include <astrotiger/hpx.hpp>
#include <astrotiger/options.hpp>
#include <boost/program_options.hpp>
#include <astrotiger/cosmos.hpp>
#include <astrotiger/fileio.hpp>

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
	("xmbnd", po::value < std::string > (&xmbnd)->default_value("periodic"), "") //
	("xpbnd", po::value < std::string > (&xpbnd)->default_value("periodic"), "") //
	("ymbnd", po::value < std::string > (&ymbnd)->default_value("periodic"), "") //
	("ypbnd", po::value < std::string > (&ypbnd)->default_value("periodic"), "") //
	("zmbnd", po::value < std::string > (&zmbnd)->default_value("periodic"), "") //
	("zpbnd", po::value < std::string > (&zpbnd)->default_value("periodic"), "") //
	("config_file", po::value < std::string > (&config_file)->default_value(""), "configuration file") //
	("problem", po::value < std::string > (&problem)->default_value("blast"), "Problem 1 - sod 2 - blast\n") //
	("hydro", po::value<bool>(&hydro)->default_value(true), "use hydro") //
	("self_gravity", po::value<bool>(&self_gravity)->default_value(false), "use self gravity") //
	("particles", po::value<bool>(&particles)->default_value(false), "use particles") //
	("gravity", po::value<bool>(&gravity)->default_value(false), "use gravity") //
	("nmulti", po::value<int>(&nmulti)->default_value(16), "multigrid solver iterations)") //
	("max_box", po::value<int>(&max_box)->default_value(32), "maximum (box volume)^(1/3)") //
	("order", po::value<int>(&order)->default_value(2), "integration order") //
	("window", po::value<int>(&window)->default_value(1), "refinement window size") //
	("min_box", po::value<int>(&min_box)->default_value(10), "minimum (box volume)^(1/3)") //
	("max_level", po::value<int>(&max_level)->default_value(3), "maximum refinement level") //
	("tmax", po::value<double>(&tmax)->default_value(1.00), "maximum simulation time") //
	("cfl", po::value<double>(&cfl)->default_value(0.2), "cfl factor") //
	("efficiency", po::value<double>(&efficiency)->default_value(0.10), "refinement efficiency") //
	("code_to_g", po::value<double>(&code_to_g)->default_value(1.988e43), "code to gram conversion") //
	("code_to_cm", po::value<double>(&code_to_cm)->default_value(3.0e25), "size of domain") //
	("code_to_cm_per_s", po::value<double>(&code_to_cm_per_s)->default_value(3e10), "code to cm/s") //
	("gamma", po::value<double>(&gamma)->default_value(5.0 / 3.0), "fluid gamma") //
	("refine_slope", po::value<double>(&refine_slope)->default_value(0.025), "refinement slope criteria") //
	("ngroup", po::value<int>(&ngroup)->default_value(1), "number of frequency groups") //
			;

	species = false;
	nspecies = 0;
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

	const auto str_to_bc_type = [](const std::string str) {
		if (str == "periodic") {
			return PERIODIC;
		} else if (str == "reflecting") {
			return REFLECTING;
		} else if (str == "outflow") {
			return OUTFLOW;
		} else {
			printf("Unknow boundary condition %s\n", str);
		}
	};
	bnd[0] = str_to_bc_type(xpbnd);
	bnd[1] = str_to_bc_type(xmbnd);
	if ( NDIM > 1) {
		bnd[2] = str_to_bc_type(ypbnd);
		bnd[3] = str_to_bc_type(ymbnd);
		if ( NDIM > 2) {
			bnd[4] = str_to_bc_type(zpbnd);
			bnd[5] = str_to_bc_type(zmbnd);
		}
	}
	G = 1;
	if (problem == "rt" && NDIM > 1) {
		bnd[2] = bnd[3] = REFLECTING;
		gamma = 7.0 / 5.0;
		gravity = true;
	} else if (problem == "sod") {
		gamma = 7.0 / 5.0;
	} else if (problem == "polytrope") {
		gravity = self_gravity = true;
		gamma = 5.0 / 3.0;
	} else if (problem == "part_test") {
		hydro = false;
		particles = true;
		gravity = self_gravity = true;
	} else if (problem == "cosmos") {
		hydro = gravity = self_gravity = particles = true;
		nspecies = 5;
		species = true;
		opts = *this;
		printf( "Reading input file\n");
		fileio_init_read();
		*this = opts;
	}

	nhydro = 4 + NDIM + nspecies;
	hbw = 2;
	gbw = 2;
	max_bw = std::max(hbw, window);
	max_bw = std::max(max_bw, gbw);
	nrk = order;
	alpha.resize(nrk);
	beta.resize(nrk);
	if (order == 1) {
		alpha[0] = beta[0] = 1.0;
		cfl = 0.999 / NDIM;
	} else if (order == 2) {
		cfl = 0.999 * (2.0 / 3.0) / NDIM;
		beta[0] = 1.0;
		beta[1] = 0.5;
		alpha[0] = 0.0;
		alpha[1] = 1.0;
	} else {
		printf("Order %i not supported\n");
		abort();
	}

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
	printf( "Sending H0 %e\n", H0);
	return true;
}
