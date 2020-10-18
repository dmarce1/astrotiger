#pragma once
#include <astrotiger/defs.hpp>

#include <string>

class options {
public:
	std::string config_file;
	int nhydro;
	int ngroup;
	int max_box;
	std::vector<int> bw;
	int ngrid;
	int window;
	int max_level;
	std::string problem;
	double refine_slope;

	template<class Arc>
	void serialize(Arc &arc, unsigned) {
		arc & bw;
		arc & config_file;
		arc & nhydro;
		arc & ngroup;
		arc & refine_slope;
		arc & max_box;
		arc & ngrid;
		arc & window;
		arc & max_level;
	}

	static void set(options);
	bool process_options(int argc, char *argv[]);
};

#ifndef OPTIONS_CPP
extern options opts;
#endif
