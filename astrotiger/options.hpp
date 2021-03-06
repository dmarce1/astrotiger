#pragma once
#include <astrotiger/defs.hpp>

#include <string>
#include <vector>

class options {
public:
	std::string config_file;
	int nhydro;
	int ngroup;
	int max_box;
	int min_box;
	int hbw;
	int max_bw;
	int window;
	int max_level;
	int order;
	int nrk;
	std::string problem;
	double refine_slope;
	double gamma;
	double cfl;
	double tmax;
	double efficiency;
	std::vector<double> alpha;
	std::vector<double> beta;

	template<class Arc>
	void serialize(Arc &arc, unsigned) {
		arc & alpha;
		arc & beta;
		arc & order;
		arc & efficiency;
		arc & min_box;
		arc & tmax;
		arc & cfl;
		arc & nrk;
		arc & config_file;
		arc & nhydro;
		arc & ngroup;
		arc & hbw;
		arc & max_bw;
		arc & refine_slope;
		arc & max_box;
		arc & window;
		arc & max_level;
		arc & problem;
		arc & gamma;
	}

	static void set(options);
	bool process_options(int argc, char *argv[]);
};

#ifndef OPTIONS_CPP
extern options opts;
#endif
