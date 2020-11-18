#pragma once
#include <astrotiger/defs.hpp>

#include <string>
#include <vector>
#include <array>

using bc_type = int;
#define REFLECTING 0
#define PERIODIC 1
#define OUTFLOW 2

class options {
public:
	double code_to_s;
	double H0;
	double clight;
	double m_tot;
	double omega_m;
	double omega_b;
	double code_to_cm;
	double code_to_cm_per_s;
	double code_to_g;
	std::string config_file;
	int nhydro;
	int ngroup;
	int max_box;
	int min_box;
	int hbw;
	int gbw;
	int max_bw;
	int window;
	int max_level;
	int order;
	int nrk;
	int nmulti;
	bool self_gravity;
	bool gravity;
	bool hydro;
	bool particles;
	std::string problem;
	double refine_slope;
	double gamma;
	double cfl;
	double G;
	double tmax;
	double efficiency;
	double part_mass;
	double z;
	std::string xpbnd;
	std::string ypbnd;
	std::string zpbnd;
	std::string xmbnd;
	std::string ymbnd;
	std::string zmbnd;
	std::array<bc_type, 2 * NDIM> bnd;

	std::vector<double> alpha;
	std::vector<double> beta;

	template<class Arc>
	void serialize(Arc &arc, unsigned) {
		arc & z;
		arc & part_mass;
		arc & omega_b;
		arc & code_to_g;
		arc & code_to_cm_per_s;
		arc & code_to_cm;
		arc & code_to_s;
		arc & H0;
		arc & clight;
		arc & m_tot;
		arc & omega_m;
		arc & particles;
		arc & self_gravity;
		arc & gravity;
		arc & hydro;
		arc & bnd;
		arc & xpbnd;
		arc & ypbnd;
		arc & zpbnd;
		arc & xmbnd;
		arc & ymbnd;
		arc & zmbnd;
		arc & nmulti;
		arc & gbw;
		arc & G;
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
