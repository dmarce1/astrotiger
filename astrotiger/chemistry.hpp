/*
 * chemsitry.hpp
 *
 *  Created on: Dec 2, 2020
 *      Author: dmarce1
 */

#ifndef ASTROTIGER_CHEMISTRY_HPP_
#define ASTROTIGER_CHEMISTRY_HPP_

#define amu 1.6605e-24

struct species {
	double H;
	double Hp;
	double Hn;
	double H2;
	double H2p;
	double He;
	double Hep;
	double Hepp;
};

struct thermo_props {
	double pressure;
	double T;
	double ne;
	double sound_speed;
	double gamma;
	double rho;
	double eion;
	double cv;
};

double ion_energy(species s);
thermo_props compute_thermo_properties(const species s, double energy);
species compute_next_species(const species s0, double energy, double Trad, double dt);

#endif /* ASTROTIGER_CHEMISTRY_HPP_ */
