/*
 * particles.hpp
 *
 *  Created on: Nov 7, 2020
 *      Author: dmarce1
 */

#ifndef ASTROTIGER_PARTICLES_HPP_
#define ASTROTIGER_PARTICLES_HPP_

#include <astrotiger/defs.hpp>
#include <astrotiger/multi_array.hpp>
#include <astrotiger/vect.hpp>

#include <cstdint>
#include <vector>
#include <array>


class energy_statistics;

struct particle {
	vect<double> x;
	vect<float> v;
	float m;
	int rung;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & x;
		arc & v;
		arc & m;
		arc & rung;
	}
};

class particles {
	std::vector<particle> parts;
	multi_array<double> rho;
	multi_range box;
	range<double> rbox;
	std::vector<range<double>> child_boxes;
	double dx;
public:

	multi_array<double> get_cic() const;
	multi_array<int> particle_count() const;
	std::vector<double> pack_cic_restrict() const;
	std::vector<double> pack_cic(const multi_range&) const;
	void unpack_cic(const std::vector<double>&, const multi_range&);
	void unpack_cic_prolong(const std::vector<double>&, const multi_range&);
	void compute_cloud_in_cell(double dt);
	std::vector<std::vector<double>> pack_output() const;
	std::vector<std::vector<double>> pack_coords() const;
	static std::vector<std::string> field_names();
	std::vector<particle> drift(double dt, double a0, double a1);
	void add_parts(const std::vector<particle> &parts);
	void initialize();
	std::vector<particle> get_particles();
	std::vector<particle> get_particles(const multi_range&);
	double max_velocity() const;
	void resize(double dx, const multi_range&);
	void zero_cic(const multi_range&);
	void set_child_boxes(const std::vector<multi_range>&);
	std::vector<particle> get_child_parts();
	void kick(int kick_level, int this_level, const std::vector<double>&, const std::vector<double>&, const std::array<multi_array<double>, NDIM> &g);
	energy_statistics get_energy_statistics(const multi_array<double>& phi) const;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & parts;
		arc & box;
		arc & rbox;
		arc & dx;
	}
};

#endif /* ASTROTIGER_PARTICLES_HPP_ */
