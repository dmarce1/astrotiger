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

struct particle {
	double next_kick;
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
	multi_range box;
	range<double> rbox;
	double dx;
public:

	multi_array<double> cloud_in_cell() const;
	std::vector<std::vector<double>> pack_output() const;
	std::vector<std::vector<double>> pack_coords() const;
	static std::vector<std::string> field_names();
	std::vector<particle> drift(double dt);
	void add_parts(const std::vector<particle> &parts);
	void initialize();
	double max_velocity() const;
	void resize(double dx, const multi_range&);

	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & parts;
		arc & box;
		arc & rbox;
		arc & dx;
	}
};

#endif /* ASTROTIGER_PARTICLES_HPP_ */
