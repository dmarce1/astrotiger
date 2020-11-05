#pragma once

#include <astrotiger/multi_array.hpp>
#include <astrotiger/options.hpp>

#include <array>

#define OUTPUT_RESID

#define PACK_POTENTIAL_REDBLACK 1
#define PACK_POTENTIAL 8
#define PACK_ACTIVE 2
#define PACK_SOURCE 4

struct gravity_return {
	double resid;
	std::vector<std::uint8_t> active;
	std::vector<double> R;
	multi_range box;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & resid;
		arc & active;
		arc & R;
	}
};

class gravity {

	multi_array<double> phi0;
	multi_array<double> phi;
	multi_array<double> R;
	multi_array<double> X;
	multi_array<std::uint8_t> active;
#ifdef OUTPUT_RESID
	multi_array<double> resid;
#endif
	double dx;
	int red;
	multi_range box;

public:

	void print() const {
//		for( int i = 0; i < box.volume(); i++) {
//			printf( "%e %e %i\n", X.data()[i], R.data()[i], active.data()[i]);
//		}
	}
	gravity() {
		red = 0;
	}

	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & red;
		arc & phi0;
		arc & phi;
		arc & X;
		arc & R;
		arc & active;
#ifdef OUTPUT_RESID
		arc & resid;
#endif
	}

	void resize(double, const multi_range&);
	void set_avg_zero();
	void store();
	void initialize_fine(const multi_array<double>&, double mtot, int level);
	void initialize_coarse(double w);
	void finish_fine();
	void set_amr_zones(const std::vector<multi_range>&, const std::vector<double>&);
	void zero() {
		for (multi_iterator i(box.pad(-opts.gbw)); !i.end(); i++) {
			X[i] = 0.0;
		}
	}
	double coord(index_type i) const;
	void set_outflow_boundaries();
	multi_array<double> get_phi() const;
	std::vector<double> pack(const multi_range&, double) const;
	std::vector<double> pack(const multi_range&, int) const;
	std::vector<double> pack_phi(const multi_range&) const;
	void unpack(const std::vector<double>&, const multi_range &bbox, int);

	void relax();
	gravity_return get_restrict(double);
	void apply_restrict(const gravity_return&);
	std::vector<double> get_prolong(const multi_range&) const;
	void apply_prolong(const std::vector<double>&);

};
