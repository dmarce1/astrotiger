#pragma once

#include <astrotiger/multi_array.hpp>
#include <astrotiger/options.hpp>

#include <array>


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
	multi_array<double> phi1;
	multi_array<double> phi;
	multi_array<double> X;
	multi_array<double> R;
	multi_array<std::uint8_t> active;
	double dx;
	bool has_phi0;
	bool has_phi1;
	multi_range box;

public:

	void print() const {
//		for( int i = 0; i < box.volume(); i++) {
//			printf( "%e %e %i\n", X.data()[i], R.data()[i], active.data()[i]);
//		}
	}
	gravity() {
		has_phi0 = false;
		has_phi1 = false;
	}

	template<class A>
	void serialize(A&& arc, unsigned ) {
		arc & has_phi0;
		arc & has_phi1;
		arc & phi0;
		arc & phi1;
		arc & phi;
		arc & X;
		arc & R;
		arc & active;
	}

	void resize(double, const multi_range&);
	void set_avg_zero();
	void store();
	void initialize_fine(const multi_array<double>&, double mtot);
	void initialize_coarse(double w);
	void finish_fine();

	void zero() {
		for( multi_iterator i(box.pad(-opts.gbw));!i.end(); i++) {
			X[i] = 0.0;
		}
	}

	multi_array<double> get_phi() const;
	std::vector<double> pack(const multi_range&) const;
	void unpack(const std::vector<double>&, const multi_range &bbox);
	std::vector<double> pack_prolong_amr(const multi_range&) const;

	void relax(bool init_zero);
	gravity_return get_restrict() ;
	void apply_restrict(const gravity_return&);
	std::vector<double> get_prolong(const multi_range&) const;
	void apply_prolong(const std::vector<double>&);

};
