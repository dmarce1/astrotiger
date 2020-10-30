#pragma once

#include <astrotiger/multi_array.hpp>
#include <astrotiger/options.hpp>

#include <array>

#define PACK_PHI 1
#define PACK_X 2

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
	multi_array<std::uint8_t> amr;
	multi_array<std::uint8_t> refined;
	std::array<multi_array<double>,NDIM> flux;
	std::array<multi_range,NDIM> fbox;
	multi_array<double> resid;
	multi_array<double> phi_c;
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
	void serialize(A &&arc, unsigned) {
		arc & fbox;
		arc & flux;
		arc & has_phi0;
		arc & has_phi1;
		arc & phi0;
		arc & phi1;
		arc & phi;
		arc & X;
		arc & R;
		arc & active;
		arc & phi_c;
		arc & amr;
		arc & refined;
	}

	void compute_flux();
	void set_refined(const std::vector<multi_range>&);
	void resize(double, const multi_range&);
	void set_avg_zero();
	void store();
	void initialize_fine(const multi_array<double>&, double mtot);
	void initialize_coarse(double w);
	void finish_fine();
	void set_amr_zones(const std::vector<multi_range>&, const std::vector<multi_range>&, const std::vector<double>&);
	void zero() {
		for (multi_iterator i(box.pad(-opts.gbw)); !i.end(); i++) {
			X[i] = 0.0;
		}
	}
	double coord(index_type i) const;
	void set_outflow_boundaries();
	multi_array<double> get_phi() const;
	std::vector<double> pack(const multi_range&, int type) const;
	std::vector<double> pack_amr(const multi_range&, double w) const;
	void unpack(const std::vector<double>&, const multi_range &bbox);

	void relax(bool init_zero);
	gravity_return get_restrict(double);
	void compute_amr_bounds(bool first_pass);
	void apply_restrict(const gravity_return&);
	std::vector<double> get_prolong(const multi_range&) const;
	void apply_prolong(const std::vector<double>&);

};
