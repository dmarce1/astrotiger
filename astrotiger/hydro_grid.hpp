/*
 * grid.hpp
 *
 *  Created on: Oct 14, 2020
 *      Author: dmarce1
 */

#ifndef ASTROTIGER_HYDRO_GRID_HPP_
#define ASTROTIGER_HYDRO_GRID_HPP_

#include <astrotiger/hpx.hpp>
#include <astrotiger/range.hpp>
#include <astrotiger/multi_array.hpp>
#include <astrotiger/options.hpp>

#include <array>

#define rho_i 0
#define egas_i 1
#define sx_i 2
#define sy_i 3
#define sz_i 4

class hydro_grid {
	multi_range box;
	std::array<multi_range,NDIM> fbox;
	std::vector<multi_array<double>> U0;
	std::vector<multi_array<double>> U;
	std::array<std::vector<multi_array<double>>, NDIM> F;
	multi_array<std::uint8_t> R;
	double dx;
public:

	hydro_grid();
	~hydro_grid();

	void resize(double dx, multi_range box_);

	double compute_flux();
	void initialize();
	void compute_refinement_criteria();
	std::vector<multi_range> refined_ranges() const;
	double coord(index_type i) const;
	std::vector<double> pack_field(int f) const;
	std::vector<double> pack_boundary(multi_range bbox) const;
	std::vector<double> pack_prolong(multi_range bbox, double w) const;
	std::vector<double> pack_restrict(multi_range bbox) const;
	void unpack(const std::vector<double>&, multi_range bbox);
	void substep_update(int rk, double dt);

	static std::vector<std::string> field_names();

	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & fbox;
		arc & box;
		arc & U0;
		arc & U;
		arc & R;
		arc & dx;
	}
};

#endif /* ASTROTIGER_HYDRO_GRID_HPP_ */
