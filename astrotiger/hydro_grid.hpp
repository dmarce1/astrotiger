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
#include <astrotiger/grid.hpp>
#include <astrotiger/multi_array.hpp>
#include <astrotiger/options.hpp>

#include <array>

#define rho_i 0
#define egas_i 1
#define sx_i 2
#define sy_i 3
#define sz_i 4

class hydro_grid: public grid {
	multi_range box;
	std::vector<multi_array<double>> U0;
	std::vector<multi_array<double>> U1;
	multi_array<std::uint8_t> R;
	double dx;
public:

	hydro_grid();
	~hydro_grid();

	void resize(double dx, multi_range box_);

	void initialize();
	void compute_refinement_criteria();
	std::vector<multi_range> refined_ranges() const;
	double coord(index_type i) const;
	std::vector<double> pack_boundary(multi_range bbox);
	void unpack_boundary(const std::vector<double>&,  multi_range bbox);

	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & box;
		arc & U0;
		arc & U1;
		arc & R;
		arc & dx;
	}
};

#endif /* ASTROTIGER_HYDRO_GRID_HPP_ */
