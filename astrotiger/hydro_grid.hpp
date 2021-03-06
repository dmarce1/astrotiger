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
	std::array<multi_range, NDIM> fbox;
	std::array<multi_range, NDIM> fcbox;
	std::vector<multi_array<double>> U0;
	std::vector<multi_array<double>> U;
	std::array<std::vector<multi_array<double>>, NDIM> F;
	std::array<std::vector<multi_array<double>>, NDIM> Fc;
	multi_array<std::uint8_t> R;
	double dx;

public:

	hydro_grid();
	~hydro_grid();


	void resize(double dx, multi_range box_);

	void sanity(const multi_range& ibox) const {
		if (U0.size()) {
			assert(box == U0[0].box);
			assert(ibox == box.pad(-opts.hbw));
		}
	}
	void reset_flux_registers();
	void reset_coarse_flux_registers();
	double compute_flux(int rk);
	void initialize();
	void compute_refinement_criteria();
	std::vector<multi_range> refined_ranges(const std::vector<multi_range>&,const std::vector<multi_range> &forced) const;
	double coord(index_type i) const;
	std::vector<double> pack_field(int f, multi_range) const;
	std::vector<double> pack(multi_range bbox) const;
	std::vector<std::uint8_t> pack_refinement(multi_range bbox) const;
	std::vector<double> pack_prolong(multi_range bbox, double w) const;
	std::vector<double> pack_restrict(multi_range bbox) const;
	std::vector<double> pack_coarse_flux();
	void unpack_coarse_flux(const std::vector<double>&, const multi_range& bbox, double dt);
	void unpack(const std::vector<double>&, multi_range bbox);
	void unpack_refinement(const std::vector<std::uint8_t>&, multi_range bbox);
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
