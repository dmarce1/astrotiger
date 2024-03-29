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

class energy_statistics;

#define rho_i 0
#define egas_i 1
#define tau_i 2
#define pot_i 3
#define sx_i 4
#define sy_i 5
#define sz_i 6
#define spc_i (sx_i+NDIM)
#define h_i (spc_i+0)
#define hp_i (spc_i+1)
#define he_i (spc_i+2)
#define hep_i (spc_i+3)
#define hepp_i (spc_i+4)

struct statistics;

class hydro_grid {
	multi_range box;
	std::array<multi_range, NDIM> fbox;
	std::array<multi_range, NDIM> fcbox;
	std::vector<multi_array<double>> U0;
	std::vector<multi_array<double>> U;
	std::vector<multi_array<double>> S;
	std::array<std::vector<multi_array<double>>, NDIM> F0;
	std::array<std::vector<multi_array<double>>, NDIM> F;
	std::array<std::vector<multi_array<double>>, NDIM> Fc;
	multi_array<std::uint8_t> R;
	multi_array<double> phi;
	multi_array<double> error;
	double dx;
	double amax;

public:

	hydro_grid();
	~hydro_grid();

	void resize(double dx, multi_range box_);
	void set_error_field(multi_array<double>&&);
	void sanity(const multi_range &ibox) const {
		if (U0.size()) {
			assert(box == U0[0].box);
			assert(ibox == box.pad(-opts.hbw));
		}
	}
	void set_phi(multi_array<double> &&phi_) {
		phi = std::move(phi_);
	}
	const multi_array<double>& get_density() const {
		return U[rho_i];
	}
	double positivity_limit() const;
	void store();
	void store_flux();
	void enforce_physical_bc(int);
	void reset_flux_registers();
	void reset_coarse_flux_registers();
	double compute_flux(int rk);
	void initialize();
	void compute_refinement_criteria(const multi_array<int> &pcount);
	std::vector<multi_range> refined_ranges(const std::vector<multi_range>&, const std::vector<multi_range> &forced) const;
	double coord(index_type i) const;
	std::vector<std::vector<double>> pack_output(const multi_array<std::uint8_t>&) const;
	std::vector<double> pack_field(int f, multi_range) const;
	std::vector<double> pack(multi_range bbox) const;
	std::vector<std::uint8_t> pack_refinement(multi_range bbox) const;
	std::vector<double> pack_prolong(multi_range bbox, double w) const;
	std::vector<double> pack_field_prolong(int f, multi_range bbox, double w) const;
	std::vector<double> pack_restrict(multi_range bbox) const;
	std::vector<double> pack_coarse_correction();
	void unpack_coarse_correction(const std::vector<double>&, const multi_range &bbox, double dt, double a0, double a1);
	void unpack(const std::vector<double>&, multi_range bbox);
	void unpack_field(int f, const std::vector<double>&, multi_range bbox);
	void unpack_refinement(const std::vector<std::uint8_t>&, multi_range bbox);
	void substep_update(int rk, double dt, double a0, double a1);
	void update_energy();
	statistics get_statistics(const std::vector<multi_range>&) const;
	double compare_analytic(const std::vector<multi_range> &cboxes, multi_array<double> &results) const;
	static std::vector<std::string> field_names();
	void to_array(multi_array<double>&, const multi_range &bbox, int, double) const;
	std::vector<double> pack_coarse_flux();
	double unpack_fine_flux(const std::vector<double>&, const multi_range&);
	energy_statistics get_energy_statistics(const multi_array<double> &phi, const std::vector<multi_range> &exclude, double rho0) const;

	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & F;
		arc & F0;
		arc & fbox;
		arc & box;
		arc & U0;
		arc & U;
		arc & R;
		arc & dx;
	}
};

#endif /* ASTROTIGER_HYDRO_GRID_HPP_ */
