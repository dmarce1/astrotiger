/*
 * tree.hpp
 *
 *  Created on: Oct 15, 2020
 *      Author: dmarce1
 */

#ifndef ASTROTIGER_TREE_HPP_
#define ASTROTIGER_TREE_HPP_

#include <astrotiger/hpx.hpp>
#include <astrotiger/levels.hpp>
#include <astrotiger/range.hpp>
#include <astrotiger/tree_client.hpp>
#include <astrotiger/hydro_grid.hpp>
#include <astrotiger/gravity.hpp>
#include <astrotiger/particles.hpp>

#include <silo.h>

#define GRAVITY_FINAL_PASS 1000000


struct statistics {
	std::vector<double> u;
	int min_level;
	int max_level;

	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & u;
		arc & min_level;
		arc & max_level;
	}
};


struct sibling {
	tree_client client;
	vect<int> shift;
	multi_range box() const {
		return client.get_box().shift(shift);
	}
	template<class A>
	void serialize(A &arc, unsigned) {
		arc & client;
		arc & shift;
	}
};

struct output_return {
	std::vector<std::vector<double>> coords;
	std::vector<std::vector<double>> data;
	std::vector<std::vector<double>> pcoords;
	std::vector<std::vector<double>> pdata;
	std::vector<int> zones;
};

class tree: public hpx::components::managed_component_base<tree> {
	multi_range box;
	hydro_grid hydro;
	gravity grav;
	particles parts;
	tree_client parent;
	tree_client self;
	std::vector<tree_client> children;
	std::vector<sibling> siblings;
	int level;
	std::vector<particle> new_parts;
	std::atomic<int> refine_step;
	std::atomic<int> energy_step;
	std::vector<multi_range> grandchild_boxes;
	mutex_type mtx;
	int gravity_step;
	int hydro_step;
	double dx;
	double t0;
	double t;
	double dt;
	double last_dt;
public:

	static hpx::future<tree_client> allocate(int, multi_range box);

	tree(tree &&other) :
			refine_step(0), energy_step(0), gravity_step(0), hydro_step(0) {
		dt = std::move(other.dt);
		box = std::move(other.box);
		hydro = std::move(other.hydro);
		grav = std::move(other.grav);
		parent = std::move(other.parent);
		self = std::move(other.self);
		children = std::move(other.children);
		siblings = std::move(other.siblings);
		level = std::move(other.level);
		grandchild_boxes = std::move(other.grandchild_boxes);
		dx = std::move(other.dx);
		t0 = std::move(other.t0);
		t = std::move(other.t);
		auto tmp = (int) other.refine_step;
		refine_step = tmp;
		last_dt = std::move(other.last_dt);
		parts = std::move(other.parts);
//		printf("Adding entry %i\n", level);
	}
	tree(const tree &other) :
			refine_step(0), energy_step(0), gravity_step(0),hydro_step(0) {
		dt = other.dt;
		last_dt = other.last_dt;
		box = other.box;
		hydro = other.hydro;
		grav = other.grav;
		parent = other.parent;
		self = other.self;
		children = other.children;
		siblings = other.siblings;
		level = other.level;
		grandchild_boxes = other.grandchild_boxes;
		dx = other.dx;
		t0 = other.t0;
		t = other.t;
		auto tmp = (int) other.refine_step;
		refine_step = tmp;
		tmp = (int) other.energy_step;
		energy_step = tmp;
		parts = other.parts;
//		printf("Adding entry %i\n", level);
	}
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & dt;
		arc & last_dt;
		arc & box;
		arc & hydro;
		arc & grav;
		arc & parent;
		arc & self;
		arc & children;
		arc & siblings;
		arc & level;
		arc & grandchild_boxes;
		arc & dx;
		arc & t0;
		arc & t;
		int tmp = (int) refine_step;
		arc & tmp;
		refine_step = tmp;
		tmp = (int) energy_step;
		arc & tmp;
		energy_step = tmp;
		arc & parts;
	}

	tree();
	~tree();
	tree(int, multi_range);
	std::shared_ptr<tree> get_ptr();
	double initialize(int);
	std::vector<multi_range> get_amr_boxes() const;
	std::vector<tree_client> get_children() const;
	void set_family(tree_client, tree_client, std::vector<sibling>);
	void delist();
	void list();
	void set_child_family();
	void set_gravity_boundary(boundary&&, const multi_range&, int);
	void set_boundary(std::vector<double>&&, const multi_range&, int);
	std::vector<double> get_energy_boundary(multi_range, int);
	std::pair<std::vector<std::uint8_t>, std::vector<multi_range>> get_refinement_boundary(multi_range, int);
	std::vector<std::vector<double>> get_hydro_prolong(std::vector<multi_range>, double);
	std::vector<std::vector<double>> get_energy_prolong(std::vector<multi_range>, double);
	std::pair<std::vector<double>, std::vector<double>> get_hydro_restrict();
	void hydro_substep(int, double, bool, double, double);
	void hydro_initialize(bool);
	output_return output(DBfile *db, int) const;
	multi_range get_box() const;
	tree_client truncate(tree_client, multi_range box);
	std::vector<multi_range> get_refinement_boundaries();
	void get_gravity_boundaries(int);
	void get_hydro_boundaries(double);
	void get_energy_boundaries(double);
	void sanity() const;
	statistics get_statistics(int lev, double t);
	std::vector<double> restrict_all();
	gravity_return gravity_solve(int pass, int level, const std::vector<double> coarse, double t, double m, int iterd);
	double compute_error();
	std::vector<double> get_fine_flux();
	double apply_fine_fluxes();
	void drift(double);
	void finish_drift(std::vector<particle>);
	void recv_parts(std::vector<particle>);
	range<double> range_int_to_double( const multi_range& box );
	multi_array<int> get_particle_count() const;
	hpx::future<void> send_child_particles();
	double max_part_velocity() const;
	void kick(int, double, std::vector<double>, std::vector<double>);
	std::vector<double> compute_cic(const std::vector<double>&, double t, int );
	void apply_coarse_correction(double a0, double a1);


	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,compute_cic);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,kick);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,max_part_velocity);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_particle_count);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,recv_parts);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,finish_drift);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,drift);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_fine_flux);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,restrict_all);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_statistics);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_refinement_boundary);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_box);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,truncate);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_ptr);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_energy_boundary);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,set_boundary);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_energy_prolong);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_hydro_prolong);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_hydro_restrict);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_children);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,set_family);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,list);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,delist);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,gravity_solve);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,initialize);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,compute_error);

};

#endif /* ASTROTIGER_TREE_HPP_ */
