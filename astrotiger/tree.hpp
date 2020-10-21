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

#include <silo.h>

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

class tree: public hpx::components::managed_component_base<tree> {
	multi_range box;
	hydro_grid hydro;
	tree_client parent;
	tree_client self;
	std::vector<tree_client> children;
	std::vector<sibling> siblings;
	int level;
	std::atomic<int> hydro_step;
	std::atomic<int> refine_step;
	std::vector<multi_range> grandchild_boxes;
	double dx;
	double t0;
	double t;
public:

	static hpx::future<tree_client> allocate(int, multi_range box);

	tree(tree &&other) {
		box = std::move(other.box);
		hydro = std::move(other.hydro);
		parent = std::move(other.parent);
		self = std::move(other.self);
		children = std::move(other.children);
		siblings = std::move(other.siblings);
		level = std::move(other.level);
		grandchild_boxes = std::move(other.grandchild_boxes);
		dx = std::move(other.dx);
		t0 = std::move(other.t0);
		t = std::move(other.t);
		levels_add_entry(level, this);
	}
	tree(const tree &other) {
		box = other.box;
		hydro = other.hydro;
		parent = other.parent;
		self = other.self;
		children = other.children;
		siblings = other.siblings;
		level = other.level;
		grandchild_boxes = other.grandchild_boxes;
		dx = other.dx;
		t0 = other.t0;
		t = other.t;
		levels_add_entry(level, this);
	}
	HPX_SERIALIZATION_SPLIT_MEMBER ();
	template<class A>
	void load(A &&arc, unsigned) {
		arc & box;
		arc & hydro;
		arc & parent;
		arc & self;
		arc & children;
		arc & siblings;
		arc & level;
		arc & grandchild_boxes;
		arc & dx;
		arc & t0;
		arc & t;
		levels_add_entry(level, this);
	}
	template<class A>
	void save(A &&arc, unsigned) const {
		arc & box;
		arc & hydro;
		arc & parent;
		arc & self;
		arc & children;
		arc & siblings;
		arc & level;
		arc & grandchild_boxes;
		arc & dx;
		arc & t0;
		arc & t;
	}

	tree();
	~tree();
	tree(int, multi_range);
	std::shared_ptr<tree> get_ptr();
	double initialize(int);
	std::vector<multi_range> get_amr_boxes() const;
	std::vector<tree_client> get_children() const;
	void set_family(tree_client, tree_client, std::vector<sibling>);
	void clear_family();
	void set_child_family();
	std::vector<double> get_hydro_boundary(multi_range, int);
	std::vector<std::vector<double>> get_hydro_prolong(std::vector<multi_range>, double);
	std::vector<double> get_hydro_restrict();
	void hydro_substep(int, double);
	double hydro_initialize(bool);
	std::string output(DBfile *db) const;
	std::vector<multi_range> get_child_boxes() const;
	std::vector<multi_range> get_grandchild_boxes(int) const;

	void get_hydro_boundaries(bool amr);

	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_ptr);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_grandchild_boxes);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_child_boxes);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_hydro_boundary);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_hydro_prolong);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_hydro_restrict);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_children);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,initialize);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,set_family);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,clear_family);

};

#endif /* ASTROTIGER_TREE_HPP_ */
