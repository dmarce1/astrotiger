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

#define SILO_CHECK(b) \
	if( b != 0 ) { \
		printf( "SILO call failed in %s on line %i\n", __FILE__, __LINE__); \
	}

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
	std::atomic<int> step;
	double dx;
	double t0;
	double t;
public:

	static hpx::future<tree_client> allocate(int, multi_range box);

	tree();
	~tree();
	tree(int, multi_range);
	double initialize(int);
	std::vector<tree_client> get_children() const;
	void set_family(tree_client, tree_client, std::vector<sibling>);
	void clear_family();
	void set_child_family();
	std::vector<double> get_hydro_boundary(multi_range, int);
	std::vector<std::vector<double>> get_hydro_prolong(std::vector<multi_range>, double);
	std::vector<double> get_hydro_restrict();
	void hydro_substep(int, double);
	double hydro_initialize();
	void output(DBfile* db) const;

	void get_hydro_boundaries(bool amr);

	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_hydro_boundary);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_hydro_prolong);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_hydro_restrict);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_children);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,initialize);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,set_family);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,clear_family);

	HPX_SERIALIZATION_SPLIT_MEMBER();
	template<class A>
	void save(A &&arc, unsigned) const {
		arc & t;
		arc & t0;
		arc & level;
		arc & dx;
		arc & box;
		arc & parent;
		arc & self;
		arc & children;
		arc & siblings;
		arc & hydro;
	}
	template<class A>
	void load(A &&arc, unsigned) {
		if (level >= 0) {
			levels_remove_entry(level, this);
		}
		arc & level;
		levels_add_entry(level, this);
		arc & dx;
		arc & box;
		arc & parent;
		arc & children;
		arc & siblings;
		arc & hydro;
	}
};

#endif /* ASTROTIGER_TREE_HPP_ */
