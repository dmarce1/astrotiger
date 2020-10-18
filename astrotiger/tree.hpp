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
#include "hydro_grid.hpp"
#include "multi_array.hpp"

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
	std::vector<std::shared_ptr<grid>> grids;
	std::shared_ptr<hydro_grid> hydro;
	tree_client parent;
	tree_client self;
	std::vector<tree_client> children;
	std::vector<sibling> siblings;
	int level;
	double dx;
public:

	static hpx::future<tree_client> allocate(int, multi_range box);

	tree();
	~tree();
	tree(int, multi_range);
	void initialize();
	std::vector<tree_client> get_children() const;
	void set_family(tree_client, tree_client, std::vector<sibling>);
	void clear_family();
	void set_child_family();

	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_children);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,initialize);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,set_family);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,clear_family);

	HPX_SERIALIZATION_SPLIT_MEMBER();
	template<class A>
	void save(A &&arc, unsigned) const {
		arc & level;
		arc & dx;
		arc & box;
		arc & parent;
		arc & self;
		arc & children;
		arc & siblings;
		arc & *hydro;
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
		arc & *hydro;
	}
};

#endif /* ASTROTIGER_TREE_HPP_ */
