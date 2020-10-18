/*
 * tree.hpp
 *
 *  Created on: Oct 15, 2020
 *      Author: dmarce1
 */

#ifndef ASTROTIGER_TREE_HPP_
#define ASTROTIGER_TREE_HPP_

#include <astrotiger/hpx.hpp>
#include <astrotiger/range.hpp>
#include <astrotiger/tree_client.hpp>
#include "hydro_grid.hpp"
#include "multi_array.hpp"

class tree: public hpx::components::managed_component_base<tree> {
	multi_range box;
	std::vector<std::shared_ptr<grid>> grids;
	std::shared_ptr<hydro_grid> hydro;
	tree_client parent;
	std::vector<tree_client> children;
	std::vector<tree_client> siblings;
	int level;
	double dx;
public:

	static hpx::future<tree_client> allocate(int, multi_range box);

	tree();
	tree(int,multi_range);
	void initialize();
	void set_family(tree_client, std::vector<tree_client>);


	HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,initialize);
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,set_family);


	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & dx;
		arc & box;
		arc & parent;
		arc & children;
		arc & siblings;
		arc & *hydro;
		arc & level;
	}
};

#endif /* ASTROTIGER_TREE_HPP_ */
