/*
 * tree_client.hpp
 *
 *  Created on: Oct 15, 2020
 *      Author: dmarce1
 */

#ifndef ASTROTIGER_TREE_CLIENT_HPP_
#define ASTROTIGER_TREE_CLIENT_HPP_

#include <astrotiger/hpx.hpp>
#include <astrotiger/multi_array.hpp>

class tree;

class tree_client {
	hpx::id_type gid;
	multi_range box;

public:
	tree_client() {
		box.set_null();
	}

	tree_client(const hpx::id_type id_, multi_range box_) {
		gid = id_;
		box = box_;
	}

	hpx::future<void> set_family(tree_client p, std::vector<tree_client> c) const;
	hpx::future<void> initialize() const;

	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & gid;
		arc & box;
	}
};

#endif /* ASTROTIGER_TREE_CLIENT_HPP_ */
