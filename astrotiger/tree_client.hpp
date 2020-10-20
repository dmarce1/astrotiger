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
class sibling;

class tree_client {
	hpx::id_type gid;
	multi_range box;

public:
	bool operator==(const tree_client &other) const {
		return gid == other.gid;
	}

	bool operator!=(const tree_client &other) const {
		return gid != other.gid;
	}

	tree_client() {
		gid = hpx::id_type();
		box.set_null();
	}

	multi_range get_box() const {
		return box;
	}

	tree_client(const hpx::id_type id_, multi_range box_) {
		gid = id_;
		box = box_;
	}

	hpx::future<void> clear_family() const;
	hpx::future<void> set_family(tree_client p, tree_client, std::vector<sibling> c) const;
	hpx::future<double> initialize(int) const;
	hpx::future<std::vector<tree_client>> get_children() const;
	hpx::future<std::vector<double>> get_hydro_boundary(multi_range b, int) const;
	hpx::future<std::vector<std::vector<double>>> get_hydro_prolong(std::vector<multi_range> b, double t) const;
	hpx::future<std::vector<double>> get_hydro_restrict() const;
	hpx::future<std::vector<multi_range>> get_child_boxes() const;


	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & gid;
		arc & box;
	}
};

#endif /* ASTROTIGER_TREE_CLIENT_HPP_ */
