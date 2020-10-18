/*
 * tree_client.cpp
 *
 *  Created on: Oct 15, 2020
 *      Author: dmarce1
 */

#include <astrotiger/tree.hpp>

hpx::future<void> tree_client::set_family(tree_client p, std::vector<tree_client> c) const {
	return hpx::async<tree::set_family_action>(gid, p, std::move(c));
}

hpx::future<void> tree_client::initialize() const {
	return hpx::async<tree::initialize_action>(gid);
}

