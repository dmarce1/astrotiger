/*
 * tree_client.cpp
 *
 *  Created on: Oct 15, 2020
 *      Author: dmarce1
 */

#include <astrotiger/tree.hpp>

hpx::future<void> tree_client::clear_family() const {
	return hpx::async<tree::clear_family_action>(gid);
}

hpx::future<void> tree_client::set_family(tree_client p, tree_client s, std::vector<sibling> c) const {
	return hpx::async<tree::set_family_action>(gid, p, s, std::move(c));
}

hpx::future<void> tree_client::initialize() const {
	return hpx::async<tree::initialize_action>(gid);
}

hpx::future<std::vector<tree_client>> tree_client::get_children() const {
	return hpx::async<tree::get_children_action>(gid);
}
