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


hpx::future<double> tree_client::initialize(int lev) const {
	return hpx::async<tree::initialize_action>(gid, lev);
}

hpx::future<std::vector<tree_client>> tree_client::get_children() const {
	return hpx::async<tree::get_children_action>(gid);
}

hpx::future<std::vector<double>> tree_client::get_hydro_boundary(multi_range b, int step) {
	return hpx::async<tree::get_hydro_boundary_action>(gid, b, step);
}

hpx::future<std::vector<std::vector<double>>> tree_client::get_hydro_prolong(std::vector<multi_range> b, double t) {
	return hpx::async<tree::get_hydro_prolong_action>(gid, std::move(b), t);
}

hpx::future<std::vector<double>> tree_client::get_hydro_restrict(multi_range b) {
	return hpx::async<tree::get_hydro_restrict_action>(gid, b);
}
