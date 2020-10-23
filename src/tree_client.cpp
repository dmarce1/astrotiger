/*
 * tree_client.cpp
 *
 *  Created on: Oct 15, 2020
 *      Author: dmarce1
 */

#include <astrotiger/tree.hpp>

hpx::future<void> tree_client::delist() const {
	assert(gid!=hpx::invalid_id);
	return hpx::async<tree::delist_action>(gid);
}

hpx::future<void> tree_client::list() const {
	assert(gid!=hpx::invalid_id);
	return hpx::async<tree::list_action>(gid);
}

hpx::future<void> tree_client::set_family(tree_client p, tree_client s, std::vector<sibling> c) const {
	assert(gid!=hpx::invalid_id);
	return hpx::async<tree::set_family_action>(gid, p, s, std::move(c));
}


hpx::future<double> tree_client::initialize(int lev) const {
	assert(gid!=hpx::invalid_id);
	return hpx::async<tree::initialize_action>(gid, lev);
}

hpx::future<std::vector<tree_client>> tree_client::get_children() const {
	assert(gid!=hpx::invalid_id);
	return hpx::async<tree::get_children_action>(gid);
}

hpx::future<std::vector<double>> tree_client::get_hydro_boundary(multi_range b, int step) const {
	assert(gid!=hpx::invalid_id);
	return hpx::async<tree::get_hydro_boundary_action>(gid, b, step);
}

hpx::future<std::pair<std::vector<std::uint8_t>,std::vector<multi_range>>> tree_client::get_refinement_boundary(multi_range b, int step) const {
	assert(gid!=hpx::invalid_id);
	return hpx::async<tree::get_refinement_boundary_action>(gid, b, step);
}

hpx::future<std::vector<std::vector<double>>> tree_client::get_hydro_prolong(std::vector<multi_range> b, double t) const {
	assert(gid!=hpx::invalid_id);
	return hpx::async<tree::get_hydro_prolong_action>(gid, std::move(b), t);
}

hpx::future<std::pair<std::vector<double>,std::vector<double>>> tree_client::get_hydro_restrict() const {
	assert(gid!=hpx::invalid_id);
	return hpx::async<tree::get_hydro_restrict_action>(gid);
}


hpx::future<tree_client> tree_client::truncate( multi_range box_) const {
	assert(gid!=hpx::invalid_id);
	return hpx::async<tree::truncate_action>(gid, *this, box_);
}


hpx::future<std::shared_ptr<tree>> tree_client::get_ptr() const {
	assert(gid!=hpx::invalid_id);
	return hpx::async<tree::get_ptr_action>(gid);
}
