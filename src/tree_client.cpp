/*
 * tree_client.cpp
 *
 *  Created on: Oct 15, 2020
 *      Author: dmarce1
 */

#include <astrotiger/tree.hpp>

hpx::future<std::vector<double>> tree_client::compute_cic(const std::vector<double> &data, double t, int lev) const {
	return hpx::async<tree::compute_cic_action>(gid, data, t, lev);
}

hpx::future<void> tree_client::kick(int rung, double t, const std::vector<double> &dt0, const std::vector<double> &dt1) const {
	return hpx::async<tree::kick_action>(gid, rung, t, dt0, dt1);
}

hpx::future<multi_array<int>> tree_client::get_particle_count() const {
	return hpx::async<tree::get_particle_count_action>(gid);
}

hpx::future<void> tree_client::delist() const {
	assert(gid != hpx::invalid_id);
	return hpx::async<tree::delist_action>(gid);
}

hpx::future<void> tree_client::send_parts(std::vector<particle> &&parts) const {
	return hpx::async<tree::recv_parts_action>(gid, std::move(parts));
}

hpx::future<double> tree_client::max_part_velocity() const {
	return hpx::async<tree::max_part_velocity_action>(gid);
}

hpx::future<void> tree_client::drift(double dt) const {
	return hpx::async<tree::drift_action>(gid, dt);
}

hpx::future<void> tree_client::finish_drift(std::vector<particle> &&parts) const {
	return hpx::async<tree::finish_drift_action>(gid, std::move(parts));
}

hpx::future<gravity_return> tree_client::gravity_solve(int pass, int level, std::vector<double> &&coarse, double t, double m) const {
	return hpx::async<tree::gravity_solve_action>(gid, pass, level, std::move(coarse), t, m);
}

hpx::future<statistics> tree_client::get_statistics(int level, double t) const {
	return hpx::async<tree::get_statistics_action>(gid, level, t);
}

hpx::future<void> tree_client::list() const {
	assert(gid != hpx::invalid_id);
	return hpx::async<tree::list_action>(gid);
}

hpx::future<void> tree_client::set_family(tree_client p, tree_client s, std::vector<sibling> c) const {
	assert(gid != hpx::invalid_id);
	return hpx::async<tree::set_family_action>(gid, p, s, std::move(c));
}

hpx::future<double> tree_client::initialize(int lev) const {
	assert(gid != hpx::invalid_id);
	return hpx::async<tree::initialize_action>(gid, lev);
}

hpx::future<std::vector<tree_client>> tree_client::get_children() const {
	assert(gid != hpx::invalid_id);
	return hpx::async<tree::get_children_action>(gid);
}

hpx::future<std::vector<double>> tree_client::restrict_all() const {
	return hpx::async<tree::restrict_all_action>(gid);
}

hpx::future<std::vector<double>> tree_client::get_fine_flux() const {
	return hpx::async<tree::get_fine_flux_action>(gid);
}

hpx::future<void> tree_client::set_boundary(std::vector<double> &&data, const multi_range &id, int step) const {
	assert(gid != hpx::invalid_id);
	return hpx::async<tree::set_boundary_action>(gid, std::move(data), id, step);
}

hpx::future<std::vector<double>> tree_client::get_energy_boundary(multi_range b, int step) const {
	assert(gid != hpx::invalid_id);
	return hpx::async<tree::get_energy_boundary_action>(gid, b, step);
}

hpx::future<double> tree_client::compute_error() const {
	return hpx::async<tree::compute_error_action>(gid);
}

hpx::future<std::pair<std::vector<std::uint8_t>, std::vector<multi_range>>> tree_client::get_refinement_boundary(multi_range b, int step) const {
	assert(gid != hpx::invalid_id);
	return hpx::async<tree::get_refinement_boundary_action>(gid, b, step);
}

hpx::future<std::vector<std::vector<double>>> tree_client::get_energy_prolong(std::vector<multi_range> b, double t) const {
	assert(gid != hpx::invalid_id);
	return hpx::async<tree::get_hydro_prolong_action>(gid, std::move(b), t);
}

hpx::future<std::vector<std::vector<double>>> tree_client::get_hydro_prolong(std::vector<multi_range> b, double t) const {
	assert(gid != hpx::invalid_id);
	return hpx::async<tree::get_hydro_prolong_action>(gid, std::move(b), t);
}

hpx::future<std::pair<std::vector<double>, std::vector<double>>> tree_client::get_hydro_restrict() const {
	assert(gid != hpx::invalid_id);
	return hpx::async<tree::get_hydro_restrict_action>(gid);
}

hpx::future<tree_client> tree_client::truncate(multi_range box_) const {
	assert(gid != hpx::invalid_id);
	return hpx::async<tree::truncate_action>(gid, *this, box_);
}

hpx::future<std::shared_ptr<tree>> tree_client::get_ptr() const {
	assert(gid != hpx::invalid_id);
	return hpx::async<tree::get_ptr_action>(gid);
}
