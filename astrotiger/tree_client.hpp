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
#include <astrotiger/hydro_grid.hpp>
#include <astrotiger/channel.hpp>

class tree;
class sibling;
class gravity_return;

class tree_client {
	hpx::id_type gid;
	multi_range box;
	channel<std::vector<double>> chan;

public:
	hpx::id_type get_id() {
		return gid;
	}
	tree_client(tree_client &&other) {
		*this = std::move(other);
	}
	tree_client(const tree_client &other) {
		*this = other;
	}
	tree_client& operator=(tree_client &&other) {
		gid = std::move(other.gid);
		box = std::move(other.box);
		return *this;
	}
	tree_client& operator=(const tree_client &other) {
		gid = other.gid;
		box = other.box;
		return *this;
	}
	bool operator==(const tree_client &other) const {
		return gid == other.gid;
	}

	bool operator!=(const tree_client &other) const {
		return gid != other.gid;
	}

	void put(std::vector<double> &&data) {
		chan.put(std::move(data));
	}

	hpx::future<std::vector<double>> get() {
		return hpx::async([this]() {
			return chan.get();
		});
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

	hpx::future<tree_client> truncate(multi_range box) const;
	hpx::future<void> delist() const;
	hpx::future<void> list() const;
	hpx::future<void> set_family(tree_client p, tree_client, std::vector<sibling> c) const;
	hpx::future<double> initialize(int) const;
	hpx::future<std::vector<tree_client>> get_children() const;
	hpx::future<void> set_boundary(std::vector<double>&& data, const multi_range& bbox) const;
	hpx::future<std::vector<double>> get_energy_boundary(multi_range b, int) const;
	hpx::future<std::vector<std::vector<double>>> get_hydro_prolong(std::vector<multi_range> b, double t) const;
	hpx::future<std::vector<std::vector<double>>> get_energy_prolong(std::vector<multi_range> b, double t) const;
	hpx::future<std::pair<std::vector<double>, std::vector<double>>> get_hydro_restrict() const;
	hpx::future<std::shared_ptr<tree>> get_ptr() const;
	hpx::future<std::pair<std::vector<std::uint8_t>, std::vector<multi_range>>> get_refinement_boundary(multi_range, int) const;
	hpx::future<gravity_return> gravity_solve(int pass, int level, const std::vector<double> coarse, double t, double m, bool b) const;
	hpx::future<statistics> get_statistics() const;
	hpx::future<std::vector<double>> restrict_all() const;
	hpx::future<std::vector<double>> get_gravity_flux() const;


	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & gid;
		arc & box;
	}
};

#endif /* ASTROTIGER_TREE_CLIENT_HPP_ */
