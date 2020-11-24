/*
 * tree.cpp
 *
 *  Created on: Oct 15, 2020
 *      Author: dmarce1
 */

#include <astrotiger/tree.hpp>

#include <astrotiger/cosmos.hpp>

#include <unordered_map>

HPX_REGISTER_COMPONENT(hpx::components::managed_component<tree>, tree);

using get_hydro_prolong_action_type = tree::get_hydro_prolong_action;
using get_hydro_restrict_action_type = tree::get_hydro_restrict_action;
using set_family_action_type = tree::set_family_action;
using delist_action_type = tree::delist_action;
using initialize_action_type = tree::initialize_action;
using get_children_action_type = tree::get_children_action;
using get_ptr_action_type = tree::get_ptr_action;
using truncate_action_type = tree::truncate_action;
using get_box_action_type = tree::get_box_action;
using get_refinement_boundary_action_type = tree::get_refinement_boundary_action;
using list_action_type = tree::list_action;
using get_energy_boundary_action_type = tree::get_energy_boundary_action;
using get_energy_prolong_action_type = tree::get_energy_prolong_action;
using gravity_solve_action_type = tree::gravity_solve_action;
using get_statistics_action_type = tree::get_statistics_action;
using restrict_all_action_type = tree::restrict_all_action;
using set_boundary_action_type = tree::set_boundary_action;
using compute_error_action_type = tree::compute_error_action;
using get_fine_flux_action_type = tree::get_fine_flux_action;
using drift_action_type = tree::drift_action;
using recv_parts_action_type = tree::recv_parts_action;
using finish_drift_action_type = tree::finish_drift_action;
using get_particle_count_action_type = tree::get_particle_count_action;
using max_part_velocity_action_type = tree::max_part_velocity_action;
using kick_action_type = tree::kick_action;
using compute_cic_action_type = tree::compute_cic_action;
using get_energy_statistics_action_type = tree::get_energy_statistics_action;
using get_average_phi_action_type = tree::get_average_phi_action;
HPX_REGISTER_ACTION (get_average_phi_action_type);
using set_average_phi_action_type = tree::set_average_phi_action;
HPX_REGISTER_ACTION (set_average_phi_action_type);
HPX_REGISTER_ACTION (get_energy_statistics_action_type);
HPX_REGISTER_ACTION (compute_cic_action_type);
HPX_REGISTER_ACTION (kick_action_type);
HPX_REGISTER_ACTION (max_part_velocity_action_type);
HPX_REGISTER_ACTION (get_particle_count_action_type);
HPX_REGISTER_ACTION (finish_drift_action_type);
HPX_REGISTER_ACTION (recv_parts_action_type);
HPX_REGISTER_ACTION (drift_action_type);
HPX_REGISTER_ACTION (get_fine_flux_action_type);
HPX_REGISTER_ACTION (compute_error_action_type);
HPX_REGISTER_ACTION (set_boundary_action_type);
HPX_REGISTER_ACTION (restrict_all_action_type);
HPX_REGISTER_ACTION (get_statistics_action_type);
HPX_REGISTER_ACTION (gravity_solve_action_type);
HPX_REGISTER_ACTION (get_energy_prolong_action_type);
HPX_REGISTER_ACTION (get_energy_boundary_action_type);
HPX_REGISTER_ACTION (list_action_type);
HPX_REGISTER_ACTION (get_refinement_boundary_action_type);
HPX_REGISTER_ACTION (get_box_action_type);
HPX_REGISTER_ACTION (truncate_action_type);
HPX_REGISTER_ACTION (get_ptr_action_type);
HPX_REGISTER_ACTION (get_hydro_prolong_action_type);
HPX_REGISTER_ACTION (get_hydro_restrict_action_type);
HPX_REGISTER_ACTION (set_family_action_type);
HPX_REGISTER_ACTION (delist_action_type);
HPX_REGISTER_ACTION (initialize_action_type);
HPX_REGISTER_ACTION (get_children_action_type);

double tree::get_average_phi(int lev) const {
	if (level == lev) {
		return grav.get_phi_tot();
	} else {
		double tot = 0.0;
		std::vector<hpx::future<double>> futs;
		for (const auto &c : children) {
			futs.push_back(c.get_average_phi(lev));
		}
		for (auto &f : futs) {
			tot += f.get();
		}
		return tot;
	}
}

void tree::set_average_phi(int lev, double dif) {
	if (level == lev) {
		grav.set_average_phi(dif);
		hydro.set_phi(grav.get_phi());
	} else {
		std::vector<hpx::future<void>> futs;
		for (const auto &c : children) {
			futs.push_back(c.set_average_phi(lev, dif));
		}
		hpx::wait_all(futs.begin(), futs.end());
	}
}

energy_statistics tree::get_energy_statistics(double rho0) const {
	energy_statistics e, tmp;
	e.ekin = e.epot = 0.0;
	std::vector<multi_range> boxes;
	std::vector<hpx::future<energy_statistics>> futs;
	for (const auto c : children) {
		boxes.push_back(c.get_box().half());
		futs.push_back(c.get_energy_statistics(rho0));
	}
	for (auto &f : futs) {
		tmp = f.get();
		e.ekin += tmp.ekin;
		e.epot += tmp.epot;
	}
	if (opts.particles) {
		tmp = parts.get_energy_statistics(grav.get_phi());
		e.ekin += tmp.ekin;
		e.epot += tmp.epot;
	}
	if (opts.hydro) {
		tmp = hydro.get_energy_statistics(grav.get_phi(), boxes, rho0);
		e.ekin += tmp.ekin;
		e.epot += tmp.epot;
	}
	return e;
}

std::vector<double> tree::compute_cic(const std::vector<double> &coarse, double this_t, int this_level) {
	gravity_step = 0;
	parts.compute_cloud_in_cell(this_t - t);
	if (level > 0) {
		parts.unpack_cic_prolong(coarse, box.half());
	}
	std::vector<hpx::future<std::vector<double>>> cfuts;
	for (const auto &c : children) {
		const auto cbox = c.get_box().half();
		const auto tmp = parts.pack_cic(cbox);
		parts.zero_cic(cbox);
		cfuts.push_back(c.compute_cic(tmp, this_t, this_level));
	}
	for (int ci = 0; ci < children.size(); ci++) {
		const auto tmp = cfuts[ci].get();
		parts.unpack_cic(tmp, children[ci].get_box().half().pad(1));
	}
	if (level == this_level) {
		gravity_step++;
		for (int i = 0; i < siblings.size(); i++) {
			const auto inter = box.pad(1).intersection(siblings[i].box());
			if (inter.volume()) {
//			printf( "packing %s\n", inter.to_string().c_str());
				siblings[i].client.set_boundary(parts.pack_cic(inter), box.shift(-siblings[i].shift), gravity_step);
			}
		}
		for (int i = 0; i < siblings.size(); i++) {
			const auto inter = box.intersection(siblings[i].box().pad(1));
			if (inter.volume()) {
//			printf( "unpacking %s\n", inter.to_string().c_str());
				parts.unpack_cic(siblings[i].client.get(gravity_step).get(), inter);
			}
		}
	}

	return parts.pack_cic_restrict();
}

void tree::kick(int rung, double tm, std::vector<double> last_dt, std::vector<double> this_dt) {
	std::vector<hpx::future<void>> futs;
	for (const auto &c : children) {
		futs.push_back(c.kick(rung, tm, last_dt, this_dt));
	}
	double w;
	if (t != t0) {
		w = (tm - t0) / (t - t0);
	} else {
		w = 1.0;
	}
	assert(w >= 0.0);
	assert(w <= 1.0);
	parts.kick(rung, level, last_dt, this_dt, grav.get_acceleration(w));
	hpx::wait_all(futs.begin(), futs.end());
}

double tree::max_part_velocity() const {
	std::vector<hpx::future<double>> futs;
	for (const auto &c : children) {
		futs.push_back(c.max_part_velocity());
	}
	double vmax = parts.max_velocity();
	for (auto &f : futs) {
		vmax = std::max(vmax, f.get());
	}
	return vmax;
}

void tree::recv_parts(std::vector<particle> these_parts) {
	std::vector<hpx::future<void>> futs;
	for (int i = 0; i < siblings.size(); i++) {
		std::vector<particle> sibparts;
		const auto rbox = range_int_to_double(siblings[i].box());
		int j = 0;
		while (j < these_parts.size()) {
			const auto &part = these_parts[j];
			if (rbox.contains(part.x)) {
				sibparts.push_back(part);
				const int sz = these_parts.size() - 1;
				these_parts[j] = these_parts[sz];
				these_parts.resize(sz);
			} else {
				j++;
			}
		}
		if (sibparts.size()) {
			if (siblings[i].client != self) {
				futs.push_back(siblings[i].client.send_parts(std::move(sibparts)));
			} else {
				these_parts.insert(these_parts.end(), sibparts.begin(), sibparts.end());
			}
		}
	}
	for (auto &p : these_parts) {
		for (int dim = 0; dim < NDIM; dim++) {
			while (p.x[dim] >= 1.0) {
				p.x[dim] -= 1.0;
			}
			while (p.x[dim] < 0.0) {
				p.x[dim] += 1.0;
			}
		}
	}
#ifndef NDEBUG
	const auto rbox = range_int_to_double(box);
	for (const auto &p : these_parts) {
		assert(rbox.contains(p.x));
	}
#endif
	{
		std::lock_guard<mutex_type> lock(mtx);
		new_parts.insert(new_parts.end(), these_parts.begin(), these_parts.end());
	}
	hpx::wait_all(futs.begin(), futs.end());
}

multi_array<int> tree::get_particle_count() const {
	std::vector<hpx::future<multi_array<int>>> futs;
	for (const auto &c : children) {
		futs.push_back(c.get_particle_count());
	}
	auto count = parts.particle_count();
	for (int i = 0; i < children.size(); i++) {
		const auto tmp = futs[i].get();
		const auto cbox = children[i].get_box();
		for (multi_iterator j(cbox); !j.end(); j++) {
			count[j.index() / 2] += tmp[j];
		}
	}
	return count;
}

double tree::compute_error() {

	double error = 0.0;

	std::vector<hpx::future<double>> futs;
	std::vector<multi_range> cboxes;
	for (const auto &c : children) {
		futs.push_back(c.compute_error());
		cboxes.push_back(c.get_box().half());
	}
	multi_array<double> results;
	if (opts.problem == "sphere") {
		error = hydro.compare_analytic(cboxes, results);
		hydro.set_error_field(std::move(results));
	}

	for (auto &f : futs) {
		error += f.get();
	}
	return error;

}

range<double> tree::range_int_to_double(const multi_range &box) {
	range<double> rbox;
	for (int dim = 0; dim < NDIM; dim++) {
		rbox.min[dim] = box.min[dim] * dx;
		rbox.max[dim] = box.max[dim] * dx;
	}
	return rbox;
}

void tree::finish_drift(std::vector<particle> parent_parts) {
	std::vector<hpx::future<void>> futs;
	new_parts.insert(new_parts.end(), parent_parts.begin(), parent_parts.end());
	for (int ci = 0; ci < children.size(); ci++) {
		std::vector<particle> cparts;
		const auto &c = children[ci];
		const auto rbox = range_int_to_double(c.get_box().half());
		int i = 0;
		while (i < new_parts.size()) {
			const auto &part = new_parts[i];
			if (rbox.contains(part.x)) {
				cparts.push_back(part);
				const int sz = new_parts.size() - 1;
				new_parts[i] = new_parts[sz];
				new_parts.resize(sz);
			} else {
				i++;
			}
		}
		futs.push_back(c.finish_drift(std::move(cparts)));
	}
	parts.add_parts(new_parts);
	new_parts = std::vector<particle>();

	hpx::wait_all(futs.begin(), futs.end());
}

void tree::drift(double dt) {
	std::vector<hpx::future<void>> futs;
	for (const auto &c : children) {
		futs.push_back(c.drift(dt));
	}
	std::vector<multi_range> child_boxes;
	for (const auto &c : children) {
		child_boxes.push_back(c.get_box().half());
	}
	parts.set_child_boxes(child_boxes);
	auto escaped = parts.drift(dt);
	for (int i = 0; i < siblings.size(); i++) {
		std::vector<particle> sibparts;
		const auto rbox = range_int_to_double(siblings[i].box());
		int j = 0;
		while (j < escaped.size()) {
			const auto &part = escaped[j];
			if (rbox.contains(part.x)) {
				sibparts.push_back(part);
				const int sz = escaped.size() - 1;
				escaped[j] = escaped[sz];
				escaped.resize(sz);
			} else {
				j++;
			}
		}
		if (sibparts.size()) {
			futs.push_back(siblings[i].client.send_parts(std::move(sibparts)));
		}
	}
	for (int i = 0; i < children.size(); i++) {
		std::vector<particle> cparts;
		const auto rbox = range_int_to_double(children[i].get_box().half());
		int j = 0;
		while (j < escaped.size()) {
			const auto &part = escaped[j];
			if (rbox.contains(part.x)) {
				cparts.push_back(part);
				const int sz = escaped.size() - 1;
				escaped[j] = escaped[sz];
				escaped.resize(sz);
			} else {
				j++;
			}
		}
		if (cparts.size()) {
			futs.push_back(children[i].send_parts(std::move(cparts)));
		}
	}
	if (escaped.size()) {
		assert(parent != tree_client());
		futs.push_back(parent.send_parts(std::move(escaped)));
	}
	hpx::wait_all(futs.begin(), futs.end());
}

statistics tree::get_statistics(int lev, double t) {
	const auto a = cosmos_a();
	statistics stats;
	std::vector<multi_range> cranges;
	std::vector<hpx::future<statistics>> futs;
	int volume = 0;
	for (const auto &c : children) {
		futs.push_back(c.get_statistics(lev, t));
		const auto this_box = c.get_box().half();
		volume += this_box.volume();
		if (level < lev) {
			cranges.push_back(box);
		}
	}
	if (level <= lev) {
		stats = hydro.get_statistics(cranges);
	}
	if (children.size() == 0) {
		stats.min_level = level;
		stats.max_level = level;
	} else {
		stats.min_level = opts.max_level;
		stats.max_level = 0;
	}
	for (auto &f : futs) {
		auto tmp = f.get();
		if (level < lev) {
			for (int field = 0; field < opts.nhydro; field++) {
				stats.u[field] += tmp.u[field];
			}
		}
		if (volume == box.volume()) {
			stats.min_level = std::min(stats.min_level, tmp.min_level);
		} else {
			stats.min_level = level;
		}
		stats.max_level = std::max(stats.max_level, tmp.max_level);
	}
	double pm = 0.0;
	if (opts.particles && level == 0) {
		compute_cic(std::vector<double>(), t, 0);
		double m = 0.0;
		const auto rho = parts.get_cic();
		for (multi_iterator i(box); !i.end(); i++) {
			m += std::pow(dx, NDIM) * rho[i];
		}
		pm += m;
//		printf( "%e %e\n", stats.u[rho_i], pm);
		stats.u[rho_i] += pm;
	}
	return stats;
}

std::atomic<int> counter(0);

gravity_return tree::gravity_solve(int pass, int fine_level, const std::vector<double> coarse_from_parent, double this_t, double mtot, int iters) {
	gravity_step = 0;
	if (level == fine_level) {
		if (pass == 0) {
			if (level == 0) {
				grav.initialize_coarse(t);
			}
			auto source = hydro.get_density();
			if (opts.particles) {
				const auto psource = parts.get_cic();
				for (multi_iterator i(box); !i.end(); i++) {
					source[i] += psource[i];
				}
			}
			grav.initialize_fine(source, mtot, level);
			if (level != 0) {
				grav.set_amr_zones(get_amr_boxes(), coarse_from_parent);
			} else {
				if (opts.problem == "sphere") {
					grav.set_outflow_boundaries();
				}
			}
		} else {
			if (level != 0) {
				grav.apply_prolong(coarse_from_parent);
			}
		}
	} else if (pass == 0) {
		grav.initialize_coarse(t);
	} else {
		if (level != 0) {
			grav.apply_prolong(coarse_from_parent);
		}
	}
	if (pass > 0 || level == fine_level) {
		for (int i = 0; i < iters; i++) {
			get_gravity_boundaries(PACK_POTENTIAL_REDBLACK);
			grav.relax();
		}
	}
	gravity_return rc;
	rc.resid = 0.0;
	double vmax = 0.0;
	if (level < fine_level) {
		std::vector<hpx::future<gravity_return>> futs;
		futs.reserve(children.size());
		double w;
		if (t0 != t) {
			w = (this_t - t0) / (t - t0);
			assert(this_t >= t0);
			assert(this_t <= t);
		} else {
			w = 1.0;
		}
		for (int i = 0; i < children.size(); i++) {
			const auto &c = children[i];
			if (pass == 0) {
				const auto this_box = c.get_box().half().pad(opts.gbw);
				futs.push_back(c.gravity_solve(pass, fine_level, grav.pack(this_box, w), this_t, mtot, iters));
			} else {
				futs.push_back(c.gravity_solve(pass, fine_level, grav.get_prolong(c.get_box().half()), this_t, mtot, iters));
			}
		}
		for (auto &f : futs) {
			auto tmp = f.get();
			grav.apply_restrict(tmp);
			rc.resid = std::max(tmp.resid, rc.resid);
			vmax = std::max(vmax, tmp.vmax);
		}
		for (int i = 0; i < iters; i++) {
			int type = PACK_POTENTIAL_REDBLACK;
			if (i == 0) {
				type |= (PACK_ACTIVE | PACK_SOURCE);
			}
			get_gravity_boundaries(type);
			grav.relax();
		}
	} else if (level == fine_level) {
		get_gravity_boundaries(PACK_POTENTIAL);
		if (pass == GRAVITY_FINAL_PASS) {
			vmax = grav.finish_fine();
			hydro.set_phi(grav.get_phi());
		}
	}
	if (level == 0 && fine_level == 0 && opts.problem != "sphere") {
		grav.set_avg_zero();
	}
	auto tmp = grav.get_restrict(mtot);
	if (level != fine_level) {
		tmp.resid = rc.resid;
	}
	tmp.vmax = vmax;
	return tmp;

}

hpx::future<tree_client> tree::allocate(int level, multi_range box) {
	auto fut1 = hpx::new_<tree>(hpx::find_here(), level, box);
	auto fut2 = fut1.then([box](hpx::future<hpx::id_type> fut) {
		auto id = fut.get();
		return tree_client(id, box);
	});
	return fut2;
}

multi_range tree::get_box() const {
	return box;
}

void tree::set_boundary(std::vector<double> &&data, const multi_range &bbox, int step) {
	bool found = false;
	for (int i = 0; i < siblings.size(); i++) {
		if (siblings[i].box() == bbox) {
			assert(!found);
			siblings[i].client.put(std::move(data), step);
			found = true;
			break;
		}
	}
	assert(found);
}

std::vector<double> tree::get_energy_boundary(multi_range b, int this_step) {
	while (energy_step % 3 != this_step % 3) {
//		printf( "%i %i %i\n", level, (int) energy_step, this_step);
		hpx::this_thread::yield();
	}
	return hydro.pack_field(tau_i, b);
}

std::pair<std::vector<std::uint8_t>, std::vector<multi_range>> tree::get_refinement_boundary(multi_range b, int this_step) {
	std::pair<std::vector<std::uint8_t>, std::vector<multi_range>> rc;
	while (refine_step % 3 != this_step % 3) {
		//	printf("%i %i\n", (int) refine_step, this_step);
		hpx::this_thread::yield();
	}
	rc.first = hydro.pack_refinement(b);
	rc.second = grandchild_boxes;
	return rc;
}

std::vector<std::vector<double>> tree::get_hydro_prolong(std::vector<multi_range> ranges, double this_t) {
	double w;
//	printf("%.14e %.14e %.14e\n", t0, this_t, t);
	if (t0 != t) {
		w = (this_t - t0) / (t - t0);
		assert(this_t >= t0);
		assert(this_t <= t);
	} else {
		w = 1.0;
	}
	std::vector<std::vector<double>> p(ranges.size());
	for (int j = 0; j < ranges.size(); j++) {
		p[j] = hydro.pack_prolong(ranges[j], w);
	}
	return p;
}

std::vector<std::vector<double>> tree::get_energy_prolong(std::vector<multi_range> ranges, double this_t) {
	double w;
	if (t0 != t) {
		w = (this_t - t0) / (t - t0);
		assert(this_t >= t0);
		assert(this_t <= t);
	} else {
		w = 1.0;
	}
	std::vector<std::vector<double>> p(ranges.size());
	for (int j = 0; j < ranges.size(); j++) {
		p[j] = hydro.pack_field_prolong(tau_i, ranges[j], w);
	}
	return p;
}

std::pair<std::vector<double>, std::vector<double>> tree::get_hydro_restrict() {
	std::pair<std::vector<double>, std::vector<double>> rc;
	const auto bbox = box.half();
//	printf( "1 %s\n", box.to_string().c_str());
//	printf( "2 %s\n", bbox.to_string().c_str());
	rc.first = hydro.pack_restrict(bbox);
	rc.second = hydro.pack_coarse_correction();
	return rc;
}

std::shared_ptr<tree> tree::get_ptr() {
	std::shared_ptr<tree> ptr(this, [](const tree*) {
	});
	return ptr;
}

tree_client tree::truncate(tree_client self, multi_range trunc_box) {
//	printf("TRUNCATE %i %s %s %i\n", level, box.to_string().c_str(), trunc_box.to_string().c_str(), children.size());
	const auto new_box = trunc_box.intersection(box);
	assert(new_box.volume());

	tree_client rclient;
	if (new_box == box) {
		std::vector<hpx::future<tree_client>> futs;
		for (const auto &c : children) {
			assert(c != tree_client());
			if (!c.get_box().intersection(new_box.double_()).empty()) {
				futs.push_back(c.truncate(new_box.double_()));
			}
		}
		std::vector<tree_client> new_children;
		for (int i = 0; i < futs.size(); i++) {
			new_children.push_back(futs[i].get());
		}
		children = std::move(new_children);
		rclient = self;
	} else {
		tree new_tree(level, new_box);
		new_tree.hydro.resize(new_tree.dx, new_box.pad(opts.hbw));
		new_tree.parts.resize(new_tree.dx, new_box);
		new_tree.grav.resize(new_tree.dx, new_box);
		new_tree.t = t;
		new_tree.t0 = t0;
		new_tree.refine_step = 1;
		new_tree.hydro.unpack(hydro.pack(new_box), new_box);
		if (opts.particles) {
			const auto these_parts = parts.get_particles();
			range<double> rbox;
			for (int dim = 0; dim < NDIM; dim++) {
				rbox.min[dim] = new_box.min[dim] * dx;
				rbox.max[dim] = new_box.max[dim] * dx;
			}
			std::vector<particle> new_parts;
			for (const auto &p : these_parts) {
				if (rbox.contains(p.x)) {
					new_parts.push_back(p);
				}
			}
			new_tree.parts.add_parts(new_parts);
		}
		for (auto &c : children) {
			if (c.get_box().intersection(new_box.double_()).volume()) {
				new_tree.children.push_back(c);
			}
		}
		auto fut = hpx::new_<tree>(hpx::find_here(), std::move(new_tree));
		rclient = tree_client(fut.get(), new_box);
		rclient.truncate(new_box).get();
	}
	return rclient;
}

void tree::sanity() const {
#ifndef NDEBUG
	const int max_x = (1 << level) * opts.max_box;
	for (int dim = 0; dim < NDIM; dim++) {
		assert(box.min[dim] >= 0 && box.max[dim] <= max_x);
	}
#endif
}

std::vector<double> tree::restrict_all() {
	std::vector<hpx::future<std::vector<double>>> futs;
	for (int i = 0; i < children.size(); i++) {
		futs.push_back(children[i].restrict_all());
	}
	for (int i = 0; i < children.size(); i++) {
		const auto tmp = futs[i].get();
		const auto cbox = children[i].get_box().half();
		hydro.unpack(tmp, cbox);
	}
	get_hydro_boundaries(t);
	hydro_step = 0;
	return hydro.pack_restrict(box.half());
}

std::vector<double> tree::get_fine_flux() {
	get_hydro_boundaries(t);
	return hydro.pack_coarse_flux();
}

void tree::apply_coarse_correction(double a0, double a1) {
	std::vector<hpx::future<std::pair<std::vector<double>, std::vector<double>>>> futs;
	std::vector<hpx::future<std::vector<tree_client>>> sfuts;
	for (const auto &s : siblings) {
		sfuts.push_back(s.client.get_children());
	}
	std::vector<multi_index> shift(children.size());
	for (int i = 0; i < children.size(); i++) {
		shift[i] = multi_index(0);
	}
	std::vector<tree_client> these_children = children;
	for (int si = 0; si < siblings.size(); si++) {
		const auto tmp = sfuts[si].get();
		for (int i = 0; i < tmp.size(); i++) {
			shift.push_back(siblings[si].shift);
		}
		these_children.insert(these_children.end(), tmp.begin(), tmp.end());
	}
	for (const auto &c : these_children) {
		futs.push_back(c.get_hydro_restrict());
	}

	std::vector<std::pair<std::vector<double>, std::vector<double>>> u;
	std::vector<multi_range> used_boxes;
	for (int i = 0; i < these_children.size(); i++) {
		const auto cbox = these_children[i].get_box().half().shift(shift[i]);
		u.push_back(futs[i].get());
		bool used = false;
	//	printf( "%s\n", cbox.to_string().c_str());
		for (int j = 0; j < used_boxes.size(); j++) {
			if (used_boxes[j] == cbox) {
				used = true;
				break;
			}
		}
		if (!used) {
			used_boxes.push_back(cbox);
			if (dt != 0.0) {
				hydro.unpack_coarse_correction(u[i].second, cbox, dt, a0, a1);
			}
		} else {
			printf( "Used\n");
		}
	}

//	hydro.update_energy();
//	energy_step++;
//	get_energy_boundaries(t);
	for (int i = 0; i < children.size(); i++) {
		const auto cbox = children[i].get_box().half();
		hydro.unpack(u[i].first, cbox);
	}
	get_hydro_boundaries(t);
	hydro_step = 0;
}

void tree::energy_update() {
	hydro.update_energy();
	energy_step++;
//	printf( "%i\n", (int) energy_step);
	get_energy_boundaries(t);
}

void tree::hydro_initialize(bool refine, bool start) {
	hydro_step = 0;
	energy_step = 0;
	std::vector<multi_range> force_refine_boxes;
	hpx::future<multi_array<int>> pfut;
	if (start) {
		hydro.reset_coarse_flux_registers();
	}
	hydro.store();
//	if (level < opts.max_level && refine) {
//		if (opts.particles) {
//			pfut = hpx::async([this]() {
//				return this->get_particle_count();
//			});
//		} else {
//			pfut = hpx::make_ready_future(multi_array<int>());
//		}
//		std::vector<std::vector<tree_client>> grandchildren;
//		std::vector<hpx::future<std::vector<tree_client>>> cfuts;
//		for (const auto &c : children) {
//			cfuts.push_back(c.get_children());
//		}
//		grandchildren.resize(children.size());
//		grandchild_boxes.resize(0);
//		for (int i = 0; i < children.size(); i++) {
//			auto &f = cfuts[i];
//			const auto tmp = f.get();
//			for (const auto &gc : tmp) {
//				grandchildren[i].push_back(gc);
//				grandchild_boxes.push_back(gc.get_box());
//			}
//		}
//		hydro.compute_refinement_criteria(pfut.get());
//		refine_step++;
//
//		force_refine_boxes = grandchild_boxes;
//		for (auto &b : force_refine_boxes) {
//			for (int dim = 0; dim < NDIM; dim++) {
//				b.min[dim] = b.min[dim] / 4;
//				b.max[dim] = (b.max[dim] + 3) / 4;
//			}
//		}
//		auto tmp1 = get_refinement_boundaries();
//		force_refine_boxes.insert(force_refine_boxes.end(), tmp1.begin(), tmp1.end());
//		auto new_boxes = hydro.refined_ranges(get_amr_boxes(), force_refine_boxes);
//		std::vector<std::vector<tree_client>> old_grandchildren;
//		std::vector<tree_client> old_children;
//		std::vector<tree_client> new_children;
//		for (auto &b : new_boxes) {
//			b = b.double_();
//		}
//		std::vector<multi_range> tmp;
//		std::vector<hpx::future<void>> void_futs;
//		std::vector < hpx::future < hpx::id_type >> new_children_futs;
//		for (int i = 0; i < children.size(); i++) {
//			auto &c = children[i];
//			bool found = false;
//			for (int j = 0; j < new_boxes.size(); j++) {
//				const auto &b = new_boxes[j];
//				if (b == c.get_box()) {
//					new_children.push_back(c);
//					found = true;
//					new_boxes[j] = new_boxes.back();
//					new_boxes.pop_back();
//					break;
//				}
//			}
//			if (!found) {
//				old_children.push_back(c);
//				void_futs.push_back(c.delist());
//				old_grandchildren.push_back(std::move(grandchildren[i]));
//			}
//		}
//		const int unchanged_cnt = new_children.size();
//		std::vector<hpx::future<std::shared_ptr<tree>>> cpfuts;
//		std::vector<std::shared_ptr<tree>> old_ptrs;
//		for (const auto &c : old_children) {
//			cpfuts.push_back(c.get_ptr());
//		}
//		for (auto &f : cpfuts) {
//			old_ptrs.push_back(f.get());
//		}
//		if (opts.particles) {
//			for (auto &op : old_ptrs) {
//				auto these_parts = op->parts.get_particles();
//				parts.add_parts(these_parts);
//			}
//		}
//		for (int i = 0; i < new_boxes.size(); i++) {
//			const auto &b = new_boxes[i];
//			auto np = std::make_shared<tree>(level + 1, b);
//			np->hydro.resize(np->dx, b.pad(opts.hbw));
//			np->parts.resize(np->dx, b);
//			np->grav.resize(np->dx, b);
//			np->refine_step = 1;
//			np->hydro.unpack(hydro.pack_prolong(b.pad(opts.hbw), 1.0), b.pad(opts.hbw));
//			np->t = np->t0 = t;
//			for (const auto &op : old_ptrs) {
//				np->t = op->t;
//				np->t0 = op->t0;
//				const auto inter = np->box.intersection(op->box);
//				if (!inter.empty()) {
//					np->hydro.unpack(op->hydro.pack(inter), inter);
//				}
//			}
//			for (int j = 0; j < old_grandchildren.size(); j++) {
//				for (auto &gc : old_grandchildren[j]) {
//					if (!gc.get_box().intersection(b.double_()).empty()) {
//						np->children.push_back(gc);
//					}
//				}
//			}
//			if (opts.particles) {
//				auto these_parts = parts.get_particles(b.half());
//				np->parts.add_parts(these_parts);
//			}
//			new_children_futs.push_back(hpx::new_<tree>(hpx::find_here(), std::move(*np)));
//		}
//		for (int i = 0; i < new_boxes.size(); i++) {
//			new_children.push_back(tree_client(new_children_futs[i].get(), new_boxes[i]));
//		}
//		hpx::wait_all(void_futs.begin(), void_futs.end());
//		children = std::move(new_children);
//		std::vector<hpx::future<tree_client>> tfuts;
//		for (int i = unchanged_cnt; i < children.size(); i++) {
//			assert(children[i] != tree_client());
//			tfuts.push_back(children[i].truncate(children[i].get_box()));
//		}
//		int j = 0;
//		void_futs.resize(0);
//		for (int i = unchanged_cnt; i < children.size(); i++) {
//			children[i] = tfuts[j].get();
//			void_futs.push_back(children[i].list());
//			j++;
//		}
//		hpx::wait_all(void_futs.begin(), void_futs.end());
//	}
}

double tree::apply_fine_fluxes() {
	get_hydro_boundaries(t);
	double amax = 0.0;
	std::vector<hpx::future<std::vector<double>>> cfuts;
	for (const auto &c : children) {
//		cfuts.push_back(c.get_fine_flux());
	}
	const auto alocal = hydro.compute_flux(0);
	double achild = 0.0;
	const auto amaxp1 = std::max(amax, hydro.positivity_limit());
//	for (int i = 0; i < children.size(); i++) {
//		achild = std::max(achild, hydro.unpack_fine_flux(cfuts[i].get(), children[i].get_box().half()));
//	}
//	const auto amaxp2 = std::max(amax, hydro.positivity_limit());
	amax = std::max(amaxp1, std::max(achild, alocal));
//	amax = std::max(amaxp2, amax);
//	printf( "%e %e %e\n", amaxp1, amaxp2, alocal);
	return amax;
}

void tree::get_hydro_boundaries(double this_time) {
//	printf( "%i\n", (int) hydro_step);
	for (const auto &sib : siblings) {
		const auto inter = sib.box().pad(opts.hbw).intersection(box);
		if (inter.volume()) {
			sib.client.set_boundary(hydro.pack(inter), box.shift(-sib.shift), hydro_step);
		}
	}

	std::vector<multi_range> bnd_ranges;
	std::vector<multi_range> parent_ranges;
	std::vector<multi_range> sib_ranges;
	bnd_ranges = box.pad(opts.hbw).subtract(box);
	for (const auto &sib : siblings) {
		multi_range inter;
		inter.set_null();
		for (int i = 0; i < bnd_ranges.size(); i++) {
			inter = sib.box().intersection(bnd_ranges[i]);
			sib_ranges.push_back(inter);
		}
	}
	if (level > 0) {
		for (int i = 0; i < bnd_ranges.size(); i++) {
			std::vector<multi_range> tmp;
			std::vector<multi_range> these_ranges(1, bnd_ranges[i]);
			for (const auto &sib : siblings) {
				tmp.resize(0);
				for (int j = 0; j < these_ranges.size(); j++) {
					auto new_ranges = these_ranges[j].subtract(sib.box());
					tmp.insert(tmp.end(), new_ranges.begin(), new_ranges.end());
				}
				these_ranges = std::move(tmp);
			}
			for (int j = 0; j < these_ranges.size(); j++) {
				assert(parent.get_box().pad(opts.max_bw).double_().contains(these_ranges[j]));
			}
			parent_ranges.insert(parent_ranges.end(), these_ranges.begin(), these_ranges.end());
		}
	}
	auto tmp = std::move(sib_ranges);
	sib_ranges.resize(0);
	for (const auto &s : siblings) {
		multi_range this_range;
		this_range.set_null();
		for (const auto &range : tmp) {
			const auto inter = range.intersection(s.box());
			this_range = this_range.union_(inter);
		}
		sib_ranges.push_back(this_range);
	}
	std::vector<hpx::future<std::vector<double>>> sib_futs;
	for (int i = 0; i < sib_ranges.size(); i++) {
		if (!sib_ranges[i].empty()) {
			auto shifted_box = sib_ranges[i].shift(-siblings[i].shift);
			sib_futs.push_back(siblings[i].client.get(hydro_step));
		}
	}
	if (parent != tree_client()) {
		assert(parent.get_box().double_().contains(box));
		auto parent_fut = parent.get_hydro_prolong(parent_ranges, this_time);
		const auto datas = parent_fut.get();
		for (int i = 0; i < datas.size(); i++) {
			hydro.unpack(datas[i], parent_ranges[i]);
		}
	}
	int j = 0;
	for (int i = 0; i < sib_ranges.size(); i++) {
		if (!sib_ranges[i].empty()) {
			const auto data = sib_futs[j].get();
			hydro.unpack(data, sib_ranges[i]);
			j++;
		}
	}
	for (int dir = 0; dir < 2 * NDIM; dir++) {
		if (dir % 2 == 0) {
			if (box.min[dir / 2] == 0) {
				hydro.enforce_physical_bc(dir);
			}
		} else {
			if (box.max[dir / 2] == opts.max_box * (1 << level)) {
				hydro.enforce_physical_bc(dir);
			}
		}
	}
	hydro_step++;
}

void tree::get_gravity_boundaries(int type) {
	if (opts.problem == "sphere" && level == 0) {
		return;
	}
	if (level == 1) {
//		printf( "%i\n", siblings.size());
	}
//	assert( !((type & PACK_POTENTIAL) && (type != PACK_POTENTIAL)));
	for (const auto &sib : siblings) {
		const auto inter = sib.box().pad(opts.gbw).intersection(box);
		if (inter.volume()) {
			sib.client.set_boundary(grav.pack(inter, type), box.shift(-sib.shift), gravity_step);
		}
	}
	for (auto &sib : siblings) {
		const auto inter = sib.box().intersection(box.pad(opts.gbw));
		if (inter.volume()) {
			grav.unpack(sib.client.get(gravity_step).get(), inter, type);
		}
	}
	gravity_step++;
}

void tree::get_energy_boundaries(double this_time) {
	std::vector<multi_range> bnd_ranges;
	std::vector<multi_range> parent_ranges;
	std::vector<multi_range> sib_ranges;
	bnd_ranges = box.pad(1).subtract(box);
	for (const auto &sib : siblings) {
		multi_range inter;
		inter.set_null();
		for (int i = 0; i < bnd_ranges.size(); i++) {
			inter = sib.box().intersection(bnd_ranges[i]);
			sib_ranges.push_back(inter);
		}
	}
	if (level > 0) {
		for (int i = 0; i < bnd_ranges.size(); i++) {
			std::vector<multi_range> tmp;
			std::vector<multi_range> these_ranges(1, bnd_ranges[i]);
			for (const auto &sib : siblings) {
				tmp.resize(0);
				for (int j = 0; j < these_ranges.size(); j++) {
					auto new_ranges = these_ranges[j].subtract(sib.box());
					tmp.insert(tmp.end(), new_ranges.begin(), new_ranges.end());
				}
				these_ranges = std::move(tmp);
			}
			for (int j = 0; j < these_ranges.size(); j++) {
				assert(parent.get_box().pad(opts.max_bw).double_().contains(these_ranges[j]));
			}
			parent_ranges.insert(parent_ranges.end(), these_ranges.begin(), these_ranges.end());
		}
	}
	auto tmp = std::move(sib_ranges);
	sib_ranges.resize(0);
	for (const auto &s : siblings) {
		multi_range this_range;
		this_range.set_null();
		for (const auto &range : tmp) {
			const auto inter = range.intersection(s.box());
			this_range = this_range.union_(inter);
		}
		sib_ranges.push_back(this_range);
	}
	std::vector<hpx::future<std::vector<double>>> sib_futs;
	for (int i = 0; i < sib_ranges.size(); i++) {
		if (!sib_ranges[i].empty()) {
			auto shifted_box = sib_ranges[i].shift(-siblings[i].shift);
			sib_futs.push_back(siblings[i].client.get_energy_boundary(shifted_box, energy_step));
		}
	}
	if (parent != tree_client()) {
		assert(parent.get_box().double_().contains(box));
		auto parent_fut = parent.get_energy_prolong(parent_ranges, this_time);
		const auto datas = parent_fut.get();
		for (int i = 0; i < datas.size(); i++) {
			hydro.unpack_field(tau_i, datas[i], parent_ranges[i]);
		}
	}
	int j = 0;
	for (int i = 0; i < sib_ranges.size(); i++) {
		if (!sib_ranges[i].empty()) {
			const auto data = sib_futs[j].get();
			hydro.unpack_field(tau_i, data, sib_ranges[i]);
			j++;
		}
	}
}

std::vector<multi_range> tree::get_refinement_boundaries() {
	std::vector<multi_range> bnd_ranges;
	std::vector<multi_range> sib_ranges;
	bnd_ranges = box.pad(opts.window).subtract(box);
	for (const auto &sib : siblings) {
		multi_range inter;
		inter.set_null();
		for (int i = 0; i < bnd_ranges.size(); i++) {
			inter = sib.box().intersection(bnd_ranges[i]);
			sib_ranges.push_back(inter);
		}
	}
	auto tmp = std::move(sib_ranges);
	sib_ranges.resize(0);
	for (const auto &s : siblings) {
		multi_range this_range;
		this_range.set_null();
		for (const auto &range : tmp) {
			const auto inter = range.intersection(s.box());
			this_range = this_range.union_(inter);
		}
		sib_ranges.push_back(this_range);
	}
	std::vector<hpx::future<std::pair<std::vector<std::uint8_t>, std::vector<multi_range>>>> sib_futs;
	for (int i = 0; i < sib_ranges.size(); i++) {
		if (!sib_ranges[i].empty()) {
			auto shifted_box = sib_ranges[i].shift(-siblings[i].shift);
			sib_futs.push_back(siblings[i].client.get_refinement_boundary(shifted_box, refine_step));
		}
	}
	int j = 0;
	std::vector<multi_range> grandchild_ranges;
	for (int i = 0; i < sib_ranges.size(); i++) {
		if (!sib_ranges[i].empty()) {
			auto data = sib_futs[j].get();
			hydro.unpack_refinement(data.first, sib_ranges[i]);
			for (auto &b : data.second) {
				for (int dim = 0; dim < NDIM; dim++) {
					b.min[dim] = b.min[dim] / 4;
					b.max[dim] = (b.max[dim] + 3) / 4;
				}
				b = b.shift(siblings[i].shift);
			}
			grandchild_ranges.insert(grandchild_ranges.end(), data.second.begin(), data.second.end());
			j++;
		}
	}
	return grandchild_ranges;
}

tree::tree() :
		refine_step(0), energy_step(0), gravity_step(0), hydro_step(0) {
	level = -1;
	dt = 0.0;
}

tree::~tree() {
}

tree::tree(int level_, multi_range box_) :
		refine_step(0), energy_step(0), gravity_step(0), hydro_step(0) {
	dt = 0.0;
	t0 = 0.0;
	t = 0.0;
	level = level_;
//	printf("Adding entry %i\n", level);
	box = box_;
	dx = 1.0 / (1 << level) / opts.max_box;
}

void tree::set_family(tree_client p, tree_client s, std::vector<sibling> sibs) {
	parent = p;
	self = s;
	siblings = std::move(sibs);
	hydro_step = 0;
	if (p != tree_client()) {
		assert(p.get_box().double_().contains(box));
	}
}

void tree::delist() {
	levels_remove_entry(level, this);
//	printf("Clearing family\n");
	std::vector<hpx::future<void>> futs;
	for (const auto &c : children) {
		futs.push_back(c.delist());
	}
	parent = tree_client();
	self = tree_client();
	siblings.resize(0);
	hpx::wait_all(futs.begin(), futs.end());
}

void tree::list() {
	std::vector<hpx::future<void>> futs;
	for (const auto &c : children) {
		futs.push_back(c.list());
	}
	levels_add_entry(level, this);
	hpx::wait_all(futs.begin(), futs.end());
}

std::vector<multi_range> tree::get_amr_boxes() const {
	std::vector<multi_range> amr_ranges;
	if (level > 0) {
		const auto bnd_ranges = box.pad(opts.max_bw).subtract(box);
		for (int i = 0; i < bnd_ranges.size(); i++) {
			static thread_local std::vector<multi_range> tmp;
			static thread_local std::vector<multi_range> these_ranges;
			these_ranges.resize(1);
			these_ranges[0] = bnd_ranges[i];
			for (const auto &sib : siblings) {
				tmp.resize(0);
				for (int j = 0; j < these_ranges.size(); j++) {
					auto new_ranges = these_ranges[j].subtract(sib.box());
					tmp.insert(tmp.end(), new_ranges.begin(), new_ranges.end());
				}
				std::swap(these_ranges, tmp);
			}
			amr_ranges.insert(amr_ranges.end(), these_ranges.begin(), these_ranges.end());
		}
	}
	return amr_ranges;
}

void tree::set_child_family() {
	std::vector<sibling> child_sibs;
	std::vector<hpx::future<std::vector<tree_client>>> child_sib_futs;
	for (const auto &s : siblings) {
		child_sib_futs.push_back(s.client.get_children());
	}
	for (int i = 0; i < siblings.size(); i++) {
		auto csibs = child_sib_futs[i].get();
		for (const auto &cs : csibs) {
			sibling sib;
			sib.client = cs;
			sib.shift = siblings[i].shift * 2;
			child_sibs.push_back(std::move(sib));
		}
	}
	for (int i = 0; i < children.size(); i++) {
		sibling sib;
		sib.client = children[i];
		sib.shift = 0;
		child_sibs.push_back(sib);
	}
	std::vector<hpx::future<void>> set_futs;
	for (const auto &c : children) {
		assert(box.double_().contains(c.get_box()));
		std::vector<sibling> these_sibs;
		auto cbox = c.get_box();
		cbox = cbox.pad(opts.max_bw);
		for (const auto &cs : child_sibs) {
			const auto inter = cbox.intersection(cs.box());
			if (!inter.empty()) {
				if (cs.client != c || cs.shift != vect<index_type>(0)) {
					these_sibs.push_back(cs);
				}
			}
		}
		set_futs.push_back(c.set_family(self, c, std::move(these_sibs)));
	}
	hpx::wait_all(set_futs.begin(), set_futs.end());
}

std::vector<tree_client> tree::get_children() const {
	return children;
}

void tree::hydro_substep(int rk, double this_dt, bool last, double a1, double a2) {
	dt = this_dt;
	if (rk == 0) {
		hydro.store();
		refine_step = 1;
		energy_step = 0;
	}
	if (rk == opts.nrk - 1) {
		hydro.compute_flux(1);
		hydro.store_flux();
	}
	hydro.substep_update(rk, dt, a1, a2);
	if (rk == opts.nrk - 1) {
//		if (!last) {
//			hydro.update_energy();
//			energy_step++;
//			get_energy_boundaries(t);
//		}
		t0 = t;
		t += dt;
		last_dt = dt;
	}
	const double next_t = ((rk != opts.nrk - 1) ? (t + opts.alpha[rk + 1] * this_dt) : t);
	get_hydro_boundaries(next_t);
	energy_step = 0;
}

double tree::initialize(int this_level) {
	double amax = 0.0;
	hpx::future<multi_array<int>> pfut;
	if (this_level == level) {
		levels_add_entry(level, this);
		hydro.resize(dx, box.pad(opts.hbw));
		parts.resize(dx, box);
		grav.resize(dx, box);
		hydro.initialize();
		if (opts.particles && level == 0) {
			parts.initialize();
		}
		hydro.reset_flux_registers();
		hydro.reset_coarse_flux_registers();
		amax = hydro.compute_flux(0);
	} else {
		if (this_level == level + 1) {
			if (opts.particles) {
				pfut = hpx::async([this]() {
					return this->get_particle_count();
				});
			} else {
				pfut = hpx::make_ready_future(multi_array<int>());
			}
			hydro.compute_refinement_criteria(pfut.get());
			refine_step++;
			get_refinement_boundaries();
			auto boxes = hydro.refined_ranges(get_amr_boxes(), std::vector<multi_range>());

			std::vector<hpx::future<tree_client>> futs(boxes.size());
			children.resize(boxes.size());
			for (int i = 0; i < boxes.size(); i++) {
				futs[i] = allocate(level + 1, boxes[i].double_());
			}
			children.resize(boxes.size());
			for (int i = 0; i < futs.size(); i++) {
				children[i] = futs[i].get();
			}
		}
		std::vector<hpx::future<double>> futs(children.size());
		for (int i = 0; i < children.size(); i++) {
			futs[i] = children[i].initialize(this_level);
		}
		for (int i = 0; i < children.size(); i++) {
			amax = std::max(amax, futs[i].get());
		}
		if (opts.particles) {
			std::vector<multi_range> child_boxes;
			for (const auto &c : children) {
				child_boxes.push_back(c.get_box().half());
			}
			parts.set_child_boxes(child_boxes);
			finish_drift(parts.get_child_parts());
		}
	}
	if (opts.particles) {
		amax = std::max(amax, max_part_velocity());
	}
	return amax;
}

struct node_hash {
	std::size_t operator()(const vect<double> &node) const {
		std::hash<double> f;
		std::size_t hash = 0;
		for (int dim = 0; dim < NDIM; dim++) {
			hash ^= f(std::exp(node[dim]));
		}
		return hash;
	}
};

output_return tree::output(DBfile *db, int node_index) const {
	int offset = node_index;
	output_return rc;
	std::unordered_map<vect<double>, int, node_hash> nodes;
	std::vector<int> zones;
	multi_array<std::uint8_t> mask(box);
	for (multi_iterator i(box); !i.end(); i++) {
		mask[i] = 1;
	}
	for (const auto &c : children) {
		const auto cbox = c.get_box().half();
		for (multi_iterator i(cbox); !i.end(); i++) {
			mask[i] = 0;
		}
	}
#if NDIM == 1
	constexpr int order[2] = {0,1};
#elif NDIM == 2
	constexpr int order[4] = { 0, 1, 3, 2 };
#else
	constexpr int order[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };
#endif

	for (multi_iterator i(box); !i.end(); i++) {
		if (mask[i]) {
			std::vector<int> these_zones(1 << NDIM);
			vect<double> node;
			for (int ci = 0; ci < (1 << NDIM); ci++) {
				for (int dim = 0; dim < NDIM; dim++) {
					if ((ci >> dim) & 1) {
						node[dim] = (i.index()[dim] + 1) * dx;
					} else {
						node[dim] = (i.index()[dim]) * dx;
					}
				}
				int this_index;
				if (nodes.find(node) == nodes.end()) {
					this_index = node_index;
					nodes[node] = this_index;
					node_index++;
				} else {
					this_index = nodes[node];
				}
				these_zones[order[ci]] = this_index;
			}
			zones.insert(zones.end(), these_zones.begin(), these_zones.end());
		}
	}
	rc.zones = std::move(zones);
	rc.coords.resize(NDIM);
	for (int dim = 0; dim < NDIM; dim++) {
		rc.coords[dim].resize(nodes.size());
	}
	for (const auto &entry : nodes) {
		for (int dim = 0; dim < NDIM; dim++) {
			rc.coords[dim][entry.second - offset] = entry.first[dim];
		}
	}
	rc.data = hydro.pack_output(mask);
	if (opts.particles) {
		rc.pdata = parts.pack_output();
		rc.pcoords = parts.pack_coords();
	}
	return rc;

}

