/*
 * tree.cpp
 *
 *  Created on: Oct 15, 2020
 *      Author: dmarce1
 */

#include <astrotiger/tree.hpp>

HPX_REGISTER_COMPONENT(hpx::components::managed_component<tree>, tree);

using get_hydro_boundary_action_type = tree::get_hydro_boundary_action;
using get_hydro_prolong_action_type = tree::get_hydro_prolong_action;
using get_hydro_restrict_action_type = tree::get_hydro_restrict_action;
using set_family_action_type = tree::set_family_action;
using clear_family_action_type = tree::clear_family_action;
using initialize_action_type = tree::initialize_action;
using get_children_action_type = tree::get_children_action;
HPX_REGISTER_ACTION (get_hydro_boundary_action_type);
HPX_REGISTER_ACTION (get_hydro_prolong_action_type);
HPX_REGISTER_ACTION (get_hydro_restrict_action_type);
HPX_REGISTER_ACTION (set_family_action_type);
HPX_REGISTER_ACTION (clear_family_action_type);
HPX_REGISTER_ACTION (initialize_action_type);
HPX_REGISTER_ACTION (get_children_action_type);

hpx::future<tree_client> tree::allocate(int level, multi_range box) {
	auto fut1 = hpx::new_<tree>(hpx::find_here(), level, box);
	auto fut2 = fut1.then([box](hpx::future<hpx::id_type> fut) {
		auto id = fut.get();
		return tree_client(id, box);
	});
	return fut2;
}

std::vector<double> tree::get_hydro_boundary(multi_range b, int this_step) {
	while (step % (opts.nrk + 1) != this_step % (opts.nrk + 1)) {
		hpx::this_thread::yield();
	}
	return hydro.pack_boundary(b);
}

std::vector<std::vector<double>> tree::get_hydro_prolong(std::vector<multi_range> ranges, double this_t) {
	double w;
	if (t0 != t) {
		w = (t - this_t) / (t - t0);
	} else {
		w = 1.0;
	}
	std::vector<std::vector<double>> p(ranges.size());
	for (int j = 0; j < ranges.size(); j++) {
		p[j] = hydro.pack_prolong(ranges[j], w);
	}
	return p;
}

std::vector<double> tree::get_hydro_restrict() {
	return hydro.pack_restrict(box.half());
}

double tree::hydro_initialize() {
	std::vector<hpx::future<std::vector<double>>> futs;
	for (int i = 0; i < children.size(); i++) {
		futs.push_back(children[i].get_hydro_restrict());
	}
	for (int i = 0; i < children.size(); i++) {
		hydro.unpack(futs[i].get(), children[i].get_box());
	}
	step++;
	get_hydro_boundaries(false);
	return hydro.compute_flux();

}

void tree::get_hydro_boundaries(bool amr) {
	std::vector<multi_range> bnd_ranges;
	std::vector<multi_range> parent_ranges;
	std::vector<multi_range> sib_ranges;
	for (int dim = 0; dim < NDIM; dim++) {
		auto this_box = box;
		this_box.max[dim] = this_box.min[dim];
		this_box.min[dim] -= opts.hbw;
		bnd_ranges.push_back(this_box);
		this_box = box;
		this_box.min[dim] = this_box.max[dim];
		this_box.max[dim] += opts.hbw;
		bnd_ranges.push_back(this_box);
	}
	for (const auto &sib : siblings) {
		multi_range inter;
		inter.set_null();
		for (int i = 0; i < bnd_ranges.size(); i++) {
			inter = sib.box().intersection(bnd_ranges[i]);
			if (!inter.empty()) {
				break;
			}
		}
		sib_ranges.push_back(inter);
	}
	if (level > 0 && amr) {
		for (int i = 0; i < bnd_ranges.size(); i++) {
			std::vector<multi_range> tmp;
			std::vector<multi_range> these_ranges(1, bnd_ranges[i]);
			//		printf( "Siblings %i\n", siblings.size());
			for (const auto &sib : siblings) {
				tmp.resize(0);
				for (int j = 0; j < these_ranges.size(); j++) {
					auto new_ranges = these_ranges[j].subtract(sib.box());
					tmp.insert(tmp.end(), new_ranges.begin(), new_ranges.end());
				}
				these_ranges = std::move(tmp);
			}
			parent_ranges.insert(parent_ranges.end(), these_ranges.begin(), these_ranges.end());
		}
	}
	std::vector<hpx::future<std::vector<double>>> sib_futs;
	for (int i = 0; i < sib_ranges.size(); i++) {
		if (!sib_ranges[i].empty()) {
			auto shifted_box = sib_ranges[i].shift(-siblings[i].shift);
			sib_futs.push_back(siblings[i].client.get_hydro_boundary(shifted_box, step));
		}
	}
	if (parent != tree_client() && amr) {
		auto parent_fut = parent.get_hydro_prolong(parent_ranges, t);
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
}

tree::tree() :
		step(0) {
	level = -1;
}

tree::~tree() {
	if (level >= 0) {
		levels_remove_entry(level, this);
	}
}

tree::tree(int level_, multi_range box_) :
		step(0) {
	t0 = 0.0;
	t = 0.0;
	level = level_;
	levels_add_entry(level, this);
	box = box_;
	dx = 1.0 / (1 << level) / opts.max_box;
	step = 0;
}

void tree::set_family(tree_client p, tree_client s, std::vector<sibling> sibs) {
	parent = p;
	self = s;
	siblings = std::move(sibs);
}

void tree::clear_family() {
	parent = tree_client();
	self = tree_client();
	siblings.resize(0);
	children.resize(0);
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

void tree::hydro_substep(int rk, double dt) {
	if (rk > 0) {
		hydro.compute_flux();
	}
	hydro.substep_update(rk, dt);
	step++;
	get_hydro_boundaries(true);
	if (rk == opts.nrk - 1) {
		t0 = t;
		t += dt;
	}
}

double tree::initialize(int this_level) {
	double amax = 0.0;
	if (this_level == level) {
		hydro.resize(dx, box.pad(opts.hbw));
		hydro.initialize();
		amax = hydro.compute_flux();
	} else {
		if (this_level == level + 1) {
			hydro.compute_refinement_criteria();
			auto boxes = hydro.refined_ranges();

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
	}
	return amax;
}

void tree::output(DBfile *db) const {
	std::array<double*, NDIM> coords;
	std::array<int, NDIM> dims1;
	std::array<int, NDIM> dims2;
	char *coordnames[NDIM];
	std::vector<std::vector<double>> vars;
	vars.resize(opts.nhydro);
	for (int dim = 0; dim < NDIM; dim++) {
		dims1[dim] = box.dims()[dim];
		dims2[dim] = dims1[dim] + 1;
		coordnames[dim] = new char[2];
		coordnames[dim][0] = 'x' + dim;
		coordnames[dim][1] = '\0';
		coords[dim] = new double[dims2[dim]];
		for (int i = box.min[dim]; i <= box.max[dim]; i++) {
			coords[dim][i - box.min[dim]] = hydro.coord(i) - 0.5 * dx;
		}
	}
	for (int f = 0; f < opts.nhydro; f++) {
		vars[f] = hydro.pack_field(f);
	}
	std::string mesh_name;
	for (int dim = 0; dim < NDIM; dim++) {
		mesh_name += std::string(coordnames[dim]) + "_";
		mesh_name += std::to_string(box.min[dim]) + "_" + std::to_string(box.max[dim]);
		if (dim != NDIM - 1) {
			mesh_name += std::string("_");
		}
	}
	SILO_CHECK(DBPutQuadmesh(db, mesh_name.c_str(), coordnames, coords.data(), dims2.data(), NDIM, DB_DOUBLE, DB_COLLINEAR, NULL));

	auto var_names = hydro_grid::field_names();
	for (int f = 0; f < opts.nhydro; f++) {
		var_names[f] += "_" + mesh_name;
	}

	for (int f = 0; f < opts.nhydro; f++) {
//		printf("%s\n", var_names[f].c_str());
		SILO_CHECK(DBPutQuadvar1(db, var_names[f].c_str(), mesh_name.c_str(), vars[f].data(), dims1.data(), NDIM, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL));
	}

	for (int dim = 0; dim < NDIM; dim++) {
		delete[] coordnames[dim];
		delete[] coords[dim];
	}
}

