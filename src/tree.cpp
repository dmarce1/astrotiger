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
using delist_action_type = tree::delist_action;
using initialize_action_type = tree::initialize_action;
using get_children_action_type = tree::get_children_action;
using get_ptr_action_type = tree::get_ptr_action;
using truncate_action_type = tree::truncate_action;
using get_box_action_type = tree::get_box_action;
using get_refinement_boundary_action_type = tree::get_refinement_boundary_action;
using list_action_type = tree::list_action;
HPX_REGISTER_ACTION (list_action_type);
HPX_REGISTER_ACTION (get_refinement_boundary_action_type);
HPX_REGISTER_ACTION (get_box_action_type);
HPX_REGISTER_ACTION (truncate_action_type);
HPX_REGISTER_ACTION (get_ptr_action_type);
HPX_REGISTER_ACTION (get_hydro_boundary_action_type);
HPX_REGISTER_ACTION (get_hydro_prolong_action_type);
HPX_REGISTER_ACTION (get_hydro_restrict_action_type);
HPX_REGISTER_ACTION (set_family_action_type);
HPX_REGISTER_ACTION (delist_action_type);
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

multi_range tree::get_box() const {
	return box;
}

std::vector<double> tree::get_hydro_boundary(multi_range b, int this_step) {
	while (hydro_step % (opts.nrk + 1) != this_step % (opts.nrk + 1)) {
		hpx::this_thread::yield();
	}
	return hydro.pack(b);
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

std::pair<std::vector<double>, std::vector<double>> tree::get_hydro_restrict() {
	std::pair<std::vector<double>, std::vector<double>> rc;
	const auto bbox = box.half();
	rc.first = hydro.pack_restrict(bbox);
	rc.second = hydro.pack_coarse_flux();
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
		new_tree.t = new_tree.t0 = t;
		new_tree.refine_step = 1;
		new_tree.hydro.unpack(hydro.pack(new_box), new_box);
		for (auto &c : children) {
			if (!c.get_box().intersection(new_box.double_()).empty()) {
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

double tree::hydro_initialize(bool refine) {
	sanity();
	std::vector<multi_range> force_refine_boxes;
	std::vector<hpx::future<std::pair<std::vector<double>, std::vector<double>>>> futs;
	for (int i = 0; i < children.size(); i++) {
		futs.push_back(children[i].get_hydro_restrict());
	}
	for (int i = 0; i < children.size(); i++) {
		const auto tmp = futs[i].get();
		const auto cbox = children[i].get_box().half();
		if (dt != 0.0) {
			hydro.unpack_coarse_flux(tmp.second, cbox, dt);
		}
		hydro.unpack(tmp.first, cbox);
	}
	hydro_step++;
	get_hydro_boundaries(true,0);
	if (refine && level < opts.max_level) {
		std::vector<std::vector<tree_client>> grandchildren;
		std::vector<hpx::future<std::vector<tree_client>>> cfuts;
		for (const auto &c : children) {
			cfuts.push_back(c.get_children());
		}
		grandchildren.resize(children.size());
		for (int i = 0; i < children.size(); i++) {
			auto &f = cfuts[i];
			const auto tmp = f.get();
			for (const auto &gc : tmp) {
				grandchildren[i].push_back(gc);
				grandchild_boxes.push_back(gc.get_box());
			}
		}
		hydro.compute_refinement_criteria();
		refine_step++;

		force_refine_boxes = grandchild_boxes;
		for (auto &b : force_refine_boxes) {
			for (int dim = 0; dim < NDIM; dim++) {
				b.min[dim] = b.min[dim] / 4;
				b.max[dim] = (b.max[dim] + 3) / 4;
			}
		}
		auto tmp1 = get_refinement_boundaries();
		force_refine_boxes.insert(force_refine_boxes.end(), tmp1.begin(), tmp1.end());
		auto new_boxes = hydro.refined_ranges(get_amr_boxes(), force_refine_boxes);
		std::vector<std::vector<tree_client>> old_grandchildren;
		std::vector<tree_client> old_children;
		std::vector<tree_client> new_children;
		for (auto &b : new_boxes) {
			b = b.double_();
		}
		std::vector<multi_range> tmp;
		std::vector<hpx::future<void>> void_futs;
		std::vector < hpx::future < hpx::id_type >> new_children_futs;
		for (int i = 0; i < children.size(); i++) {
			auto &c = children[i];
			bool found = false;
			for (int j = 0; j < new_boxes.size(); j++) {
				const auto &b = new_boxes[j];
				if (b == c.get_box()) {
					new_children.push_back(c);
					found = true;
					new_boxes[j] = new_boxes.back();
					new_boxes.pop_back();
					break;
				}
			}
			if (!found) {
				old_children.push_back(c);
				void_futs.push_back(c.delist());
				old_grandchildren.push_back(std::move(grandchildren[i]));
			}
		}
		const int unchanged_cnt = new_children.size();
		std::vector<hpx::future<std::shared_ptr<tree>>> cpfuts;
		std::vector<std::shared_ptr<tree>> old_ptrs;
		for (const auto &c : old_children) {
			cpfuts.push_back(c.get_ptr());
		}
		for (auto &f : cpfuts) {
			old_ptrs.push_back(f.get());
		}
		for (int i = 0; i < new_boxes.size(); i++) {
			const auto &b = new_boxes[i];
			auto np = std::make_shared<tree>(level + 1, b);
			np->hydro.resize(np->dx, b.pad(opts.hbw));
			np->t = np->t0 = t;
			np->refine_step = 1;
			np->hydro.unpack(hydro.pack_prolong(b, 1.0), b);
			for (const auto &op : old_ptrs) {
				const auto inter = np->box.intersection(op->box);
				if (!inter.empty()) {
					np->hydro.unpack(op->hydro.pack(inter), inter);
				}
			}
			for (int j = 0; j < old_grandchildren.size(); j++) {
				for (auto &gc : old_grandchildren[j]) {
					if (!gc.get_box().intersection(b.double_()).empty()) {
						np->children.push_back(gc);
					}
				}
			}
			new_children_futs.push_back(hpx::new_<tree>(hpx::find_here(), std::move(*np)));
		}
		for (int i = 0; i < new_boxes.size(); i++) {
			new_children.push_back(tree_client(new_children_futs[i].get(), new_boxes[i]));
		}
		hpx::wait_all(void_futs.begin(), void_futs.end());
		children = std::move(new_children);
		std::vector<hpx::future<tree_client>> tfuts;
		for (int i = unchanged_cnt; i < children.size(); i++) {
			assert(children[i] != tree_client());
			tfuts.push_back(children[i].truncate(children[i].get_box()));
		}
		int j = 0;
		void_futs.resize(0);
		for (int i = unchanged_cnt; i < children.size(); i++) {
			children[i] = tfuts[j].get();
			void_futs.push_back(children[i].list());
			j++;
		}
		hpx::wait_all(void_futs.begin(), void_futs.end());
	}

	return hydro.compute_flux(0);

}

void tree::get_hydro_boundaries(bool amr, int rk) {
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
	if (level > 0 && amr) {
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
			sib_futs.push_back(siblings[i].client.get_hydro_boundary(shifted_box, hydro_step));
		}
	}
	if (parent != tree_client() && amr) {
		assert(parent.get_box().double_().contains(box));
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
		hydro_step(0), refine_step(0) {
	level = -1;
	dt = 0.0;
}

tree::~tree() {
}

tree::tree(int level_, multi_range box_) :
		hydro_step(0), refine_step(0) {
	dt = 0.0;
	t0 = 0.0;
	t = 0.0;
	level = level_;
//	printf("Adding entry %i\n", level);
	box = box_;
	dx = 1.0 / (1 << level) / opts.max_box;
	hydro_step = 0;
}

void tree::set_family(tree_client p, tree_client s, std::vector<sibling> sibs) {
	parent = p;
	self = s;
	siblings = std::move(sibs);
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
			std::vector<multi_range> tmp;
			std::vector<multi_range> these_ranges(1, bnd_ranges[i]);
			//		printf( "%i\n", siblings.size());
			for (const auto &sib : siblings) {
				tmp.resize(0);
				for (int j = 0; j < these_ranges.size(); j++) {
					auto new_ranges = these_ranges[j].subtract(sib.box());
					for (const auto &b : new_ranges) {
						assert(b.volume());
					}
					tmp.insert(tmp.end(), new_ranges.begin(), new_ranges.end());
				}
				these_ranges = std::move(tmp);
			}
			amr_ranges.insert(amr_ranges.end(), these_ranges.begin(), these_ranges.end());
		}
	}
//	for( int i = 0; i < amr_ranges.size(); i++) {
//		printf( "%s %s\n",  box.to_string().c_str(), amr_ranges[i].to_string().c_str());
//	}
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

void tree::hydro_substep(int rk, double this_dt) {
	dt = this_dt;
	if (rk == 0) {
		refine_step = 1;
	} else {
		hydro.compute_flux(rk);
	}
	hydro.substep_update(rk, dt);
	hydro_step++;
	if (rk != opts.nrk - 1) {
		get_hydro_boundaries(true, rk + 1);
	} else {
		t0 = t;
		t += dt;
	}
}

double tree::initialize(int this_level) {
	double amax = 0.0;
	if (this_level == level) {
		levels_add_entry(level, this);
		hydro.resize(dx, box.pad(opts.hbw));
		hydro.initialize();
		amax = hydro.compute_flux(0);
		hydro.reset_flux_registers();
		hydro.reset_coarse_flux_registers();
	} else {
		if (this_level == level + 1) {
			hydro.compute_refinement_criteria();
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
	}
	return amax;
}

std::string tree::output(DBfile *db) const {
	auto options = DBMakeOptlist(1);
	int one = 1;
	DBAddOption(options, DBOPT_HIDE_FROM_GUI, &one);
	std::array<double*, NDIM> coords;
	std::array<int, NDIM> dims1;
	std::array<int, NDIM> dims2;
	char *coordnames[NDIM];
	std::vector<std::vector<double>> vars;
	vars.resize(opts.nhydro);
	for (int dim = 0; dim < NDIM; dim++) {
		dims1[dim] = box.dims()[dim]/* + 2* opts.hbw*/;
		dims2[dim] = dims1[dim] + 1;
		coordnames[dim] = new char[2];
		coordnames[dim][0] = 'x' + dim;
		coordnames[dim][1] = '\0';
		coords[dim] = new double[dims2[dim]];
		for (int i = box.min[dim]/* - opts.hbw*/; i <= box.max[dim] /*+ opts.hbw*/; i++) {
			coords[dim][i - box.min[dim]/* + opts.hbw*/] = hydro.coord(i) - 0.5 * dx;
		}
	}
	for (int f = 0; f < opts.nhydro; f++) {
		vars[f] = hydro.pack_field(f, box/*.pad(opts.hbw)*/);
	}
	std::string mesh_name;
	for (int dim = 0; dim < NDIM; dim++) {
		mesh_name += std::to_string(level) + std::string("_") + std::string(coordnames[dim]) + "_";
		mesh_name += std::to_string(box.min[dim]) + "_" + std::to_string(box.max[dim]);
		if (dim != NDIM - 1) {
			mesh_name += std::string("_");
		}
	}
	SILO_CHECK(DBPutQuadmesh(db, mesh_name.c_str(), coordnames, coords.data(), dims2.data(), NDIM, DB_DOUBLE, DB_COLLINEAR, options));

	auto var_names = hydro_grid::field_names();
	for (int f = 0; f < opts.nhydro; f++) {
		var_names[f] += "_" + mesh_name;
	}

	for (int f = 0; f < opts.nhydro; f++) {
		SILO_CHECK(DBPutQuadvar1(db, var_names[f].c_str(), mesh_name.c_str(), vars[f].data(), dims1.data(), NDIM, NULL, 0, DB_DOUBLE, DB_ZONECENT,options));
	}

	for (int dim = 0; dim < NDIM; dim++) {
		delete[] coordnames[dim];
		delete[] coords[dim];
	}
	DBFreeOptlist(options);
	return mesh_name;
}

