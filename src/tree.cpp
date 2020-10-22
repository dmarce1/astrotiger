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
using get_grandchild_boxes_action_type = tree::get_grandchild_boxes_action;
using get_ptr_action_type = tree::get_ptr_action;
using truncate_action_type = tree::truncate_action;
HPX_REGISTER_ACTION (truncate_action_type);
HPX_REGISTER_ACTION (get_ptr_action_type);
HPX_REGISTER_ACTION (get_grandchild_boxes_action_type);
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
	while (hydro_step % (opts.nrk + 1) != this_step % (opts.nrk + 1)) {
		hpx::this_thread::yield();
	}
	return hydro.pack(b);
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

std::vector<double> tree::get_hydro_restrict() {
	return hydro.pack_restrict(box.half());
}

std::vector<multi_range> tree::get_grandchild_boxes(int step) const {
	while (step % 2 != refine_step % 2) {
//		printf( "%i %i\n", (int) step, (int) refine_step);
		hpx::this_thread::yield();
	}
	return grandchild_boxes;
}

std::shared_ptr<tree> tree::get_ptr() {
	std::shared_ptr<tree> ptr(this, [](const tree*) {
	});
	return ptr;
}

tree_client tree::truncate(tree_client self, multi_range trunc_box) {
//	printf("%s %s\n", box.to_string().c_str(), trunc_box.to_string().c_str());
	std::vector<hpx::future<tree_client>> futs;
	const auto new_box = trunc_box.intersection(box);
	for (const auto &c : children) {
		if (!c.get_box().intersection(new_box).empty()) {
			futs.push_back(c.truncate(trunc_box));
		} else {
			futs.push_back(hpx::make_ready_future(c));
		}
	}
	for (int i = 0; i < children.size(); i++) {
		if (!children[i].get_box().intersection(new_box).empty()) {
			children[i] = futs[i].get();
		}
	}
	tree_client rclient;
	if (new_box == box) {
		rclient = self;
	} else {
		clear_family();
		tree new_tree(level + 1, new_box);
		new_tree.hydro.resize(new_tree.dx, new_box.pad(opts.hbw));
		new_tree.t = new_tree.t0 = t;
		new_tree.hydro.unpack(hydro.pack(new_box), new_box);
		for (auto &c : children) {
			if (!c.get_box().intersection(new_box).empty()) {
				new_tree.children.push_back(c);
			}
		}
		auto fut = hpx::new_<tree>(hpx::find_here(), std::move(new_tree));
		rclient = tree_client(fut.get(), new_box);
	}
	return rclient;
}

double tree::hydro_initialize(bool refine) {
	std::vector<multi_range> force_refine_boxes;
	std::vector<hpx::future<std::vector<double>>> futs;
	for (int i = 0; i < children.size(); i++) {
		futs.push_back(children[i].get_hydro_restrict());
	}
	for (int i = 0; i < children.size(); i++) {
		hydro.unpack(futs[i].get(), children[i].get_box().half());
	}
	hydro.sanity(box);
	if (refine && level < opts.max_level) {
		std::vector<std::vector<tree_client>> grandchildren;
		std::vector<hpx::future<std::vector<tree_client>>> cfuts;
		std::vector<hpx::future<std::vector<multi_range>>> sfuts;
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
//		printf("refine_step++\n");
		refine_step++;
		for (const auto &s : siblings) {
			sfuts.push_back(s.client.get_grandchild_boxes(refine_step % 2));
		}
		force_refine_boxes = grandchild_boxes;
//		printf("sibling grandchildren\n");
		for (int i = 0; i < siblings.size(); i++) {
			auto tmp = sfuts[i].get();
			for (auto &b : tmp) {
				b = b.shift(siblings[i].shift);
			}
			force_refine_boxes.insert(force_refine_boxes.begin(), tmp.begin(), tmp.end());
		}
//		printf("sibling grandchildren done\n");
		for (auto &b : force_refine_boxes) {
			for (int dim = 0; dim < NDIM; dim++) {
				b.min[dim] = b.min[dim] / 4;
				b.max[dim] = (b.max[dim] + 3) / 4;
			}
		}
		hydro.compute_refinement_criteria(force_refine_boxes);
		auto new_boxes = hydro.refined_ranges(get_amr_boxes());
		std::vector<std::vector<tree_client>> old_grandchildren;
		std::vector<tree_client> old_children;
		std::vector<tree_client> new_children;
		for (auto &b : new_boxes) {
			b = b.double_();
		}
		std::vector<multi_range> tmp;
		std::vector<hpx::future<void>> old_children_futs;
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
				old_children_futs.push_back(c.clear_family());
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
//		printf("%i %i\n", unchanged_cnt, new_boxes.size());
//		printf("Making\n");
		for (int i = 0; i < new_boxes.size(); i++) {
			const auto &b = new_boxes[i];
			auto np = std::make_shared<tree>(level + 1, b);
			np->hydro.resize(np->dx, b.pad(opts.hbw));
			np->t = np->t0 = t;
			np->hydro.unpack(hydro.pack_prolong(b, 1.0), b);
			for (const auto &op : old_ptrs) {
				const auto inter = np->box.intersection(op->box);
				if (!inter.empty()) {
					np->hydro.unpack(op->hydro.pack(inter), inter);
				}
			}
			for (int j = 0; j < old_grandchildren.size(); j++) {
				for (auto &gc : old_grandchildren[j]) {
		//			printf( "%s %s\n", gc.get_box().half().to_string().c_str(),new_boxes[i].to_string().c_str());
					if (!gc.get_box().double_().intersection(new_boxes[i]).empty()) {
						np->children.push_back(gc);
					}
				}
			}
			new_children_futs.push_back(hpx::new_<tree>(hpx::find_here(), std::move(*np)));
		}
//		printf("made\n");
		for (int i = 0; i < new_boxes.size(); i++) {
			new_children.push_back(tree_client(new_children_futs[i].get(), new_boxes[i]));
		}
//		printf("moving\n");
		children = std::move(new_children);
		std::vector<hpx::future<tree_client>> tfuts;
		for (int i = unchanged_cnt; i < children.size(); i++) {
			tfuts.push_back(children[i].truncate(children[i].get_box()));
		}
		int j = 0;
		for (int i = unchanged_cnt; i < children.size(); i++) {
			children[i] = tfuts[j].get();
			j++;
		}
//		printf("Waiting\n");
		hpx::wait_all(old_children_futs.begin(), old_children_futs.end());
		printf("Done refining level %i\n", level);
	}

	hydro_step++;
	get_hydro_boundaries(true);

	return hydro.compute_flux();

}

void tree::get_hydro_boundaries(bool amr) {
	std::vector<multi_range> bnd_ranges;
	std::vector<multi_range> parent_ranges;
	std::vector<vect<index_type>> parent_shifts;
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
		const int nx = (1 << level) * opts.max_box;
		parent_shifts.resize(parent_ranges.size());
		int j = 0;
		for (auto &r : parent_ranges) {
			parent_shifts[j] = vect<index_type>(0);
			for (int dim = 0; dim < NDIM; dim++) {
				if (r.min[dim] < 0) {
					r.min[dim] += nx;
					r.max[dim] += nx;
					parent_shifts[j][dim] -= nx;
				} else if (r.max[dim] > nx) {
					r.min[dim] -= nx;
					r.max[dim] -= nx;
					parent_shifts[j][dim] += nx;
				}
			}
			//		printf("%s %i %i\n", r.to_string().c_str(), parent_shifts[j][0], parent_shifts[j][1]);
			j++;
		}
	}
	std::vector<hpx::future<std::vector<double>>> sib_futs;
	for (int i = 0; i < sib_ranges.size(); i++) {
		if (!sib_ranges[i].empty()) {
			auto shifted_box = sib_ranges[i].shift(-siblings[i].shift);
			sib_futs.push_back(siblings[i].client.get_hydro_boundary(shifted_box, hydro_step));
		}
	}
	if (parent != tree_client() && amr) {
		auto parent_fut = parent.get_hydro_prolong(parent_ranges, t);
		const auto datas = parent_fut.get();
		for (int i = 0; i < datas.size(); i++) {
			hydro.unpack(datas[i], parent_ranges[i].shift(parent_shifts[i]));
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
		hydro_step(0), refine_step(0) {
	level = -1;
}

tree::~tree() {
	if (level >= 0) {
//		printf("Removing entry %i\n", level);
		levels_remove_entry(level, this);
	}
}

tree::tree(int level_, multi_range box_) :
		hydro_step(0), refine_step(0) {
	t0 = 0.0;
	t = 0.0;
	level = level_;
//	printf("Adding entry %i\n", level);
	levels_add_entry(level, this);
	box = box_;
	dx = 1.0 / (1 << level) / opts.max_box;
	hydro_step = 0;
}

void tree::set_family(tree_client p, tree_client s, std::vector<sibling> sibs) {
	parent = p;
	self = s;
	siblings = std::move(sibs);
}

void tree::clear_family() {
//	printf("Clearing family\n");
	std::vector<hpx::future<void>> futs;
	for (const auto &c : children) {
		futs.push_back(c.clear_family());
	}
	{
		std::lock_guard<mutex_type> lock(mtx);
		parent = tree_client();
		self = tree_client();
		siblings.resize(0);
	}
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
	if (rk == 0) {
		if (refine_step % 2 == 1) {
			refine_step++;
		}
	}
	if (rk > 0) {
		hydro.compute_flux();
	}
	hydro.substep_update(rk, dt);
	hydro_step++;
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
			hydro.compute_refinement_criteria(std::vector<multi_range>());
			auto boxes = hydro.refined_ranges(get_amr_boxes());

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
//		printf("%s\n", var_names[f].c_str());
		SILO_CHECK(DBPutQuadvar1(db, var_names[f].c_str(), mesh_name.c_str(), vars[f].data(), dims1.data(), NDIM, NULL, 0, DB_DOUBLE, DB_ZONECENT,options));
	}

	for (int dim = 0; dim < NDIM; dim++) {
		delete[] coordnames[dim];
		delete[] coords[dim];
	}
	DBFreeOptlist(options);
	return mesh_name;
}

