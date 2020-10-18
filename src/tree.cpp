/*
 * tree.cpp
 *
 *  Created on: Oct 15, 2020
 *      Author: dmarce1
 */

#include <astrotiger/tree.hpp>

HPX_REGISTER_COMPONENT(hpx::components::managed_component<tree>, tree);

using set_family_action_type = tree::set_family_action;
HPX_REGISTER_ACTION (set_family_action_type);

using clear_family_action_type = tree::clear_family_action;
HPX_REGISTER_ACTION (clear_family_action_type);

using initialize_action_type = tree::initialize_action;
HPX_REGISTER_ACTION (initialize_action_type);

using get_children_action_type = tree::get_children_action;
HPX_REGISTER_ACTION (get_children_action_type);

hpx::future<tree_client> tree::allocate(int level, multi_range box) {
	auto fut1 = hpx::new_<tree>(hpx::find_here(), level, box);
	auto fut2 = fut1.then([box](hpx::future<hpx::id_type> fut) {
		auto id = fut.get();
		return tree_client(id, box);
	});
	return fut2;
}

tree::tree() {
	hydro = std::make_shared<hydro_grid>();
	grids.resize(opts.ngrid);
	grids[hydro_i] = hydro;
	level = -1;
}

tree::~tree() {
	if (level >= 0) {
		levels_remove_entry(level, this);
	}
}

tree::tree(int level_, multi_range box_) {
	level = level_;
	levels_add_entry(level, this);
	box = box_;
	hydro = std::make_shared<hydro_grid>();
	grids.resize(opts.ngrid);
	grids[hydro_i] = hydro;
	dx = 1.0 / (1 << level) / opts.max_box;
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
				if (cs.client != c) {
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

void tree::initialize() {
	for (int i = 0; i < opts.ngrid; i++) {
		grids[i]->resize(dx, box.pad(opts.bw[i]));
	}
	hydro->initialize();
//	printf("Initializing grid at level %i\n", level);
//	for (int dim = 0; dim < NDIM; dim++) {
//		printf("%i %i\n", box.min[dim], box.max[dim]);
//	}
	if (opts.max_level > level) {
		hydro->compute_refinement_criteria();
		auto boxes = hydro->refined_ranges();

		std::vector<hpx::future<void>> futs(boxes.size());
		children.resize(boxes.size());
//		printf("number of boxes = %i\n", boxes.size());
		for (int i = 0; i < boxes.size(); i++) {
			futs[i] = allocate(level + 1, boxes[i].double_()).then([this, i](hpx::future<tree_client> fut) {
				children[i] = fut.get();
				children[i].initialize().get();
			});
		}
		hpx::wait_all(futs.begin(), futs.end());
	}

}
