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
using initialize_action_type = tree::initialize_action;
HPX_REGISTER_ACTION (initialize_action_type);

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
}

tree::tree(int level_, multi_range box_) {
	level = level_;
	box = box_;
	hydro = std::make_shared<hydro_grid>();
	grids.resize(opts.ngrid);
	grids[hydro_i] = hydro;
	dx = 1.0 / (1 << level) / opts.max_box;
}

void tree::set_family(tree_client p, std::vector<tree_client> c) {
	parent = p;
	children = std::move(c);
}

void tree::initialize() {
	for (int i = 0; i < opts.ngrid; i++) {
		grids[i]->resize(dx, box.pad(opts.bw[i]));
	}
	hydro->initialize();
	printf("Initializing grid at level %i\n", level);
	for (int dim = 0; dim < NDIM; dim++) {
		printf("%i %i\n", box.min[dim], box.max[dim]);
	}
	if (opts.max_level > level) {
		hydro->compute_refinement_criteria();
		auto boxes = hydro->refined_ranges();

		std::vector<hpx::future<void>> futs(boxes.size());
		children.resize(boxes.size());
		printf("number of boxes = %i\n", boxes.size());
		for (int i = 0; i < boxes.size(); i++) {
			futs[i] = allocate(level + 1, boxes[i].double_()).then([this, i](hpx::future<tree_client> fut) {
				children[i] = fut.get();
				children[i].initialize().get();
			});
		}
		hpx::wait_all(futs.begin(), futs.end());
	}

}
