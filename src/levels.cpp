#include <astrotiger/levels.hpp>

#include <astrotiger/options.hpp>
#include <astrotiger/tree.hpp>
#include <astrotiger/hpx.hpp>

#include <unordered_set>

HPX_PLAIN_ACTION (levels_set_child_families);
HPX_PLAIN_ACTION (levels_hydro_substep);

static std::vector<std::unordered_set<tree*>> levels;
static mutex_type mtx;
static std::vector<hpx::id_type> other_localities;

void levels_init() {
	levels.resize(opts.max_level + 1);
	if (hpx::get_locality_id() == 0) {
		auto localities = hpx::find_all_localities();
		other_localities.reserve(localities.size() - 1);
		for (int i = 1; i < localities.size(); i++) {
			other_localities.push_back(localities[i]);
		}
	}
}

void levels_add_entry(int level, tree *ptr) {
	std::lock_guard<mutex_type> lock(mtx);
	levels[level].insert(ptr);
}

void levels_remove_entry(int level, tree *ptr) {
	std::lock_guard<mutex_type> lock(mtx);
	levels[level].erase(ptr);
}

void levels_set_child_families(int level) {
	std::vector<hpx::future<void>> futs;
	if (hpx::get_locality_id() == 0 && other_localities.size()) {
		futs.push_back(hpx::lcos::broadcast < levels_set_child_families_action > (other_localities, level));
	}
	for (auto *ptr : levels[level]) {
		futs.push_back(hpx::async([ptr]() {
			ptr->set_child_family();
		}));
	}
	hpx::wait_all(futs.begin(), futs.end());
}

double levels_hydro_substep(int level, int rk, double dt) {
	std::vector<hpx::future<double>> futs;
	hpx::future<std::vector<double>> fut;
	if (hpx::get_locality_id() == 0 && other_localities.size()) {
		fut = hpx::lcos::broadcast < levels_hydro_substep_action > (other_localities, level, rk, dt);
	}
	for (auto *ptr : levels[level]) {
		futs.push_back(hpx::async([rk, ptr, dt]() {
			return ptr->hydro_substep(rk, dt);
		}));
	}
	double a = 0.0;
	for (int i = 0; i < futs.size(); i++) {
		a = std::max(a, futs[i].get());
	}
	if (fut.valid()) {
		auto others = fut.get();
		for (int i = 0; i < others.size(); i++) {
			a = std::max(a, others[i]);
		}
	}
	return a;
}

