#include <astrotiger/levels.hpp>

#include <astrotiger/options.hpp>
#include <astrotiger/tree.hpp>
#include <astrotiger/hpx.hpp>

#include <unordered_set>

HPX_PLAIN_ACTION (levels_set_child_families);
HPX_PLAIN_ACTION (levels_hydro_initialize);
HPX_PLAIN_ACTION (levels_hydro_substep);
HPX_PLAIN_ACTION (levels_output_silo);

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
	assert(levels[level].find(ptr) == levels[level].end());
	levels[level].insert(ptr);
}

void levels_remove_entry(int level, tree *ptr) {
	std::lock_guard<mutex_type> lock(mtx);
	assert(levels[level].find(ptr) != levels[level].end());
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

void levels_hydro_substep(int level, int rk, double dt) {
	std::vector<hpx::future<void>> futs;
	if (hpx::get_locality_id() == 0 && other_localities.size()) {
		futs.push_back(hpx::lcos::broadcast < levels_hydro_substep_action > (other_localities, level, rk, dt));
	}
	for (auto *ptr : levels[level]) {
		futs.push_back(hpx::async([ptr, rk, dt]() {
			ptr->hydro_substep(rk, dt);
		}));
	}
	hpx::wait_all(futs.begin(), futs.end());
}

void levels_show() {
	for (int l = 0; l <= opts.max_level; l++) {
		printf("%i on level %i\n", levels[l].size(), l);
	}
}

double levels_hydro_initialize(int level, bool refine) {
	std::vector<hpx::future<double>> futs;
	hpx::future<std::vector<double>> fut;
	if (hpx::get_locality_id() == 0 && other_localities.size()) {
		fut = hpx::lcos::broadcast < levels_hydro_initialize_action > (other_localities, level, refine);
	}
	for (auto *ptr : levels[level]) {
		futs.push_back(hpx::async([ptr, refine]() {
			return ptr->hydro_initialize(refine);
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

std::vector<std::string> levels_output_silo(int level, const std::string filename) {
	DBfile *db;
	std::vector<std::string> mnames;
	if (level == 0 && hpx::get_locality_id() == 0) {
		db = DBCreateReal(filename.c_str(), DB_CLOBBER, DB_LOCAL, "Astro-Tiger", DB_PDB);
	} else {
		db = DBOpenReal(filename.c_str(), DB_PDB, DB_APPEND);
	}
	if (db == nullptr) {
		printf("Unable to open %s for writing\n", filename.c_str());
	} else {
		for (auto *ptr : levels[level]) {
			mnames.push_back(ptr->output(db));
		}
		SILO_CHECK(DBClose(db));
		auto localities = hpx::find_all_localities();
		if (hpx::get_locality_id() != localities.size() - 1) {
			auto other_names = levels_output_silo_action()(localities[hpx::get_locality_id() + 1], level, filename);
			mnames.insert(mnames.end(), other_names.begin(), other_names.end());
		}
	}
	return mnames;
}

