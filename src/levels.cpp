#include <astrotiger/levels.hpp>

#include <astrotiger/options.hpp>
#include <astrotiger/tree.hpp>
#include <astrotiger/hpx.hpp>

#include <unordered_set>

HPX_PLAIN_ACTION (levels_get_hydro_boundaries);
HPX_PLAIN_ACTION (levels_set_child_families);
HPX_PLAIN_ACTION (levels_hydro_initialize);
HPX_PLAIN_ACTION (levels_hydro_substep);
HPX_PLAIN_ACTION (levels_output_silo);
HPX_PLAIN_ACTION (levels_fine_fluxes);
HPX_PLAIN_ACTION (levels_max_level);
HPX_PLAIN_ACTION (levels_energy_update);
HPX_PLAIN_ACTION (levels_apply_coarse_correction);

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

void levels_get_hydro_boundaries(int level, double t) {
	auto these_levels = levels;
	std::vector<hpx::future<void>> futs;
	if (hpx::get_locality_id() == 0 && other_localities.size()) {
		futs.push_back(hpx::lcos::broadcast < levels_get_hydro_boundaries_action > (other_localities, level, t));
	}
	for (auto *ptr : these_levels[level]) {
		futs.push_back(hpx::async([ptr, t]() {
			ptr->get_hydro_boundaries(t);
		}));
	}
	hpx::wait_all(futs.begin(), futs.end());
}

int levels_max_level() {
	hpx::future<std::vector<int>> fut;
	if (hpx::get_locality_id() == 0 && other_localities.size()) {
		fut = hpx::lcos::broadcast < levels_max_level_action > (other_localities);
	}
	int l;
	for (l = 0; l <= opts.max_level; l++) {
		if (levels[l].size() == 0) {
			break;
		}
	}
	l--;
	if (fut.valid()) {
		const auto tmp = fut.get();
		for (int i = 0; i < tmp.size(); i++) {
			l = std::max(l, tmp[i]);
		}
	}
	return l;
}

void levels_set_child_families(int level) {
	auto these_levels = levels;
	std::vector<hpx::future<void>> futs;
	if (hpx::get_locality_id() == 0 && other_localities.size()) {
		futs.push_back(hpx::lcos::broadcast < levels_set_child_families_action > (other_localities, level));
	}
	for (auto *ptr : these_levels[level]) {
		futs.push_back(hpx::async([ptr]() {
			ptr->set_child_family();
		}));
	}
	hpx::wait_all(futs.begin(), futs.end());
}
void levels_energy_update(int level) {
	auto these_levels = levels;
	std::vector<hpx::future<void>> futs;
	if (hpx::get_locality_id() == 0 && other_localities.size()) {
		futs.push_back(hpx::lcos::broadcast < levels_energy_update_action > (other_localities, level));
	}
	for (auto *ptr : these_levels[level]) {
		futs.push_back(hpx::async([ptr]() {
			ptr->energy_update();
		}));
	}
	hpx::wait_all(futs.begin(), futs.end());
}

void levels_apply_coarse_correction(int level, double a0, double a1) {
	auto these_levels = levels;
	std::vector<hpx::future<void>> futs;
	if (hpx::get_locality_id() == 0 && other_localities.size()) {
		futs.push_back(hpx::lcos::broadcast < levels_apply_coarse_correction_action > (other_localities, level, a0, a1));
	}
	for (auto *ptr : these_levels[level]) {
		futs.push_back(hpx::async([a0, a1, ptr]() {
			ptr->apply_coarse_correction(a0, a1);
		}));
	}
	hpx::wait_all(futs.begin(), futs.end());
}

void levels_hydro_substep(int level, int rk, double dt, bool last, double a0, double a1) {
	auto these_levels = levels;
	std::vector<hpx::future<void>> futs;
	if (hpx::get_locality_id() == 0 && other_localities.size()) {
		futs.push_back(hpx::lcos::broadcast < levels_hydro_substep_action > (other_localities, level, rk, dt, last, a0, a1));
	}
	for (auto *ptr : these_levels[level]) {
		futs.push_back(hpx::async([ptr, rk, dt, last, a0, a1]() {
			ptr->hydro_substep(rk, dt, last, a0, a1);
		}));
	}
	hpx::wait_all(futs.begin(), futs.end());
}

void levels_show() {
	for (int l = 0; l <= opts.max_level; l++) {
		printf("%i on level %i\n", levels[l].size(), l);
	}
}

void levels_hydro_initialize(int level, bool refine) {
	auto these_levels = levels;
	std::vector<hpx::future<void>> futs;
	hpx::future<void> fut;
	if (hpx::get_locality_id() == 0 && other_localities.size()) {
		fut = hpx::lcos::broadcast < levels_hydro_initialize_action > (other_localities, level, refine);
	}
	for (auto *ptr : these_levels[level]) {
		futs.push_back(hpx::async([ptr, refine]() {
			return ptr->hydro_initialize(refine);
		}));
	}
	for (int i = 0; i < futs.size(); i++) {
		futs[i].get();
	}
	if (fut.valid()) {
		fut.get();
	}
}

double levels_fine_fluxes(int level) {
	auto these_levels = levels;
	std::vector<hpx::future<double>> futs;
	hpx::future<std::vector<double>> fut;
	if (hpx::get_locality_id() == 0 && other_localities.size()) {
		fut = hpx::lcos::broadcast < levels_fine_fluxes_action > (other_localities, level);
	}
	for (auto *ptr : these_levels[level]) {
		futs.push_back(hpx::async([ptr]() {
			return ptr->apply_fine_fluxes();
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

void levels_output_silo(const std::string filename) {
	DBfile *db = DBCreateReal(filename.c_str(), DB_CLOBBER, DB_LOCAL, "Astro-Tiger", DB_PDB);
	if (db == nullptr) {
		printf("Unable to open %s for writing\n", filename.c_str());
	} else {
		auto these_levels = levels;
		output_return output;
		output.coords.resize(NDIM);
		output.data.resize(opts.nhydro + 2);
		if (opts.particles) {
			output.pcoords.resize(NDIM);
			output.pdata.resize(NDIM + 2);
		}
		int node_index = 0;
		for (int level = 0; level <= opts.max_level; level++) {
			for (auto *ptr : these_levels[level]) {
				const output_return this_output = ptr->output(db, node_index);
				for (int dim = 0; dim < NDIM; dim++) {
					output.coords[dim].insert(output.coords[dim].end(), this_output.coords[dim].begin(), this_output.coords[dim].end());
				}
				node_index += this_output.coords[0].size();
				for (int f = 0; f < opts.nhydro + 2; f++) {
					output.data[f].insert(output.data[f].end(), this_output.data[f].begin(), this_output.data[f].end());
				}
				output.zones.insert(output.zones.end(), this_output.zones.begin(), this_output.zones.end());
				if (opts.particles) {
					for (int dim = 0; dim < NDIM; dim++) {
						output.pcoords[dim].insert(output.pcoords[dim].end(), this_output.pcoords[dim].begin(), this_output.pcoords[dim].end());
					}
					for (int f = 0; f < NDIM + 2; f++) {
						output.pdata[f].insert(output.pdata[f].end(), this_output.pdata[f].begin(), this_output.pdata[f].end());
					}
				}
			}
		}
#if NDIM==1
		const int shape = DB_ZONETYPE_BEAM;
#elif NDIM==2
		const int shape = DB_ZONETYPE_QUAD;
#else
		const int shape = DB_ZONETYPE_HEX;
#endif
		const std::vector<int> shapes(1, shape);
		const std::vector<int> shapesizes(1, 1 << NDIM);
		const std::vector<int> shapecnts(1, output.zones.size() / (1 << NDIM));
		const char *coordnames[] = { "x", "y", "z" };
		const double *coords[NDIM];
		for (int dim = 0; dim < NDIM; dim++) {
			coords[dim] = output.coords[dim].data();
		}
		const auto nnodes = output.coords[0].size();
		const auto nzones = output.zones.size() / (1 << NDIM);
		SILO_CHECK(
				DBPutZonelist2(db, "zones", nzones, NDIM, output.zones.data(), nzones * (1<<NDIM), 0, 0, 0, shapes.data(), shapesizes.data(), shapecnts.data(), 1, NULL));
		SILO_CHECK(DBPutUcdmesh (db, "amr_mesh", NDIM, coordnames, coords, nnodes, nzones, "zones", NULL, DB_DOUBLE, NULL));
		const auto names = hydro_grid::field_names();
		for (int f = 0; f < opts.nhydro + 2; f++) {
			SILO_CHECK(DBPutUcdvar1 (db, names[f].c_str(),"amr_mesh", output.data[f].data(),output.data[f].size(),NULL,0,DB_DOUBLE, DB_ZONECENT, NULL));
		}
		if (opts.particles) {
			for (int dim = 0; dim < NDIM; dim++) {
				coords[dim] = output.pcoords[dim].data();
			}
			SILO_CHECK(DBPutPointmesh(db, "point_mesh", NDIM, coords, output.pcoords[0].size(), DB_DOUBLE, NULL));
			const auto names = particles::field_names();
			for (int f = 0; f < NDIM + 2; f++) {
				SILO_CHECK(DBPutPointvar1(db,names[f].c_str(), "point_mesh", output.pdata[f].data(), output.pdata[f].size(), DB_DOUBLE, NULL));
			}
		}
		SILO_CHECK(DBClose(db));
	}
}

