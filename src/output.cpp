#include <astrotiger/output.hpp>
#include <astrotiger/defs.hpp>
#include <astrotiger/levels.hpp>
#include <astrotiger/hydro_grid.hpp>
#include <astrotiger/options.hpp>

#include <cstring>

#include <silo.h>

void output_silo(const std::string &filename) {
	DBfile *db;
	auto options = DBMakeOptlist(1);
	char *main_mesh = "main_mesh";
	DBAddOption(options, DBOPT_MMESH_NAME, main_mesh);
	std::vector<std::string> tmp_mnames;
	std::vector<char*> mesh_names;
	for (int l = 0; l <= opts.max_level; l++) {
		auto tmp = levels_output_silo(l, filename);
		tmp_mnames.insert(tmp_mnames.end(), tmp.begin(), tmp.end());
	}
	mesh_names.resize(tmp_mnames.size());
	std::vector<int> mesh_types(tmp_mnames.size(), DB_QUADMESH);
	for (int i = 0; i < tmp_mnames.size(); i++) {
		mesh_names[i] = new char[tmp_mnames[i].size() + 1];
		strcpy(mesh_names[i], tmp_mnames[i].c_str());
	}
	db = DBOpenReal(filename.c_str(), DB_PDB, DB_APPEND);
	SILO_CHECK(DBPutMultimesh(db, "main_mesh", mesh_names.size(), mesh_names.data(), mesh_types.data(), options));
	const auto field_names = hydro_grid::field_names();
	std::vector<int> var_types(mesh_names.size(), DB_QUADVAR);
	for (int f = 0; f < opts.nhydro; f++) {
		std::vector<char*> var_names(mesh_names.size());
		for (int i = 0; i < mesh_names.size(); i++) {
			std::string varname = field_names[f] + std::string("_") + mesh_names[i];
			var_names[i] = new char[varname.size() + 1];
			strcpy(var_names[i], varname.c_str());
		}
		DBPutMultivar(db, field_names[f].c_str(), var_names.size(), var_names.data(), var_types.data(), options);
		for (int i = 0; i < var_names.size(); i++) {
			delete[] var_names[i];
		}
	}
	for (int i = 0; i < tmp_mnames.size(); i++) {
		delete[] mesh_names[i];
	}
	SILO_CHECK(DBClose(db));
	DBFreeOptlist(options);
}
