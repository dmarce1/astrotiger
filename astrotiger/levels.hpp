#pragma once

#include <string>
#include <vector>

class tree;

void levels_init();
void levels_add_entry(int level, tree*);
void levels_remove_entry(int level, tree*);

void levels_get_hydro_boundaries(int level, double dt);
void levels_set_child_families(int level);
void levels_hydro_initialize(int, bool);
void levels_hydro_substep(int,int rk, double dt, bool, double, double);
void levels_apply_fine_fluxes(int);
void levels_output_silo(const std::string filename);
void levels_show();
void levels_apply_coarse_correction(int, double, double);
int levels_max_level();
double levels_fine_fluxes(int);





