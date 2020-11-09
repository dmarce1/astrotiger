#pragma once

#include <string>
#include <vector>

class tree;

void levels_init();
void levels_add_entry(int level, tree*);
void levels_remove_entry(int level, tree*);

void levels_set_child_families(int level);
double levels_hydro_initialize(int, bool);
void levels_hydro_substep(int,int rk, double dt);
std::vector<std::string> levels_output_silo(int level, const std::string filename);
void levels_show();
double levels_fine_fluxes(int);





