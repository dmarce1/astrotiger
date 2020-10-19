#pragma once

class tree;

void levels_init();
void levels_add_entry(int level, tree*);
void levels_remove_entry(int level, tree*);

void levels_set_child_families(int level);
double levels_hydro_substep(int,int rk, double dt);






