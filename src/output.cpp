#include <astrotiger/output.hpp>
#include <astrotiger/defs.hpp>
#include <astrotiger/levels.hpp>
#include <astrotiger/hydro_grid.hpp>
#include <astrotiger/options.hpp>

#include <cstring>

#include <silo.h>

void output_silo(const std::string &filename) {
	printf( "Output\n");
	levels_output_silo(filename);
	printf( "Done output\n");
}
