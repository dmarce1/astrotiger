#include <astrotiger/geometry.hpp>

geometry::geometry(const multi_range &box_) {
	box = box_;

	multi_range dirs(multi_range(0).pad(1));
	for (multi_iterator i(dirs); !i.end(); i++) {
		bool cont = false;
		bool all_zero = true;
		const auto this_dir = i.index();
		for (int dim = 0; dim < NDIM; dim++) {
			if (this_dir[dim] != 0) {
				all_zero = false;
			}
			if (all_zero) {
				continue;
			}
			bool found = false;
			for (const auto &d : directions) {
				if (d == -this_dir) {
					found = true;
					break;
				}
			}
			if (!found) {
				directions.push_back(this_dir);
			}
		}
	}




}
