#include <astrotiger/boxes.hpp>
#include <astrotiger/options.hpp>
#include <astrotiger/hpx.hpp>

#include <unordered_map>

#define TABLE_SIZE 64

struct dim_hash {
	std::size_t operator()(const vect<index_type> &dims) const {
		std::size_t hash = 0;
		for (int dim = 0; dim < NDIM; dim++) {
			hash = hash * opts.max_box + dims[dim];
		}
		return hash / TABLE_SIZE;
	}
};

const auto table_index(const vect<index_type> &dims) {
	std::size_t hash = 0;
	for (int dim = 0; dim < NDIM; dim++) {
		hash = hash * opts.max_box + dims[dim];
	}
	return hash % TABLE_SIZE;
}

struct entry_t {
	std::shared_ptr<std::vector<int>> red;
	std::shared_ptr<std::vector<int>> black;

	entry_t() = default;
	entry_t(const vect<index_type> &dims) {
		red = std::make_shared<std::vector<int>>();
		black = std::make_shared<std::vector<int>>();
		int sz = 1;
		for (int dim = 0; dim < NDIM; dim++) {
			sz *= dims[dim];
		}
		multi_range box;
		box.min = 0;
		box.max = dims;
		red->reserve(sz / 2);
		black->reserve(sz / 2);
		int index = 0;
		const auto ibox = box.pad(-1);
		for (multi_iterator i(box); !i.end(); i++) {
			if (ibox.contains(i)) {
				if (i.index().sum() % 2 == 0) {
					red->push_back(index);
				} else {
					black->push_back(index);
				}
			}
			index++;
		}
	}
};

static std::unordered_map<vect<index_type>, entry_t, dim_hash> table[TABLE_SIZE];
static mutex_type mtx[TABLE_SIZE];

static const entry_t& get_entry(const multi_index &dims, int index) {
	auto &this_table = table[index];
	auto iter = this_table.find(dims);
	if (iter == this_table.end()) {
		auto tmp = this_table.emplace(std::make_pair(dims, entry_t(dims)));
		assert(tmp.second);
		iter = tmp.first;
	}
	return iter->second;
}

std::shared_ptr<std::vector<int>> get_red_indices(const multi_range &box) {
	const auto index = table_index(box.dims());
	std::lock_guard<mutex_type> lock(mtx[index]);
	return get_entry(box.dims(), index).red;
}

std::shared_ptr<std::vector<int>> get_black_indices(const multi_range &box) {
	const auto index = table_index(box.dims());
	std::lock_guard<mutex_type> lock(mtx[index]);
	return get_entry(box.dims(), index).black;
}

