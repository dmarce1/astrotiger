/*
 * id_type_impl.hpp
 *
 *  Created on: Dec 15, 2015
 *      Author: dmarce1
 */

#ifndef ID_TYPE_IMPL_HPP_
#define ID_TYPE_IMPL_HPP_

#include "hpx_fwd.hpp"

namespace hpx {

template<class Archive>
void id_type::load(const Archive &arc, const unsigned) {
	std::uintptr_t ptr;
	arc >> rank;
	arc >> ptr;
	if (ptr) {
		action_add_ref_count act;
		auto ptrptr = new (detail::agas_entry_t*)(reinterpret_cast<detail::agas_entry_t*>(ptr));
		int remote_rank = rank;
		set_address(ptrptr, [=](detail::agas_entry_t **ptrptr) {
			id_type remote;
			remote.set_rank(remote_rank);
			act(remote, reinterpret_cast<std::uintptr_t>(*ptrptr), -1);
			delete ptrptr;
		});
	} else {
		address = nullptr;
	}
}

template<class Archive>
void id_type::save(Archive &arc, const unsigned) const {
	arc << rank;
	if (address != nullptr) {
		auto ptr = reinterpret_cast<std::uintptr_t>(*(address));
		arc << ptr;
		action_add_ref_count act;
		id_type remote;
		remote.set_rank(rank);
		act(remote, ptr, +1);
	} else {
		arc << std::uintptr_t(0);
	}
}

template<class Function, class ... Args>
future<typename Function::return_type> async(const id_type &id, Args &&... args) {
	Function function;
	future<typename Function::return_type> future;
	if (id.get_rank() == hpx::detail::mpi_comm_rank()) {
		future = hpx::async(function, id, std::forward<Args>(args)...);
	} else {
		future = function.get_action_future(id, std::forward<Args>(args)...);
	}
	return future;
}

namespace lcos {
template<class Function, class ... Args>
future<typename std::enable_if<std::is_void<typename Function::return_type>::value, void>::type> broadcast(const std::vector<id_type> &ids, Args &&... args) {
	std::vector<hpx::future<void>> futs;
	for (int i = 0; i < ids.size(); i++) {
		futs.push_back(hpx::async<Function>(ids[i], std::forward<Args>(args)...));
	}
	return hpx::when_all(futs.begin(), futs.end()).then([](hpx::future<std::vector<hpx::future<void>>> fut) {
		auto futs = fut.get();
		hpx::wait_all(futs.begin(), futs.end());
	});
}

template<class Function, class ... Args>
future<std::vector<typename std::enable_if<!std::is_void<typename Function::return_type>::value, typename Function::return_type>::type>> broadcast(
		const std::vector<id_type> &ids, Args &&... args) {
	using type = typename Function::return_type;
	std::vector < hpx::future < type >> futs;
	for (int i = 0; i < ids.size(); i++) {
		futs.push_back(hpx::async<Function>(ids[i], std::forward<Args>(args)...));
	}
	return hpx::when_all(futs.begin(), futs.end()).then([](hpx::future<std::vector<hpx::future<type>>> fut) {
		auto futs = fut.get();
		std::vector<type> res(futs.size());
		for (int i = 0; i < res.size(); i++) {
			res[i] = futs[i].get();
		}
		return res;
	});
}

}

}

#endif /* ID_TYPE_IMPL_HPP_ */
