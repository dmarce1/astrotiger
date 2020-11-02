#pragma once

#include <astrotiger/hpx.hpp>
#include <unordered_map>

template<class T>
class channel {
	std::unordered_map<int,T> q;
	mutex_type mtx;

public:
	void put(T &&data, int step) {
		std::lock_guard<mutex_type> lock(mtx);
		q[step] = std::move(data);
	}
	T get(int step) {
		std::unique_lock<mutex_type> lock(mtx);
		while(1) {
			if(q.find(step) != q.end()) {
				auto  data = std::move(q[step]);
				return std::move(data);
			}
			lock.unlock();
			hpx::this_thread::yield();
			lock.lock();
		}
	}
};
