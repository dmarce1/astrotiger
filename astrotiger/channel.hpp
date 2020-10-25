#pragma once

#include <astrotiger/hpx.hpp>
#include <queue>

template<class T>
class channel {
	std::queue<T> q;
	mutex_type mtx;

public:
	void put(T &&data) {
		std::lock_guard<mutex_type> lock(mtx);
		q.push(std::move(data));
	}
	T get() {
		bool ready = false;
		do {
			std::unique_lock<mutex_type> lock(mtx);
			if (!q.empty()) {
				ready = true;
				break;
			} else {
				lock.unlock();
				hpx::this_thread::yield();
			}
		} while (!ready);
		T data = std::move(q.front());
		q.pop();
		return data;
	}
};
