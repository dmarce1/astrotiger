#pragma once

#include <astrotiger/hpx.hpp>
#include <unordered_map>

#ifdef HPX_LITE

template<class T>
class channel {
	std::unordered_map<int, T> q;
	mutex_type mtx;

public:
	void put(T &&data, int step) {
		std::lock_guard<mutex_type> lock(mtx);
		q[step] = std::move(data);
	}
	hpx::future<T> get(int step) {
		return hpx::async([step, this]() {
			std::unique_lock<mutex_type> lock(mtx);
			while (1) {
				if (q.find(step) != q.end()) {
					auto data = std::move(q[step]);
					q.erase(step);
					return std::move(data);
				}
				lock.unlock();
				hpx::this_thread::yield();
				lock.lock();
			}
		});
	}
};

#else

#include <hpx/local_lcos/receive_buffer.hpp>

#include <boost/atomic.hpp>

#include <atomic>
#include <cstddef>

template<class T>
class channel {
private:
	hpx::lcos::local::receive_buffer<T> buffer;
public:
	channel() {
	}
	~channel() = default;
	channel(const channel&) = delete;
	channel(channel &&other) = delete;
	channel& operator=(channel &&other) = delete;
	channel& operator=(const channel &other) = delete;

	void put(T value, std::size_t cycle) {
		buffer.store_received(cycle, std::move(value));
	}

	hpx::future<T> get(std::size_t cycle) {
		return buffer.receive(cycle);
	}
};

#endif
