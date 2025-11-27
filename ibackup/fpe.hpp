#pragma once

#include <csignal>

struct FpeEnabler {
	FpeEnabler() noexcept;
	~FpeEnabler() noexcept;
private:
	static void fpeHandler(int, siginfo_t*, void*) noexcept;
	struct sigaction savedSignalAction_;
	sigset_t savedSignalSet_;
	static thread_local bool enabled;
};


#define enableFPE() FpeEnabler fpeEnabler{}
