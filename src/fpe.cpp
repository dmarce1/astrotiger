#include "fpe.hpp"

#include <cfenv>
#include <cstring>
#include <ostream>

void FpeEnabler::fpeHandler(int, siginfo_t *info, void*) noexcept {
	char const *reason = "Unknown FPE";
	switch (info->si_code) {
	case FPE_INTDIV:
		reason = "Integer divide by zero";
		break;
	case FPE_INTOVF:
		reason = "Integer overflow";
		break;
	case FPE_FLTDIV:
		reason = "Floating-point divide by zero";
		break;
	case FPE_FLTOVF:
		reason = "Floating-point overflow";
		break;
	case FPE_FLTUND:
		reason = "Floating-point underflow";
		break;
	case FPE_FLTRES:
		reason = "Floating-point inexact result";
		break;
	case FPE_FLTINV:
		reason = "Invalid floating-point operation";
		break;
	case FPE_FLTSUB:
		reason = "Subscript out of range";
		break;
	}
	write(STDERR_FILENO, reason, strlen(reason));
	write(STDERR_FILENO, "\n", 1);
	exit(EXIT_FAILURE);
}

FpeEnabler::FpeEnabler() noexcept :
		savedSignalAction_(), savedSignalSet_() {
	sigset_t signalSet;
	struct sigaction signalAction;
	sigemptyset(&signalSet);
	sigaddset(&signalSet, SIGFPE);
	pthread_sigmask(SIG_UNBLOCK, &signalSet, &savedSignalSet_);
	signalAction.sa_flags = SA_SIGINFO;
	signalAction.sa_sigaction = &FpeEnabler::fpeHandler;
	sigemptyset(&signalAction.sa_mask);
	feclearexcept(FE_ALL_EXCEPT);
	feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
	if (sigaction(SIGFPE, &signalAction, &savedSignalAction_) != 0) {
		perror("FpeEnabler failed to install");
		exit(1);
	}
}

FpeEnabler::~FpeEnabler() noexcept {
	fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
	pthread_sigmask(SIG_UNBLOCK, &savedSignalSet_, nullptr);
	if (sigaction(SIGFPE, &savedSignalAction_, nullptr)) {
		perror("FpeEnabler failed to install");
		exit(1);
	}
}

