#include "Util.hpp"

#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stacktrace>

#include <fenv.h>
#include <signal.h>

void enableFPE() {
	feclearexcept(FE_ALL_EXCEPT);
	feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}

void disableFPE() {
	feclearexcept(FE_ALL_EXCEPT);
	fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}

void stringToFile(std::string const &content, std::filesystem::path const &filePath) {
	namespace fs = std::filesystem;
	auto parentPath = filePath.parent_path();
	if (!parentPath.empty() && !fs::exists(parentPath)) {
		fs::create_directories(parentPath);
	}
	std::ofstream ofs(filePath, std::ios::out | std::ios::trunc);
	if (!ofs) {
		throw std::runtime_error("Failed to open file: " + filePath.string());
	}
	ofs << content;
}

#ifndef NDEBUG

void fpeHandler(int, siginfo_t*, void*) {
	std::cerr << "SIGFPE!\n" << std::stacktrace::current() << std::endl;
	std::_Exit(1);
}

void installFpeHandler() {
	struct sigaction sa;
	memset(&sa, 0, sizeof(sa));
	sa.sa_flags = SA_SIGINFO;
	sa.sa_sigaction = fpeHandler;
	sigemptyset(&sa.sa_mask);
	if (sigaction(SIGFPE, &sa, nullptr) != 0) {
		perror("sigaction");
		std::exit(1);
	}
}



#endif
