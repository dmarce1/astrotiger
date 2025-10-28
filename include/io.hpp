/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
*******************************************************************************/
#pragma once

#include <filesystem>

void toFile(std::string const &content, std::filesystem::path const &filePath);

template<typename ...Args>
std::string print2string(char const *format, Args &&...args) {
	char *buffer;
	asprintf(&buffer, format, std::forward<Args>(args)...);
	std::string str(buffer);
	free(buffer);
	return str;
}

bool writeList(std::string const &, std::string const &, std::string const &);
