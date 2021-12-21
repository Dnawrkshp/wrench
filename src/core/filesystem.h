/*
	wrench - A set of modding tools for the Ratchet & Clank PS2 games.
	Copyright (C) 2019-2021 chaoticgd

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef CORE_FILESYSTEM_H
#define CORE_FILESYSTEM_H

#include <filesystem>
namespace fs = std::filesystem;

#include "buffer.h"

// These functions all call exit on error.
std::vector<u8> read_file(FILE* file, s64 offset, s64 size);
std::vector<u8> read_file(fs::path path, const char* open_mode = "rb");
std::string write_file(fs::path dest_dir, fs::path rel_path, Buffer buffer, const char* open_mode = "wb");
void extract_file(fs::path dest_path, FILE* dest, FILE* src, s64 offset, s64 size);

#endif
