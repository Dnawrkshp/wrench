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

#ifndef LZ_COMPRESSION_H
#define LZ_COMPRESSION_H

#include "../stream.h"

#include <map>
#include <cstring>
#include <utility>

// Decompress and recompress WAD segments used by the games to store various
// assets. Not to be confused with WAD archives.

packed_struct(wad_header,
	char magic[3]; // "WAD"
	uint32_t total_size; // Including header.
	uint8_t pad[9];
)

// Check the magic bytes.
bool validate_wad(char* magic);

// Throws stream_io_error, stream_format_error.
void decompress_wad(array_stream& dest, array_stream& src);

void compress_wad(array_stream& dest, array_stream& src, int thread_count);

#endif
