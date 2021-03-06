/*
	wrench - A set of modding tools for the Ratchet & Clank PS2 games.
	Copyright (C) 2019 chaoticgd

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

#ifndef FORMATS_RACPAK_H
#define FORMATS_RACPAK_H

#include <vector>

#include "../stream.h"
#include "../iso_stream.h"
#include "wad.h"

# /*
#	Tools to open and modify racpak (*.WAD) archive files.
# */

struct racpak_entry {
	std::size_t offset;
	std::size_t size;
};

class racpak {
public:
	racpak(stream* backing, std::size_t offset, std::size_t size);

	std::size_t num_entries();
	std::size_t base();
	racpak_entry entry(std::size_t index);
	stream* open(racpak_entry file);
	bool is_compressed(racpak_entry entry);

private:
	proxy_stream _backing;
	std::size_t _base;
	std::vector<std::unique_ptr<proxy_stream>> _open_segments;
};

#endif
