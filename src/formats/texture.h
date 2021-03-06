/*
	wrench - A set of modding tools for the Ratchet & Clank PS2 games.
	Copyright (C) 2019-2020 chaoticgd

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

#ifndef FORMATS_TEXTURE_H
#define FORMATS_TEXTURE_H

#include <array>
#include <vector>
#include <stdint.h>
#include <glm/glm.hpp>

#include "../stream.h"
#include "../gl_includes.h"

# /*
#	Stream-backed indexed texture.
# */

struct colour {
	uint8_t r, g, b, a;
};

struct vec2i {
	std::size_t x, y;
	
	bool operator==(vec2i& rhs) {
		return x == rhs.x && y == rhs.y;
	}
	
	bool operator!=(vec2i& rhs) {
		return x != rhs.x || y != rhs.y;
	}
};

class texture {
public:
	texture(
		stream* pixel_backing,
		std::size_t pixel_data_offset,
		stream* palette_backing,
		std::size_t palette_offset,
		vec2i size);
	texture(const texture&) = delete;
	texture(texture&&) = default;

	vec2i size() const;

	std::array<colour, 256> palette() const;
	void set_palette(std::array<colour, 256> palette_);

	std::vector<uint8_t> pixel_data() const;
	void set_pixel_data(std::vector<uint8_t> pixel_data_);

	std::string palette_path() const;
	std::string pixel_data_path() const;
	
#ifdef WRENCH_EDITOR
	void upload_to_opengl();
	GLuint opengl_id() const;
#else
	// Dummy to get the randomiser linking.
	void upload_to_opengl() {}
#endif
	
	std::string name;
	
private:
	stream* _pixel_backing;
	std::size_t _pixel_data_offset;
	stream* _palette_backing;
	std::size_t _palette_offset;
	vec2i _size;
#ifdef WRENCH_EDITOR
	gl_texture _opengl_texture;
#endif
};

// Won't affect the position indicator of backing.
std::optional<texture> create_fip_texture(stream* backing, std::size_t offset);

// Read a list of textures in the following format:
//  uint32_t count;
//  uint32_t offsets[count];
//  ... PIF textures ...
std::vector<texture> read_pif_list(stream* backing, std::size_t offset);

#endif
