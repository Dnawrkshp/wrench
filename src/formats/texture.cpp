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

#include "texture.h"

#include "../util.h" // int_to_hex
#include "../iso_stream.h"
#include "fip.h" // decode_palette_index

texture::texture(
	stream* pixel_backing,
	std::size_t pixel_data_offset,
	stream* palette_backing,
	std::size_t palette_offset,
	vec2i size)
	: _pixel_backing(pixel_backing),
	  _pixel_data_offset(pixel_data_offset),
	  _palette_backing(palette_backing),
	  _palette_offset(palette_offset),
	  _size(size) {}
	
vec2i texture::size() const {
	return _size;
}

std::array<colour, 256> texture::palette() const {
	char data[1024];
	_palette_backing->peek_n(data, _palette_offset, 1024);

	std::array<colour, 256> result;
	for(int i = 0; i < 256; i++) {
		result[decode_palette_index(i)] = {
			static_cast<uint8_t>(data[i * 4 + 0]),
			static_cast<uint8_t>(data[i * 4 + 1]),
			static_cast<uint8_t>(data[i * 4 + 2]),
			static_cast<uint8_t>(data[i * 4 + 3])
		};
	}
	return result;
}

void texture::set_palette(std::array<colour, 256> palette_) {
	std::array<char, 1024> colours;
	for(int i =0 ; i < 256; i++) {
		colour c = palette_[decode_palette_index(i)];
		colours[i * 4 + 0] = c.r;
		colours[i * 4 + 1] = c.g;
		colours[i * 4 + 2] = c.b;
		colours[i * 4 + 3] = c.a;
	}
	_palette_backing->seek(_palette_offset);
	_palette_backing->write_n(colours.data(), colours.size());
}

std::vector<uint8_t> texture::pixel_data() const {
	vec2i size_ = size();
	std::vector<uint8_t> result(size_.x * size_.y);
	_pixel_backing->peek_n(reinterpret_cast<char*>(result.data()), _pixel_data_offset, result.size());
	return result;
}

void texture::set_pixel_data(std::vector<uint8_t> pixel_data_) {
	_pixel_backing->seek(_pixel_data_offset);
	_pixel_backing->write_n(reinterpret_cast<char*>(pixel_data_.data()), pixel_data_.size());
}

std::string texture::palette_path() const {
	return _palette_backing->resource_path() + "+0x" + int_to_hex(_palette_offset);
}

std::string texture::pixel_data_path() const {
	return _pixel_backing->resource_path() + "+0x" + int_to_hex(_pixel_data_offset);
}

#ifdef WRENCH_EDITOR

void texture::upload_to_opengl() {
	std::vector<uint8_t> indexed_pixel_data = pixel_data();
	std::vector<uint8_t> colour_data(indexed_pixel_data.size() * 4);
	auto palette_data = palette();
	for(std::size_t i = 0; i < indexed_pixel_data.size(); i++) {
		colour c = palette_data[indexed_pixel_data[i]];
		colour_data[i * 4] = c.r;
		colour_data[i * 4 + 1] = c.g;
		colour_data[i * 4 + 2] = c.b;
		colour_data[i * 4 + 3] = static_cast<int>(c.a) * 2 - 1;
	}
	
	glDeleteTextures(1, &_opengl_texture());
	glGenTextures(1, &_opengl_texture());
	glBindTexture(GL_TEXTURE_2D, _opengl_texture());
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, _size.x, _size.y, 0, GL_RGBA, GL_UNSIGNED_BYTE, colour_data.data());
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
}

GLuint texture::opengl_id() const {
	return _opengl_texture();
}

#endif

std::optional<texture> create_fip_texture(stream* backing, std::size_t offset) {
	fip_header header = backing->peek<fip_header>(offset);
	if(!validate_fip(header.magic)) {
		return {};
	}
	vec2i size { header.width, header.height };
	std::size_t pixel_offset = offset + sizeof(fip_header);
	std::size_t palette_offset = offset + offsetof(fip_header, palette);
	return texture(backing, pixel_offset, backing, palette_offset, size);
}

std::vector<texture> read_pif_list(stream* backing, std::size_t offset) {
	uint32_t count = backing->read<uint32_t>(0);
	
	std::vector<uint32_t> offsets(count);
	backing->read_v(offsets);
	
	std::vector<texture> textures;
	for(uint32_t texture_offset : offsets) {
		auto texture = create_fip_texture(backing, offset + texture_offset);
		if(texture) {
			textures.emplace_back(std::move(*texture));
		}
	}
	
	return textures;
}
