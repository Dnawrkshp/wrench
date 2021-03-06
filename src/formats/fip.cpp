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

#include "fip.h"

#include <cstring>
#include <stdexcept>

#include "bmp.h"

bool validate_fip(char* magic) {
	return std::memcmp(magic, "2FIP", 4) == 0;
}

void fip_to_bmp(stream& dest, stream& src) {

	auto src_header = src.read<fip_header>(0);
	if(!validate_fip(src_header.magic)) {
		throw stream_format_error("Tried to read invalid FIP segment.");
	}

	bmp_file_header header;
	std::memcpy(header.magic, "BM", 2);
	header.pixel_data =
		sizeof(bmp_file_header) +
		sizeof(bmp_info_header) +
		sizeof(bmp_colour_table_entry) * 256;
	header.file_size =
		header.pixel_data.value +
		src_header.width * src_header.height;
	header.reserved = 1337;
	dest.write<bmp_file_header>(0, header);

	bmp_info_header info;
	info.info_header_size      = 40;
	info.width                 = src_header.width;
	info.height                = src_header.height;
	info.num_colour_planes     = 1;
	info.bits_per_pixel        = 8;
	info.compression_method    = 0;
	info.pixel_data_size       = info.width * info.height;
	info.horizontal_resolution = 0;
	info.vertical_resolution   = 0;
	info.num_colours           = 256;
	info.num_important_colours = 0;
	dest.write<bmp_info_header>(info);

	for(int i = 0; i < 256; i++) {
		auto src_pixel = src_header.palette[i];
		bmp_colour_table_entry pixel;
		pixel.b = src_pixel.b;
		pixel.g = src_pixel.g;
		pixel.r = src_pixel.r;
		pixel.pad = 0;
		dest.write<bmp_colour_table_entry>(pixel);
	}

	uint32_t row_size = ((info.bits_per_pixel * info.width + 31) / 32) * 4;
	std::size_t pixel_data = dest.tell();

	for(int y = info.height - 1; y >= 0; y--) {
		dest.seek(pixel_data + y * row_size);
		for(int x = 0; x < info.width; x++) {
			uint8_t palette_index = src.read<uint8_t>();
			dest.write<uint8_t>(decode_palette_index(palette_index));
		}
	}
}

void bmp_to_fip(stream& dest, stream& src) {
	auto file_header = src.read<bmp_file_header>(0);

	if(!validate_bmp(file_header)) {
		throw stream_format_error("Invalid BMP header.");
	}

	std::size_t secondary_header_offset = src.tell();
	auto info_header = src.read<bmp_info_header>();

	if(info_header.bits_per_pixel != 8) {
		throw stream_format_error("The BMP file must use indexed colour (with at most 256 colours).");
	}

	if(info_header.num_colours > 256) {
		throw stream_format_error("The BMP colour palette must contain at most 256 colours.");
	}

	fip_header header;
	std::memcpy(header.magic, "2FIP", 4);
	std::memset(header.unknown1, 0, sizeof(header.unknown1));
	header.width = info_header.width;
	header.height = info_header.height;
	std::memset(header.unknown2, 0, sizeof(header.unknown2));
	// Some BMP files have a larger header.
	src.seek(secondary_header_offset + info_header.info_header_size);
	uint32_t i;
	for(i = 0; i < info_header.num_colours; i++) {
		auto src_pixel = src.read<bmp_colour_table_entry>();
		auto& dest_pixel = header.palette[i];
		dest_pixel.r = src_pixel.r;
		dest_pixel.g = src_pixel.g;
		dest_pixel.b = src_pixel.b;
		dest_pixel.a = 0x80;
	}
	for(; i < 256; i++) {
		// Set unused palette entries to black.
		auto& dest_pixel = header.palette[i];
		dest_pixel.r = 0;
		dest_pixel.g = 0;
		dest_pixel.b = 0;
		dest_pixel.a = 0x80;
	}
	dest.write<fip_header>(0, header);

	uint32_t row_size = ((info_header.bits_per_pixel * info_header.width + 31) / 32) * 4;
	uint32_t pixel_data = file_header.pixel_data.value;

	for(int y = info_header.height - 1; y >= 0; y--) {
		src.seek(pixel_data + y * row_size);
		for(int x = 0; x < info_header.width; x++) {
			uint8_t palette_index = src.read<uint8_t>();
			dest.write<uint8_t>(decode_palette_index(palette_index));
		}
	}
}

uint8_t decode_palette_index(uint8_t index) {
	// Swap middle two bits
	//  e.g. 00010000 becomes 00001000.
	int a = index & 8;
	int b = index & 16;
	if(a && !b) {
		index += 8;
	} else if(!a && b) {
		index -= 8;
	}
	return index;
}
