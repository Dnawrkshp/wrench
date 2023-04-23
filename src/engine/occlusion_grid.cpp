/*
	wrench - A set of modding tools for the Ratchet & Clank PS2 games.
	Copyright (C) 2019-2023 chaoticgd

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

#include "occlusion_grid.h"

#include <core/algorithm.h>

std::vector<OcclusionOctant> read_occlusion_grid(Buffer src) {
	ERROR_CONTEXT("reading occlusion grid");
	
	std::vector<OcclusionOctant> octants;
	
	s32 masks_offset = src.read<s32>(0);
	u16 z_coord = src.read<u16>(4);
	u16 z_count = src.read<u16>(6);
	auto z_offsets = src.read_multiple<u16>(8, z_count, "z offsets");
	for(s32 i = 0; i < z_count; i++) {
		if(z_offsets[i] == 0) {
			continue;
		}
		s32 z_offset = z_offsets[i] * 4;
		u16 y_coord = src.read<u16>(z_offset + 0);
		u16 y_count = src.read<u16>(z_offset + 2);
		auto y_offsets = src.read_multiple<u16>(z_offset + 4, y_count, "y offsets");
		for(s32 j = 0; j < y_count; j++) {
			if(y_offsets[j] == 0) {
				continue;
			}
			s32 x_offset = y_offsets[j] * 4;
			u16 x_coord = src.read<u16>(x_offset + 0);
			u16 x_count = src.read<u16>(x_offset + 2);
			auto x_indices = src.read_multiple<u16>(x_offset + 4, x_count, "x offsets");
			for(s32 k = 0; k < x_count; k++) {
				if(x_indices[k] == 0xffff) {
					continue;
				}
				
				OcclusionOctant& octant = octants.emplace_back();
				octant.x = x_coord + k;
				octant.y = y_coord + j;
				octant.z = z_coord + i;
				octant.index = x_indices[k];
				
				auto mask = src.read_multiple<u8>(masks_offset + x_indices[k] * 128, 128, "octant mask");
				verify_fatal(mask.size() == sizeof(octant.mask));
				memcpy(octant.mask, mask.lo, sizeof(octant.mask));
			}
		}
	}
	
	// This is to try and support multiple octants referencing a single mask.
	std::sort(BEGIN_END(octants), [&](auto& lhs, auto& rhs) { return lhs.index < rhs.index; });
	
	return octants;
}

void write_occlusion_grid(OutBuffer dest, std::vector<OcclusionOctant>& octants) {
	ERROR_CONTEXT("writing occlusion grid");
	
	s64 begin_offset = dest.alloc<s32>();
	
	// Mark duplicate masks.
	mark_duplicates(octants,
		[](OcclusionOctant& lhs, OcclusionOctant& rhs) { return memcmp(lhs.mask, rhs.mask, sizeof(OcclusionOctant::mask)); },
		[&](s32 index, s32 canonical) { octants[index].index = canonical; });
	
	for(s32 i = 0; i < (s32) octants.size(); i++) {
		octants[i].index2 = i;
	}
	
	std::stable_sort(BEGIN_END(octants), [&](auto& lhs, auto& rhs) { return lhs.x < rhs.x; });
	std::stable_sort(BEGIN_END(octants), [&](auto& lhs, auto& rhs) { return lhs.y < rhs.y; });
	std::stable_sort(BEGIN_END(octants), [&](auto& lhs, auto& rhs) { return lhs.z < rhs.z; });
	
	// Write out the tree.
	if(!octants.empty()) {
		u16 z_coord = checked_int_cast<u16>(octants.front().z);
		u16 z_count = checked_int_cast<u16>(octants.back().z - z_coord + 1);
		
		std::vector<u16> z_offsets(z_count, 0);
		
		printf("Z offsets = %lx\n", dest.tell());
		
		// Allocate the Z offsets.
		dest.pad(4);
		dest.write(z_coord);
		dest.write(z_count);
		dest.alloc_multiple<u16>(z_count);
		
		s32 i_start, j_start;
		
		printf("Y offsets = %lx\n", dest.tell());
		
		// Allocate the Y offsets, write out the Z offsets.
		i_start = 0;
		for(s32 i = 0; i < (s32) octants.size(); i++) {
			if(i == (s32) octants.size() - 1 || octants[i].z != octants[i + 1].z) {
				u16 y_coord = checked_int_cast<u16>(octants[i_start].y);
				u16 y_count = checked_int_cast<u16>(octants[i].y - y_coord + 1);
				
				dest.pad(4);
				u16 offset = checked_int_cast<u16>(dest.tell() / 4);
				
				dest.write(y_coord);
				dest.write(y_count);
				dest.alloc_multiple<u16>(y_count);
				
				// Write out Z offset.
				dest.write(begin_offset + 8 + (octants[i].z - z_coord) * 2, offset);
				z_offsets.at(octants[i].z - z_coord) = offset;
				
				i_start = i + 1;
			}
		}
		
		printf("Z offsets = %lx\n", dest.tell());
		
		// Write out the X offsets and Y offsets.
		i_start = 0;
		for(s32 i = 0; i < (s32) octants.size(); i++) {
			if(i == (s32) octants.size() - 1 || octants[i].z != octants[i + 1].z) {
				u16 y_coord = checked_int_cast<u16>(octants[i_start].y);
				
				j_start = i_start;
				for(s32 j = i_start; j <= i; j++) {
					if(j == i || octants[j].y != octants[j + 1].y) {
						u16 x_coord = checked_int_cast<u16>(octants[j_start].x);
						u16 x_count = checked_int_cast<u16>(octants[j].x - x_coord + 1);
						
						dest.pad(4);
						u16 offset = checked_int_cast<u16>(dest.tell() / 4);
						dest.write(x_coord);
						dest.write(x_count);
						s64 offsets = dest.alloc_multiple<u16>(x_count, 0xff);
						
						// Fill in mask indices.
						for(s32 k = j_start; k <= j; k++) {
							dest.write<u16>(offsets + (octants[k].x - x_coord) * 2, octants[k].index2);
						}
						
						// Write out Y offset.
						verify_fatal(z_offsets.at(octants[i].z - z_coord) != 0);
						dest.write(z_offsets.at(octants[j].z - z_coord) * 4 + 4 + (octants[j].y - y_coord) * 2, offset);
						
						j_start = j + 1;
					}
				}
				
				i_start = i + 1;
			}
		}
	}
	
	printf("Masks = %lx\n", dest.tell());
	
	std::stable_sort(BEGIN_END(octants), [&](auto& lhs, auto& rhs) { return lhs.index2 < rhs.index2; });
	
	// Write out the masks.
	dest.pad(0x10, 0);
	s64 masks_offset = dest.tell();
	for(s32 i = 0; i < (s32) octants.size(); i++) {
		if(octants[i].index == octants[i].index2) {
			dest.vec.insert(dest.vec.end(), octants[i].mask, octants[i].mask + sizeof(OcclusionOctant::mask));
		}
	}
	
	// Fill in pointer to masks.
	dest.write(begin_offset, (s32) masks_offset);
}

std::vector<OcclusionVector> read_occlusion_octants(std::string& str) {
	std::vector<OcclusionVector> octants;
	
	char line[128];
	
	const char* ptr = str.c_str();
	while(*ptr != '\0') {
		OcclusionVector& octant = octants.emplace_back();
		int total_read = sscanf(ptr, "%d %d %d\n", &octant.x, &octant.y, &octant.z);
		verify(octant.x != -1 && octant.y != -1 && octant.z != -1, "Failed to parse octants list.");
		ptr += total_read;
	}
	
	return octants;
}

void write_occlusion_octants(OutBuffer dest, const std::vector<OcclusionVector>& octants) {
	for(const OcclusionVector& octant : octants) {
		dest.writelf("%d %d %d", octant.x, octant.y, octant.z);
	}
}

void swap_occlusion(std::vector<OcclusionOctant>& grid, std::vector<OcclusionVector>& vectors) {
	verify_fatal(grid.size() == vectors.size());
	for(size_t i = 0; i < grid.size(); i++) {
		OcclusionOctant& octant = grid[i];
		OcclusionVector& vector = vectors[i];
		std::swap(octant.x, vector.x);
		std::swap(octant.y, vector.y);
		std::swap(octant.z, vector.z);
	}
}
