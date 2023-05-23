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

#ifndef LEVEL_H
#define LEVEL_H

#include <map>

#include <core/util.h>
#include <core/collada.h>
#include <core/buffer.h>
#include <core/texture.h>
#include <core/filesystem.h>
#include <core/build_config.h>
#include <instancemgr/instance.h>

// A plane that defines the bounds of a chunk. Everything on the side of the
// plane in the direction that the normal is pointing is inside the chunk.
struct ChunkPlane {
	glm::vec3 point;
	glm::vec3 normal;
};

packed_struct(LevelSettingsThirdPart,
	s32 unknown_0;
	s32 unknown_4;
	s32 unknown_8;
	s32 unknown_c;
)

packed_struct(LevelSettingsRewardStats,
	f32 xp_decay_rate;
	f32 xp_decay_min;
	f32 bolt_decay_rate;
	f32 bolt_decay_min;
	s32 unknown_10;
	s32 unknown_14;
)

packed_struct(LevelSettingsFifthPart,
	/* 0x00 */ s32 unknown_0;
	/* 0x04 */ s32 moby_inst_count;
	/* 0x08 */ s32 unknown_8;
	/* 0x0c */ s32 unknown_c;
	/* 0x10 */ s32 unknown_10;
	/* 0x14 */ s32 dbg_hit_points;
)

struct LevelSettings {
	Opt<glm::vec3> background_colour;
	Opt<glm::vec3> fog_colour;
	f32 fog_near_dist = 0.f;
	f32 fog_far_dist = 0.f;
	f32 fog_near_intensity = 0.f;
	f32 fog_far_intensity = 0.f;
	f32 death_height = 0.f;
	bool is_spherical_world;
	glm::vec3 sphere_pos = {0.f, 0.f, 0.f};
	glm::vec3 ship_pos = {0.f, 0.f, 0.f};
	f32 ship_rot_z = 0.f;
	PathLink ship_path;
	CuboidLink ship_camera_cuboid_start;
	CuboidLink ship_camera_cuboid_end;
	// Planes specifying the volumes of the level chunks. The first element
	// represents the second chunk, and the second element represents the third
	// chunk. If both tests fail, you can assume it's the first chunk (chunk 0).
	std::vector<ChunkPlane> chunk_planes;
	Opt<s32> core_sounds_count;
	Opt<s32> rac3_third_part;
	Opt<std::vector<LevelSettingsThirdPart>> third_part;
	Opt<LevelSettingsRewardStats> reward_stats;
	Opt<LevelSettingsFifthPart> fifth_part;
	Opt<std::vector<u8>> dbg_attack_damage;
};

struct WtfNode;
struct WtfWriter;

LevelSettings read_level_settings(const WtfNode* node);
void write_level_settings(WtfWriter* ctx, const LevelSettings& settings);

enum PvarFieldDescriptor {
	PVAR_INTEGERS_BEGIN = 0,
	PVAR_S8 = 1, PVAR_S16 = 2, PVAR_S32 = 3,
	PVAR_U8 = 4, PVAR_U16 = 5, PVAR_U32 = 6,
	PVAR_INTEGERS_END = 7,
	PVAR_F32 = 8,
	PVAR_POINTERS_BEGIN = 100,
	PVAR_RUNTIME_POINTER = 101,
	PVAR_RELATIVE_POINTER = 102,
	PVAR_SCRATCHPAD_POINTER = 103,
	PVAR_GLOBAL_PVAR_POINTER = 104,
	PVAR_POINTERS_END = 105,
	PVAR_STRUCT = 106
};

std::string pvar_descriptor_to_string(PvarFieldDescriptor descriptor);
PvarFieldDescriptor pvar_string_to_descriptor(std::string str);

struct PvarField {
	s32 offset;
	std::string name;
	PvarFieldDescriptor descriptor = PVAR_U8;
	std::string value_type; // Only set for pointer types.
	
	s32 size() const;
	
	template <typename T>
	void enumerate_fields(T& t) {
		DEF_FIELD(offset);
		DEF_FIELD(name);
		std::string type = pvar_descriptor_to_string(descriptor);
		DEF_FIELD(type);
		descriptor = pvar_string_to_descriptor(type);
		if(descriptor > PVAR_POINTERS_BEGIN || descriptor < PVAR_POINTERS_END) {
			DEF_FIELD(value_type);
		}
	}
};

struct PvarType {
	std::vector<PvarField> fields;
	
	bool insert_field(PvarField to_insert, bool sort);
	
	template <typename T>
	void enumerate_fields(T& t) {
		DEF_FIELD(fields);
	}
};

struct PvarTypes {
	std::map<s32, PvarType> moby;
	std::map<s32, PvarType> camera;
	std::map<s32, PvarType> sound;
};

#endif
