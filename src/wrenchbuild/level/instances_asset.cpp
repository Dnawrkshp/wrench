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

#include "instances_asset.h"

#include <wrenchbuild/asset_unpacker.h>
#include <wrenchbuild/asset_packer.h>
#include <wrenchbuild/tests.h>

static void unpack_instances_asset(InstancesAsset& dest, InputStream& src, BuildConfig config, const char* hint);
static void unpack_help_messages(LevelWadAsset& dest, const HelpMessages& src, BuildConfig config, const char* hint);
static bool test_instances_asset(std::vector<u8>& src, AssetType type, BuildConfig config, const char* hint, AssetTestMode mode);
static const std::vector<GameplayBlockDescription>* get_gameplay_block_descriptions(Game game, const char* hint);

on_load(Instances, []() {
	InstancesAsset::funcs.unpack_rac1 = wrap_hint_unpacker_func<InstancesAsset>(unpack_instances_asset);
	InstancesAsset::funcs.unpack_rac2 = wrap_hint_unpacker_func<InstancesAsset>(unpack_instances_asset);
	InstancesAsset::funcs.unpack_rac3 = wrap_hint_unpacker_func<InstancesAsset>(unpack_instances_asset);
	InstancesAsset::funcs.unpack_dl = wrap_hint_unpacker_func<InstancesAsset>(unpack_instances_asset);
	
	InstancesAsset::funcs.test_rac = new AssetTestFunc(test_instances_asset);
	InstancesAsset::funcs.test_gc = new AssetTestFunc(test_instances_asset);
	InstancesAsset::funcs.test_uya = new AssetTestFunc(test_instances_asset);
	InstancesAsset::funcs.test_dl = new AssetTestFunc(test_instances_asset);
})

static void unpack_instances_asset(InstancesAsset& dest, InputStream& src, BuildConfig config, const char* hint) {
	std::vector<u8> buffer = src.read_multiple<u8>(0, src.size());
	unpack_instances(dest, nullptr, buffer, nullptr, config, hint);
}

void unpack_instances(InstancesAsset& dest, LevelWadAsset* help_dest, const std::vector<u8>& main, const std::vector<u8>* art, BuildConfig config, const char* hint) {
	Gameplay gameplay;
	read_gameplay(gameplay, main, config.game(), *get_gameplay_block_descriptions(config.game(), hint));
	
	Instances instances;
	HelpMessages help;
	std::map<s32, std::string> pvar_types;
	move_gameplay_to_instances(instances, &help, nullptr, pvar_types, gameplay, config.game());
	
	if(art) {
		Gameplay art_instances;
		read_gameplay(art_instances, *art, config.game(), DL_ART_INSTANCE_BLOCKS);
		instances.dir_lights = std::move(opt_iterator(art_instances.dir_lights));
		instances.tie_instances = std::move(opt_iterator(art_instances.tie_instances));
		instances.tie_groups = std::move(opt_iterator(art_instances.tie_groups));
		instances.shrub_instances = std::move(opt_iterator(art_instances.shrub_instances));
		instances.shrub_groups = std::move(opt_iterator(art_instances.shrub_groups));
	}
	
	std::string text = write_instances(instances);
	FileReference ref = dest.file().write_text_file(stringf("%s.instances", hint), text.c_str());
	dest.set_src(ref);
	
	if(help_dest) {
		unpack_help_messages(*help_dest, help, config, hint);
	}
	
	// Write types.
	AssetFile& build_file = dest.bank().asset_file("build.asset");
	for(auto& [class_id, cpp] : pvar_types) {
		std::string header_path = stringf("src/game_%s/update/moby%d.h", game_to_string(config.game()).c_str(), class_id);
		build_file.write_text_file(fs::path(header_path), cpp.c_str());
	}
}

static void unpack_help_messages(LevelWadAsset& dest, const HelpMessages& src, BuildConfig config, const char* hint) {
	if(src.us_english.has_value()) {
		MemoryInputStream us_english_stream(*src.us_english);
		unpack_asset_impl(dest.us_english_help_messages(), us_english_stream, nullptr, config, hint);
	}
	if(src.uk_english.has_value()) {
		MemoryInputStream uk_english_stream(*src.uk_english);
		unpack_asset_impl(dest.uk_english_help_messages(), uk_english_stream, nullptr, config, hint);
	}
	if(src.french.has_value()) {
		MemoryInputStream french_stream(*src.french);
		unpack_asset_impl(dest.french_help_messages(), french_stream, nullptr, config, hint);
	}
	if(src.german.has_value()) {
		MemoryInputStream german_stream(*src.german);
		unpack_asset_impl(dest.german_help_messages(), german_stream, nullptr, config, hint);
	}
	if(src.spanish.has_value()) {
		MemoryInputStream spanish_stream(*src.spanish);
		unpack_asset_impl(dest.spanish_help_messages(), spanish_stream, nullptr, config, hint);
	}
	if(src.italian.has_value()) {
		MemoryInputStream italian_stream(*src.italian);
		unpack_asset_impl(dest.italian_help_messages(), italian_stream, nullptr, config, hint);
	}
	if(src.japanese.has_value()) {
		MemoryInputStream japanese_stream(*src.japanese);
		unpack_asset_impl(dest.japanese_help_messages(), japanese_stream, nullptr, config, hint);
	}
	if(src.korean.has_value()) {
		MemoryInputStream korean_stream(*src.korean);
		unpack_asset_impl(dest.korean_help_messages(), korean_stream, nullptr, config, hint);
	}
}

Gameplay load_instances(const Asset& src, const BuildConfig& config, const char* hint) {
	std::vector <u8> gameplay_buffer;
	if(const InstancesAsset* asset = src.maybe_as<InstancesAsset>()) {
		std::unique_ptr<InputStream> gameplay_stream = asset->file().open_binary_file_for_reading(asset->src());
		gameplay_buffer = gameplay_stream->read_multiple<u8>(gameplay_stream->size());
	} else if(const BinaryAsset* asset = src.maybe_as<BinaryAsset>()) {
		std::unique_ptr<InputStream> gameplay_stream = asset->file().open_binary_file_for_reading(asset->src());
		gameplay_buffer = gameplay_stream->read_multiple<u8>(gameplay_stream->size());
	}
	Gameplay gameplay;
	read_gameplay(gameplay, gameplay_buffer, config.game(), *gameplay_block_descriptions_from_game(config.game()));
	return gameplay;
}

static bool test_instances_asset(std::vector<u8>& src, AssetType type, BuildConfig config, const char* hint, AssetTestMode mode) {
	const std::vector<GameplayBlockDescription>* blocks = get_gameplay_block_descriptions(config.game(), hint);
	
	// Parse gameplay file.
	Gameplay gameplay_in;
	read_gameplay(gameplay_in, src, config.game(), *blocks);
	
	// Separate out the different parts of the file.
	Instances instances_in;
	HelpMessages help_messages;
	OcclusionMappings occlusion;
	std::map<s32, std::string> pvar_types;
	move_gameplay_to_instances(instances_in, &help_messages, &occlusion, pvar_types, gameplay_in, config.game());
	
	// Write out instances file and read it back.
	std::string instances_text = write_instances(instances_in);
	write_file("/tmp/instances.txt", instances_text);
	Instances instances_out = read_instances(instances_text);
	instances_out.global_pvar = std::move(instances_in.global_pvar);
	
	// Write out new gameplay file.
	Gameplay gameplay_out;
	move_instances_to_gameplay(gameplay_out, instances_out, &help_messages, &occlusion);
	gameplay_out.pvar_table = gameplay_in.pvar_table;
	gameplay_out.pvar_scratchpad = gameplay_in.pvar_scratchpad;
	gameplay_out.pvar_relatives = gameplay_in.pvar_relatives;
	gameplay_out.global_pvar = gameplay_in.global_pvar;
	gameplay_out.global_pvar_table = gameplay_in.global_pvar_table;
	std::vector<u8> dest = write_gameplay(gameplay_out, config.game(), *blocks);
		
	// Compare the new file against the original.
	strip_trailing_padding_from_lhs(src, dest);
	bool headers_equal = diff_buffers(src, dest, 0, 0x100, mode == AssetTestMode::PRINT_DIFF_ON_FAIL);
	bool data_equal = diff_buffers(src, dest, 0x100, DIFF_REST_OF_BUFFER, mode == AssetTestMode::PRINT_DIFF_ON_FAIL);
	return headers_equal && data_equal;
}

static const std::vector<GameplayBlockDescription>* get_gameplay_block_descriptions(Game game, const char* hint) {
	const std::vector<GameplayBlockDescription>* blocks = nullptr;
	switch(game) {
		case Game::RAC:
			blocks = &RAC_GAMEPLAY_BLOCKS;
			break;
		case Game::GC:
		case Game::UYA:
			blocks = &GC_UYA_GAMEPLAY_BLOCKS;
			break;
		case Game::DL:
			if(strcmp(hint, FMT_INSTANCES_GAMEPLAY) == 0) {
				blocks = &DL_GAMEPLAY_CORE_BLOCKS;
			} else if(strcmp(hint, FMT_INSTANCES_ART) == 0) {
				blocks = &DL_ART_INSTANCE_BLOCKS;
			} else if(strcmp(hint, FMT_INSTANCES_MISSION) == 0) {
				blocks = &DL_GAMEPLAY_MISSION_INSTANCE_BLOCKS;
			} else {
				verify_not_reached("Invalid hint. Must be '%s', '%s' or '%s'.",
					FMT_INSTANCES_GAMEPLAY, FMT_INSTANCES_ART, FMT_INSTANCES_MISSION);
			}
			break;
		default: verify_fatal(0);
	}
	return blocks;
}
