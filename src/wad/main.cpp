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

#include <fstream>

#include <core/util.h>
#include <core/timer.h>
#include <assetmgr/asset.h>
#include "moby.h"
#include "tests.h"
#include "wad_file.h"
#include "collision.h"
#include "global_wads.h"

static void extract(fs::path input_path, fs::path output_path);
static void build(fs::path input_path, fs::path output_path);
static void extract_collision(fs::path input_path, fs::path output_path);
static void build_collision(fs::path input_path, fs::path output_path);
static void extract_moby(const char* input_path, const char* output_path);
static void build_moby(const char* input_path, const char* output_path);
static void print_usage(char* argv0);

#define require_args(arg_count) verify(argc == arg_count, "Incorrect number of arguments.");

int main(int argc, char** argv) {
	if(argc < 2) {
		print_usage(argv[0]);
		return 1;
	}
	
	std::string mode = argv[1];
	
	if(mode == "unpack") {
		require_args(4);
		AssetForest forest;
		AssetPack& pack = forest.mount<LooseAssetPack>(argv[3]);
		FILE* src = fopen(argv[2], "rb");
		std::vector<u8> header_bytes = read_file(src, 0, 0x800);
		unpack_misc_wad(pack, src, header_bytes);
		pack.write();
	} if(mode == "extract") {
		require_args(4);
		extract(argv[2], argv[3]);
	} else if(mode == "build") {
		require_args(4);
		build(argv[2], argv[3]);
	} else if(mode == "extract_collision") {
		require_args(4);
		extract_collision(argv[2], argv[3]);
	} else if(mode == "build_collision") {
		require_args(4);
		build_collision(argv[2], argv[3]);
	} else if(mode == "extract_moby") {
		require_args(4);
		extract_moby(argv[2], argv[3]);
	} else if(mode == "build_moby") {
		require_args(4);
		build_moby(argv[2], argv[3]);
	} else if(mode == "test") {
		require_args(3);
		run_tests(argv[2]);
	} else {
		print_usage(argv[0]);
		return 1;
	}
	return 0;
}

static void extract(fs::path input_path, fs::path output_path) {
	start_timer("Extract WAD");
	FILE* wad_file = fopen(input_path.string().c_str(), "rb");
	verify(wad_file, "Failed to open input file.");
	std::unique_ptr<Wad> wad = read_wad(wad_file);
	write_wad_json(output_path, wad.get());
	fclose(wad_file);
	stop_timer();
}

static void build(fs::path input_path, fs::path output_path) {
	start_timer("Building WAD");
	start_timer("Collecting files");
	std::unique_ptr<Wad> wad = read_wad_json(input_path);
	stop_timer();
	FILE* wad_file = fopen(output_path.string().c_str(), "wb");
	verify(wad_file, "Failed to open output file for writing.");
	write_wad(wad_file, wad.get());
	fclose(wad_file);
	stop_timer();
}

static void extract_collision(fs::path input_path, fs::path output_path) {
	auto collision = read_file(input_path);
	write_file("/", output_path, write_collada(read_collision(collision)), "w");
}

static void build_collision(fs::path input_path, fs::path output_path) {
	auto collision = read_file(input_path, "r");
	std::vector<u8> bin;
	write_collision(bin, read_collada(collision));
	write_file("/", output_path, bin);
}

static void extract_moby(const char* input_path, const char* output_path) {
	auto bin = read_file(input_path);
	MobyClassData moby = read_moby_class(bin, Game::RAC2);
	ColladaScene scene = recover_moby_class(moby, 0, 0);
	auto xml = write_collada(scene);
	write_file("/", output_path, xml, "w");
}

static void build_moby(const char* input_path, const char* output_path) {
	auto xml = read_file(input_path, "r");
	ColladaScene scene = read_collada(xml);
	MobyClassData moby = build_moby_class(scene);
	std::vector<u8> buffer;
	write_moby_class(buffer, moby, Game::RAC2);
	write_file("/", output_path, buffer);
}

static void print_usage(char* argv0) {
	printf("~* The Wrench WAD Utility *~\n");
	printf("\n");
	printf("An asset packer/unpacker for the Ratchet & Clank PS2 games intended for modding.\n");
	printf("\n");
	printf("usage: \n");
	printf("  %s extract <input wad> <output dir>\n", argv0);
	printf("  %s build <input level json> <output wad>\n", argv0);
	printf("  %s extract_collision <input bin> <output dae>\n", argv0);
	printf("  %s build_collision <input dae> <output bin>\n", argv0);
	printf("  %s test <level wads dir>\n", argv0);
}
