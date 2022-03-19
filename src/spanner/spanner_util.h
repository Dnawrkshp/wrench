/*
	wrench - A set of modding tools for the Ratchet & Clank PS2 games.
	Copyright (C) 2019-2022 chaoticgd

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

#ifndef SPANNER_UTIL_H
#define SPANNER_UTIL_H

#include <assetmgr/asset.h>
#include <assetmgr/asset_types.h>

template <typename Range>
Asset& unpack_binary(Asset& parent, const FileHandle& src, Range range, const char* child, fs::path path) {
	std::vector<u8> bytes = src.read_binary(range.bytes());
	BinaryAsset& binary = parent.child<BinaryAsset>(child);
	if(range.bytes().size > 0) {
		binary.set_src(parent.file().write_binary_file(path, bytes));
	}
	return binary;
}

template <typename Range>
std::vector<Asset*> unpack_binaries(Asset& parent, const FileHandle& src, Range* ranges, s32 count, const char* child, const char* extension = ".bin") {
	fs::path path = fs::path(child)/child;
	CollectionAsset& collection = parent.asset_file(path).child<CollectionAsset>(child);
	
	std::vector<Asset*> assets;
	for(s32 i = 0; i < count; i++) {
		std::string name = std::to_string(i);
		assets.emplace_back(&unpack_binary(collection, src, ranges[i], name.c_str(), name + extension));
	}
	
	return assets;
}

Asset& unpack_binary_from_memory(Asset& parent, Buffer src, ByteRange range, const char* child, const char* extension = ".bin");
Asset& unpack_compressed_binary(Asset& parent, Buffer src, ByteRange range, const char* child, const char* extension = ".bin");
std::vector<Asset*> unpack_compressed_binaries(Asset& parent, Buffer src, ByteRange* ranges, s32 count, const char* child);

template <typename Header>
std::pair<FileHandle, Header> open_wad_file(BinaryAsset& asset) {
	FileHandle file = asset.file().open_binary_file_for_reading(asset.src());
	s32 header_size = Buffer(file.read_binary(ByteRange64{0, 4})).read<s32>(0, "header");
	std::vector<u8> header_bytes = file.read_binary(ByteRange64{0, header_size});
	Header header = Buffer(header_bytes).read<Header>(0, "file header");
	return {std::move(file), header};
}

#endif
