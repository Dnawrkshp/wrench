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

#ifndef WAD_GAMEPLAY_H
#define WAD_GAMEPLAY_H

#include "util.h"
#include "level.h"

struct GameplayBlockPos {
	s32 offset; // The offset of the pointer in the header.
	s32 order; // The ordering of the blocks in the file.
};

using GameplayBlockReadFunc = std::function<void(Gameplay& gameplay, Buffer src)>;
using GameplayBlockWriteFunc = std::function<void(OutBuffer dest, const Gameplay& gameplay)>;

struct GameplayBlockFuncs {
	GameplayBlockReadFunc read;
	GameplayBlockWriteFunc write;
};

struct GameplayBlockDescription {
	GameplayBlockPos rac1;
	GameplayBlockPos rac23;
	GameplayBlockPos rac4;
	GameplayBlockFuncs funcs;
	const char* name;
};

extern const std::vector<GameplayBlockDescription> gameplay_blocks;

void read_gameplay(Gameplay& gameplay, Buffer src);
std::vector<u8> write_gameplay(const Gameplay& gameplay);

#endif
