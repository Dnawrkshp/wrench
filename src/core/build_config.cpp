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

#include "build_config.h"

BuildConfig::BuildConfig(Game game, Region region)
	: value((u8) game | ((u8) region << 4)) {}

BuildConfig::BuildConfig(const std::string& game, const std::string& region)
	: BuildConfig(game_from_string(game), region_from_string(region)) {}

Game BuildConfig::game() const {
	return (Game) (value & 0xf);
}

Region BuildConfig::region() const {
	return (Region) (value >> 4);
}

Game game_from_string(const std::string& game) {
	if(game == "rac") return Game::RAC;
	if(game == "gc") return Game::GC;
	if(game == "uya") return Game::UYA;
	if(game == "dl") return Game::DL;
	return Game::UNKNOWN;
}

std::string game_to_string(Game game) {
	switch(game) {
		case Game::RAC: return "rac";
		case Game::GC: return "gc";
		case Game::UYA: return "uya";
		case Game::DL: return "dl";
		default: return "unknown";
	}
}

Region region_from_string(const std::string& region) {
	if(region == "us") return Region::US;
	if(region == "eu") return Region::EU;
	if(region == "japan") return Region::JAPAN;
}

std::string region_to_string(Region region) {
	switch(region) {
		case Region::US: return "us";
		case Region::EU: return "eu";
		case Region::JAPAN: return "japan";
		default: return "unknown";
	}
}
