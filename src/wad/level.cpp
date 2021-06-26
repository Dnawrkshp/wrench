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

#include "level.h"

Json write_gameplay_json(Gameplay& gameplay) {
	Json json;
	
	json["metadata"] = Json {
		{"format", "gameplay"},
		{"format_version", 0},
		{"application", "Wrench WAD Extractor"},
		{"application_version", "0"}
	};
	
	Json data = to_json(gameplay);
	for(auto& object : data.items()) {
		json[object.key()] = object.value();
	}
	
	return json;
}
