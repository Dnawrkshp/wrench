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

#include "level.h"

level::level()
	: _history_index(0) { }

const texture_provider* level::get_texture_provider() const {
	return const_cast<level*>(this)->get_texture_provider();
}

std::vector<const point_object*> level::point_objects() const {
	std::vector<const point_object*> result;
	for(auto shrub : shrubs()) {
		result.push_back(shrub);
	}
	for(auto moby : mobies()) {
		result.push_back(moby.second);
	}
	return result;
}

std::vector<const shrub*> level::shrubs() const {
	auto result = const_cast<level*>(this)->shrubs();
	return std::vector<const shrub*>(result.begin(), result.end());
}

std::map<int32_t, const moby*> level::mobies() const {
	auto moby_map = const_cast<level*>(this)->mobies();
	std::map<int32_t, const moby*> result;
	for(auto& [uid, moby] : moby_map) {
		result.emplace(uid, moby);
	}
	return result;
}

bool level::is_selected(const game_object* obj) const {
	return std::find(selection.begin(), selection.end(), obj) != selection.end();
}

void level::undo() {
	if(_history_index <= 0) {
		throw command_error("Nothing to undo.");
	}
	_history_stack[--_history_index]->undo();
}

void level::redo() {
	if(_history_index >= _history_stack.size()) {
		throw command_error("Nothing to redo.");
	}
	_history_stack[_history_index++]->apply();
}
