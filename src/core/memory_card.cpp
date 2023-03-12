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

#include "memory_card.h"

namespace memory_card {

packed_struct(SaveSlotFileHeader,
	/* 0x0 */ s32 game_data_size;
	/* 0x4 */ s32 level_data_size;
)

packed_struct(ChecksumHeader,
	/* 0x0 */ s32 size;
	/* 0x4 */ s32 checksum;
)

packed_struct(SectionHeader,
	/* 0x0 */ s32 type;
	/* 0x4 */ s32 size;
)

File read(Buffer src, const fs::path& path) {
	File file;
	s64 pos = 0;
	
	file.path = path;
	file.type = identify(path.filename().string());
	
	switch(file.type) {
		case FileType::MAIN: {
			file.main.data = src.read_all<u8>().copy();
			break;
		}
		case FileType::NET: {
			file.net.sections = read_sections(&file.checksum_does_not_match, src, pos);
			break;
		}
		case FileType::PATCH: {
			file.patch.data = src.read_all<u8>().copy();
			break;
		}
		case FileType::SLOT: {
			const SaveSlotFileHeader& file_header = src.read<SaveSlotFileHeader>(pos, "file header");
			pos += sizeof(SaveSlotFileHeader);
			
			file.slot.sections = read_sections(&file.checksum_does_not_match, src, pos);
			while(pos + 3 < src.size()) {
				file.slot.levels.emplace_back(read_sections(&file.checksum_does_not_match, src, pos));
			}
			break;
		}
		case FileType::SYS: {
			file.sys.data = src.read_all<u8>().copy();
			break;
		}
	}
	
	return file;
}

FileType identify(std::string filename) {
	for(char& c : filename) c = tolower(c);
	if(filename.find("ratchet") != std::string::npos) {
		return FileType::MAIN;
	} if(filename.starts_with("net")) {
		return FileType::NET;
	} else if(filename.starts_with("patch")) {
		return FileType::PATCH;
	} else if(filename.starts_with("save")) {
		return FileType::SLOT;
	} else if(filename.starts_with("icon")) {
		return FileType::SYS;
	} else {
		verify_not_reached("Unable to identify file type.");
	}
}

std::vector<Section> read_sections(bool* checksum_does_not_match_out, Buffer src, s64& pos) {
	std::vector<Section> sections;
	
	const ChecksumHeader& checksum_header = src.read<ChecksumHeader>(pos, "checksum header");
	pos += sizeof(ChecksumHeader);
	
	u32 check_value = checksum(src.subbuf(pos, checksum_header.size));
	if(check_value != checksum_header.checksum) {
		*checksum_does_not_match_out |= true;
	}
	
	for(;;) {
		const SectionHeader& section_header = src.read<SectionHeader>(pos, "section header");
		pos += sizeof(SectionHeader);
		if(section_header.type == -1) {
			break;
		}
		
		// Preserve uninitialised padding.
		s32 read_size = align64(section_header.size, 4);
		
		Section& section = sections.emplace_back();
		section.offset = (s32) pos;
		section.type = section_header.type;
		section.unpadded_size = section_header.size;
		section.data = src.read_bytes(pos, read_size, "section data");
		pos += read_size;
	}
	
	return sections;
}

void write(OutBuffer dest, File& file) {
	switch(file.type) {
		case FileType::MAIN: {
			dest.write_multiple(file.main.data);
			break;
		}
		case FileType::NET: {
			write_sections(dest, file.net.sections);
			break;
		}
		case FileType::PATCH: {
			dest.write_multiple(file.patch.data);
			break;
		}
		case FileType::SLOT: {
			s64 file_header_ofs = dest.alloc<SaveSlotFileHeader>();
			SaveSlotFileHeader file_header;
			file_header.game_data_size = (s32) write_sections(dest, file.slot.sections);
			file_header.level_data_size = 0;
			for(std::vector<Section>& sections : file.slot.levels) {
				u32 data_size = (s32) write_sections(dest, sections);
				if(file_header.level_data_size == 0) {
					file_header.level_data_size = data_size;
				} else {
					assert(data_size == file_header.level_data_size);
				}
			}
			dest.write(file_header_ofs, file_header);
			break;
		}
		case FileType::SYS: {
			dest.write_multiple(file.sys.data);
			break;
		}
	}
}

s64 write_sections(OutBuffer dest, std::vector<Section>& sections) {
	s64 checksum_header_ofs = dest.alloc<ChecksumHeader>();
	s64 checksum_start_ofs = dest.tell();
	
	for(Section& section : sections) {
		SectionHeader header;
		header.type = section.type;
		header.size = section.unpadded_size;
		dest.write(header);
		section.offset = (s32) dest.tell();
		dest.write_multiple(section.data);
		dest.pad(4);
	}
	SectionHeader terminator;
	terminator.type = -1;
	terminator.size = 0;
	dest.write(terminator);
	
	s64 checksum_end_ofs = dest.tell();
	
	ChecksumHeader checksum_header;
	checksum_header.size = (s32) (checksum_end_ofs - checksum_start_ofs);
	checksum_header.checksum = checksum(Buffer(dest.vec).subbuf(checksum_start_ofs, checksum_header.size));
	dest.write(checksum_header_ofs, checksum_header);
	
	return dest.tell() - checksum_header_ofs;
}

u32 checksum(Buffer src) {
	const u8* ptr = src.lo;
	u32 value = 0xedb88320;
	for(const u8* ptr = src.lo; ptr < src.hi; ptr++) {
		value ^= (u32) *ptr << 8;
		for(s32 repeat = 0; repeat < 8; repeat++) {
			if((value & 0x8000) == 0) {
				value <<= 1;
			} else {
				value = (value << 1) ^ 0x1f45;
			}
		}
	}
	return value & 0xffff;
}

SaveGame parse(const File& file) {
	SaveGame save;
	switch(file.type) {
		case FileType::MAIN: break;
		case FileType::NET: save = parse_net(file); break;
		case FileType::PATCH: break;
		case FileType::SLOT: save = parse_slot(file); break;
		case FileType::SYS: break;
	}
	save.type = file.type;
	return save;
}

void update(File& dest, const SaveGame& save) {
	switch(dest.type) {
		case FileType::MAIN: break;
		case FileType::NET: update_net(dest, save); break;
		case FileType::PATCH: break;
		case FileType::SLOT: update_slot(dest, save); break;
		case FileType::SYS: break;
	}
}

template <typename T>
static void parse_section(const memory_card::Section& section, Opt<T>& dest) {
	ERROR_CONTEXT("%s section", section_type(section.type));
	s64 size_difference = section.data.size() - sizeof(T);
	verify(size_difference > -1 && size_difference < 4, "Section has unexpected size.");
	dest = Buffer(section.data).read<T>(0);
}

template <typename T>
static void update_section(OutBuffer dest, const Opt<T>& src) {
	if(src.has_value()) {
		dest.write(0, *src);
	}
}

template <typename T>
static void parse_section_array(const memory_card::Section& section, Opt<T>& dest) {
	ERROR_CONTEXT("%s section", section_type(section.type));
	s64 size_difference = section.data.size() - sizeof(T);
	verify(size_difference > -1 && size_difference < 4, "Array section has unexpected size.");
	dest = Buffer(section.data).read_multiple<typename T::value_type>(0, T::element_count, "array");
}

template <typename T>
static void update_section_array(OutBuffer dest, const Opt<T>& src) {
	if(src.has_value()) {
		dest.write_multiple(0, *src);
	}
}

SaveGame parse_net(const File& file) {
	assert(file.type == FileType::NET);
	SaveGame save;
	save.loaded = true;
	for(const Section& section : file.net.sections) {
		Buffer buffer = section.data;
		switch(section.type) {
			case ST_GAMEMODEOPTIONS: save.game_mode_options = buffer.read<GameModeStruct>(0);            break;
			case ST_MPPROFILES:      save.mp_profiles       = buffer.read_multiple<ProfileStruct>(0, 8); break;
		}
	}
	return save;
}

void update_net(File& dest, const SaveGame& save) {
	for(Section& section : dest.net.sections) {
		OutBuffer buffer = section.data;
		switch(section.type) {
			case ST_GAMEMODEOPTIONS: update_section      (buffer, save.game_mode_options); break;
			case ST_MPPROFILES:      update_section_array(buffer, save.mp_profiles);       break;
		}
	}
}

SaveGame parse_slot(const File& file) {
	assert(file.type == FileType::SLOT);
	SaveGame save;
	save.loaded = true;
	for(const Section& section : file.slot.sections) {
		switch(section.type) {
			case ST_LEVEL:                parse_section      (section, save.level);                     break;
			case ST_ELAPSEDTIME:          parse_section      (section, save.elapsed_time);              break;
			case ST_LASTSAVETIME:         parse_section      (section, save.last_save_time);            break;
			case ST_GLOBALFLAGS:          parse_section_array(section, save.global_flags);              break;
			case ST_CHEATSACTIVATED:      parse_section_array(section, save.cheats_activated);          break;
			case ST_SKILLPOINTS:          parse_section_array(section, save.skill_points);              break;
			case ST_HELPDATAMESSAGES:     parse_section_array(section, save.help_data_messages);        break;
			case ST_HELPDATAMISC:         parse_section_array(section, save.help_data_misc);            break;
			case ST_HELPDATAGADGETS:      parse_section_array(section, save.help_data_gadgets);         break;
			case ST_CHEATSEVERACTIVATED:  parse_section_array(section, save.cheats_ever_activated);     break;
			case ST_SETTINGS:             parse_section      (section, save.settings);                  break;
			case ST_HEROSAVE:             parse_section      (section, save.hero_save);                 break;
			case ST_MOVIESPLAYEDRECORD:   parse_section_array(section, save.movies_played_record);      break;
			case ST_TOTALPLAYTIME:        parse_section      (section, save.total_play_time);           break;
			case ST_TOTALDEATHS:          parse_section      (section, save.total_deaths);              break;
			case ST_HELPLOG:              parse_section_array(section, save.help_log);                  break;
			case ST_HELPLOGPOS:           parse_section      (section, save.help_log_pos);              break;
			case ST_HEROGADGETBOX:        parse_section      (section, save.hero_gadget_box);           break;
			case ST_PURCHASEABLEGADGETS:  parse_section      (section, save.purchaseable_gadgets);      break;
			case ST_PURCHASEABLEBOTUPGRD: parse_section      (section, save.purchaseable_bot_upgrades); break;
			case ST_PURCHASEABLEWRENCH:   parse_section      (section, save.purchaseable_wrench_level); break;
			case ST_PURCHASEABLEPOSTMODS: parse_section_array(section, save.purchaseable_post_fx_mods); break;
			case ST_BOTSAVE:              parse_section      (section, save.bot_save);                  break;
			case ST_FIRSTPERSONMODE:      parse_section_array(section, save.first_person_desired_mode); break;
			case ST_SAVEDDIFFICULTYLEVEL: parse_section      (section, save.saved_difficulty_level);    break;
			case ST_PLAYERSTATISTICS:     parse_section_array(section, save.player_statistics);         break;
			case ST_BATTLEDOMEWINSLOSSES: parse_section_array(section, save.battledome_wins_and_losses); break;
			case ST_ENEMYKILLS:           parse_section_array(section, save.enemy_kills);               break;
			case ST_QUICKSWITCHGADGETS:   parse_section      (section, save.quick_switch_gadgets);      break;
		}
	}
	for(const std::vector<Section>& sections : file.slot.levels) {
		LevelSaveGame& level_save_game = save.levels.emplace_back();
		for(const Section& section : sections) {
			Buffer buffer = section.data;
			switch(section.type) {
				case ST_LEVELSAVEDATA: level_save_game.level = buffer.read<LevelSave>(0); break;
			}
		}
	}
	return save;
}

void update_slot(File& dest, const SaveGame& save) {
	for(Section& section : dest.slot.sections) {
		OutBuffer buffer = section.data;
		switch(section.type) {
			case ST_LEVEL:                update_section      (buffer, save.level);                      break;
			case ST_ELAPSEDTIME:          update_section      (buffer, save.elapsed_time);               break;
			case ST_LASTSAVETIME:         update_section      (buffer, save.last_save_time);             break;
			case ST_GLOBALFLAGS:          update_section_array(buffer, save.global_flags);               break;
			case ST_CHEATSACTIVATED:      update_section_array(buffer, save.cheats_activated);           break;
			case ST_SKILLPOINTS:          update_section_array(buffer, save.skill_points);               break;
			case ST_HELPDATAMESSAGES:     update_section_array(buffer, save.help_data_messages);         break;
			case ST_HELPDATAMISC:         update_section_array(buffer, save.help_data_misc);             break;
			case ST_HELPDATAGADGETS:      update_section_array(buffer, save.help_data_gadgets);          break;
			case ST_CHEATSEVERACTIVATED:  update_section_array(buffer, save.cheats_ever_activated);      break;
			case ST_SETTINGS:             update_section      (buffer, save.settings);                   break;
			case ST_HEROSAVE:             update_section      (buffer, save.hero_save);                  break;
			case ST_MOVIESPLAYEDRECORD:   update_section_array(buffer, save.movies_played_record);       break;
			case ST_TOTALPLAYTIME:        update_section      (buffer, save.total_play_time);            break;
			case ST_TOTALDEATHS:          update_section      (buffer, save.total_deaths);               break;
			case ST_HELPLOG:              update_section_array(buffer, save.help_log);                   break;
			case ST_HELPLOGPOS:           update_section      (buffer, save.help_log_pos);               break;
			case ST_HEROGADGETBOX:        update_section      (buffer, save.hero_gadget_box);            break;
			case ST_PURCHASEABLEGADGETS:  update_section_array(buffer, save.purchaseable_gadgets);       break;
			case ST_PURCHASEABLEBOTUPGRD: update_section_array(buffer, save.purchaseable_bot_upgrades);  break;
			case ST_PURCHASEABLEWRENCH:   update_section      (buffer, save.purchaseable_wrench_level);  break;
			case ST_PURCHASEABLEPOSTMODS: update_section_array(buffer, save.purchaseable_post_fx_mods);  break;
			case ST_BOTSAVE:              update_section      (buffer, save.bot_save);                   break;
			case ST_FIRSTPERSONMODE:      update_section_array(buffer, save.first_person_desired_mode);  break;
			case ST_SAVEDDIFFICULTYLEVEL: update_section      (buffer, save.saved_difficulty_level);     break;
			case ST_PLAYERSTATISTICS:     update_section_array(buffer, save.player_statistics);          break;
			case ST_BATTLEDOMEWINSLOSSES: update_section_array(buffer, save.battledome_wins_and_losses); break;
			case ST_ENEMYKILLS:           update_section_array(buffer, save.enemy_kills);                break;
			case ST_QUICKSWITCHGADGETS:   update_section      (buffer, save.quick_switch_gadgets);       break;
		}
	}
	assert(dest.slot.levels.size() == save.levels.size());
	for(size_t i = 0; i < dest.slot.levels.size(); i++) {
		for(Section& section : dest.slot.levels[i]) {
			OutBuffer buffer = section.data;
			switch(section.type) {
				case ST_LEVELSAVEDATA: update_section(buffer, save.levels[i].level); break;
			}
		}
	}
}

const char* section_type(u32 type) {
	switch(type) {
		case ST_LEVEL: return "level ID";
		case ST_ELAPSEDTIME: return "elapsed time";
		case ST_LASTSAVETIME: return "last save time";
		case ST_GLOBALFLAGS: return "global flags";
		case ST_CHEATSACTIVATED: return "cheats activated";
		case ST_SKILLPOINTS: return "skill points";
		case ST_HELPDATAMESSAGES: return "help data messages";
		case ST_HELPDATAMISC: return "help data misc";
		case ST_HELPDATAGADGETS: return "help data gadgets";
		case ST_CHEATSEVERACTIVATED: return "cheats ever activated";
		case ST_SETTINGS: return "settings";
		case ST_HEROSAVE: return "hero save";
		case ST_MOVIESPLAYEDRECORD: return "movies played record";
		case ST_TOTALPLAYTIME: return "total play time";
		case ST_TOTALDEATHS: return "total deaths";
		case ST_HELPLOG: return "help log";
		case ST_HELPLOGPOS: return "help log pos";
		case ST_GAMEMODEOPTIONS: return "game mode options";
		case ST_MPPROFILES: return "mp profiles";
		case ST_HEROGADGETBOX: return "hero gadget box";
		case ST_LEVELSAVEDATA: return "level save data";
		case ST_PURCHASEABLEGADGETS: return "purchaseable gadgets";
		case ST_PURCHASEABLEBOTUPGRD: return "purchaseable bot upgrades";
		case ST_PURCHASEABLEWRENCH: return "purchaseable wrench";
		case ST_PURCHASEABLEPOSTMODS: return "purchaseable post fx mods";
		case ST_BOTSAVE: return "bot save";
		case ST_FIRSTPERSONMODE: return "first person mode";
		case ST_SAVEDDIFFICULTYLEVEL: return "saved difficulty level";
		case ST_PLAYERSTATISTICS: return "player statistics";
		case ST_BATTLEDOMEWINSLOSSES: return "battle dome wins and losses";
		case ST_ENEMYKILLS: return "enemy kills";
		case ST_QUICKSWITCHGADGETS: return "quick select gadgets";
	}
	return "???";
}

const std::vector<FileFormat> FILE_FORMATS = {
	{Game::DL, FileType::SLOT, {
		ST_LEVEL,
		ST_HEROSAVE,
		ST_ELAPSEDTIME,
		ST_LASTSAVETIME,
		ST_TOTALPLAYTIME,
		ST_SAVEDDIFFICULTYLEVEL,
		ST_GLOBALFLAGS,
		ST_CHEATSACTIVATED,
		ST_CHEATSEVERACTIVATED,
		ST_SKILLPOINTS,
		ST_HEROGADGETBOX,
		ST_HELPDATAMESSAGES,
		ST_HELPDATAMISC,
		ST_HELPDATAGADGETS,
		ST_SETTINGS,
		ST_FIRSTPERSONMODE,
		ST_MOVIESPLAYEDRECORD,
		ST_TOTALDEATHS,
		ST_HELPLOG,
		ST_HELPLOGPOS,
		ST_PURCHASEABLEGADGETS
	}, {
		ST_LEVELSAVEDATA
	}}
};


}
