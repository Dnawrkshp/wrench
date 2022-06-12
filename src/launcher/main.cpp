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

#include <nfd.h>
#include <core/png.h>
#include <core/shell.h>
#include <engine/compression.h>
#include <gui/gui.h>
#include <gui/about.h>
#include <gui/settings.h>
#include <gui/command_output.h>
#include <launcher/oobe.h>
#include <launcher/mod_list.h>
#include <launcher/game_list.h>
#include <launcher/global_state.h>
#include <launcher/new_mod_screen.h>

static void update_gui(f32 delta_time);
static Mod* mod_list_window();
static void details_window(Mod* mod);
static void not_specified();
static void buttons_window(f32 buttons_window_height);
static void create_dock_layout();
static void begin_docking(f32 buttons_window_height);
static ImFont* load_font_from_gui_wad(SectorRange range);
static Texture load_image_from_launcher_wad(SectorRange range);

int main(int argc, char** argv) {
	g_launcher.mode = LauncherMode::DRAWING_GUI;
	verify(g_launcher.wad.open("data/launcher.wad"), "Failed to open 'launcher.wad'.");
	verify(g_launcher.buildwad.open("data/build.wad"), "Failed to open 'build.wad'.");
	
	g_launcher.bin_paths.wrenchbuild = "./bin/wrenchbuild";
	
	if(!gui::config_file_exists()) {
		if(!run_oobe()) {
			return 0;
		}
	} else {
		g_config.read();
	}
	
	for(;;) {
		switch(g_launcher.mode) {
			case LauncherMode::DRAWING_GUI: {
				g_launcher.window = gui::startup("Wrench Launcher", 960, 600);
				
				g_launcher.font_regular = load_font_from_gui_wad(wadinfo.gui.fonts[0]);
				g_launcher.font_italic = load_font_from_gui_wad(wadinfo.gui.fonts[1]);
				g_launcher.placeholder_image = load_image_from_launcher_wad(wadinfo.launcher.placeholder_images[0]);
				
				load_game_list(g_config.folders.games_folder);
				load_mod_list(g_config.folders.mods_folders);
				
				while(g_launcher.mode == LauncherMode::DRAWING_GUI) {
					gui::run_frame(g_launcher.window, update_gui);
					
					if(glfwWindowShouldClose(g_launcher.window)) {
						g_launcher.mode = LauncherMode::EXIT;
					}
				}
				
				g_launcher.font_buffers.clear();
				g_launcher.placeholder_image.destroy();
				
				free_game_list();
				free_mod_list();
				
				gui::shutdown(g_launcher.window);
				break;
			}
			case LauncherMode::RUNNING_EMULATOR: {
				const char* args[] = {
					g_launcher.launch_params.emulator_executable.c_str(),
					g_launcher.launch_params.iso_path.c_str()
				};
				execute_command(ARRAY_SIZE(args), args, nullptr);
				g_launcher.mode = LauncherMode::DRAWING_GUI;
				break;
			}
			case LauncherMode::EXIT: {
				return 0;
			}
		}
	}
}

void update_gui(f32 delta_time) {
	f32 button_height = GImGui->Font->FontSize + GImGui->Style.FramePadding.y * 2.f;
	f32 buttons_window_height = button_height + GImGui->Style.WindowPadding.y * 2.f;
	
	begin_docking(buttons_window_height);
	
	static bool first_frame = true;
	if(first_frame) {
		create_dock_layout();
		first_frame = false;
	}
	
	Mod* mod = mod_list_window();
	details_window(mod);
	
	ImGui::End(); // docking
	
	buttons_window(buttons_window_height);
}

static Mod* mod_list_window() {
	ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
	ImGui::Begin("Mod List");
	
	char greeting[64];
	BuildWadHeader& build = wadinfo.build;
	if(build.version_major > -1 && build.version_minor > -1) {
		snprintf(greeting, sizeof(greeting), "Wrench Modding Toolset v%hd.%hd", build.version_major, build.version_minor);
	} else {
		snprintf(greeting, sizeof(greeting), "Wrench Modding Toolset");
	}
	
	ImGui::NewLine();
	f32 greeting_width = ImGui::CalcTextSize(greeting).x;
	ImGui::SetCursorPosX((ImGui::GetWindowSize().x - greeting_width) / 2);
	ImGui::Text("%s", greeting);
	ImGui::NewLine();
	
	static std::string filter;
	
	if(ImGui::BeginTable("inputs", 2)) {
		ImGui::TableSetupColumn("labels", ImGuiTableColumnFlags_WidthFixed);
		ImGui::TableSetupColumn("inputs", ImGuiTableColumnFlags_WidthStretch);
		
		ImGui::TableNextRow();
		ImGui::TableNextColumn();
		ImGui::AlignTextToFramePadding();
		ImGui::Text(" Game");
		ImGui::TableNextColumn();
		ImGui::SetNextItemWidth(-1);
		game_list();
		
		ImGui::TableNextRow();
		ImGui::TableNextColumn();
		ImGui::AlignTextToFramePadding();
		ImGui::Text(" Filter");
		ImGui::TableNextColumn();
		ImGui::SetNextItemWidth(-1);
		ImGui::InputText("##filter", &filter);
		
		ImGui::EndTable();
	}
	
	ImGui::BeginChild("table");
	Mod* mod = mod_list(filter);
	ImGui::EndChild();
	
	ImGui::End();
	ImGui::PopStyleVar(); // ImGuiStyleVar_WindowPadding
	
	return mod;
}

static void details_window(Mod* mod) {
	ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
	ImGui::Begin("Details");
	
	assert(g_mod_images.size() >= 1);
	ModImage& image = g_mod_images[0];
	
	f32 display_height = 320.f;
	f32 display_width = display_height * image.width / (f32) image.height;
	
	ImGui::SetCursorPosX(ImGui::GetWindowWidth() / 2 - display_width / 2);
	ImGui::Image((void*) (intptr_t) image.texture.id, ImVec2(display_width, display_height));
	
	if(mod) {
		if(ImGui::BeginTable("attributes", 2)) {
			ImGui::TableSetupColumn("key", ImGuiTableColumnFlags_WidthFixed);
			ImGui::TableSetupColumn("value", ImGuiTableColumnFlags_WidthStretch);
			
			ImGui::TableNextRow();
			ImGui::TableNextColumn();
			ImGui::Text("Author");
			ImGui::TableNextColumn();
			if(!mod->info.author.empty()) {
				ImGui::TextWrapped("%s", mod->info.author.c_str());
			} else {
				not_specified();
			}
			
			ImGui::TableNextRow();
			ImGui::TableNextColumn();
			ImGui::Text("Description");
			ImGui::TableNextColumn();
			if(!mod->info.description.empty()) {
				ImGui::TextWrapped("%s", mod->info.description.c_str());
			} else {
				not_specified();
			}
			
			ImGui::TableNextRow();
			ImGui::TableNextColumn();
			ImGui::Text("Version");
			ImGui::TableNextColumn();
			if(!mod->info.version.empty()) {
				ImGui::TextWrapped("%s", mod->info.version.c_str());
			} else {
				not_specified();
			}
			
			ImGui::EndTable();
		}
	}
	
	ImGui::End();
	ImGui::PopStyleVar(); // ImGuiStyleVar_WindowPadding
}

static void not_specified() {
	ImGui::PushFont(g_launcher.font_italic);
	ImGui::TextWrapped("Not specified.");
	ImGui::PopFont();
}

static void buttons_window(f32 buttons_window_height) {
	ImGuiWindowFlags flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize;
	
	ImVec2 viewport_size = ImGui::GetMainViewport()->Size;
	ImGui::SetNextWindowPos(ImVec2(-1, viewport_size.y - buttons_window_height));
	ImGui::SetNextWindowSize(ImVec2(viewport_size.x + 2, buttons_window_height + 1));
	
	bool p_open;
	ImGui::Begin("Buttons", &p_open, flags);
	
	if(ImGui::Button("Import ISO")) {
		nfdchar_t* path = nullptr;
		nfdresult_t result = NFD_OpenDialog("iso", nullptr, &path);
		if(result == NFD_OKAY) {
			std::vector<std::string> args = {
				g_launcher.bin_paths.wrenchbuild,
				"unpack",
				path,
				"-o",
				g_config.folders.games_folder.c_str(),
				"-s" // Unpack it into a subdirectory.
			};
			spawn_command_thread(args, &g_launcher.import_iso_command);
			free(path);
			
			ImGui::OpenPopup("Import ISO");
		} else if(result != NFD_CANCEL) {
			printf("error: %s\n", NFD_GetError());
		}
	}
	
	gui::command_output_screen("Import ISO", g_launcher.import_iso_command, []() {
		load_game_list(g_config.folders.games_folder);
	});
	
	ImGui::SameLine();
	if(ImGui::Button("Open Mods Folder")) {
		if(g_config.folders.mods_folders.size() == 1) {
			open_in_file_manager(fs::absolute(g_config.folders.mods_folders[0]).string().data());
		} else {
			ImGui::OpenPopup("Mods Folder Selector");
		}
	}
	
	if(ImGui::BeginPopup("Mods Folder Selector")) {
		for(const std::string& mods_folder : g_config.folders.mods_folders) {
			if(ImGui::Selectable(mods_folder.c_str())) {
				open_in_file_manager(fs::absolute(mods_folder).string().data());
			}
		}
		ImGui::EndPopup();
	}
	
	ImGui::SameLine();
	if(ImGui::Button("Refresh")) {
		load_game_list(g_config.folders.games_folder);
		load_mod_list(g_config.folders.mods_folders);
	}
	
	if(new_mod_screen()) {
		load_mod_list(g_config.folders.mods_folders);
	}
	
	ImGui::SameLine();
	if(ImGui::Button("Open in Editor")) {
		
	}
	
	ImGui::SameLine();
	if(ImGui::Button("···")) {
		ImGui::OpenPopup("More Buttons");
	}
	
	static bool open_new_mod = false;
	static bool open_about = false;
	static bool open_settings = false;
	static bool show_the_demo = false;
	
	if(ImGui::BeginPopup("More Buttons")) {
		if(ImGui::Selectable("New Mod##the_button")) {
			open_new_mod = true;
		}
		ImGui::Separator();
		if(ImGui::Selectable("About##the_button")) {
			open_about = true;
		}
		if(ImGui::Selectable("Settings##the_button")) {
			open_settings = true;
		}
		if(g_config.ui.developer) {
			if(ImGui::BeginMenu("Developer")) {
				if(ImGui::Selectable("The Demo")) {
					show_the_demo = !show_the_demo;
				}
				ImGui::EndMenu();
			}
		}
		ImGui::EndPopup();
	}
	
	if(open_new_mod) {
		ImGui::OpenPopup("New Mod##the_popup");
		open_new_mod = false;
	}
	
	if(open_about) {
		ImGui::OpenPopup("About##the_popup");
		open_about = false;
	}
	
	gui::about_screen();
	
	if(open_settings) {
		ImGui::OpenPopup("Settings##the_popup");
		open_settings = false;
	}
	
	update_settings_popup();
	
	if(show_the_demo) {
		ImGui::ShowDemoWindow();
	}
	
	ImGui::SameLine();
	f32 build_run_text_width = ImGui::CalcTextSize("Build & Run").x;
	ImGuiStyle& s = ImGui::GetStyle();
	f32 build_run_button_width = s.FramePadding.x + build_run_text_width + s.FramePadding.x;
	f32 build_area_width = 192 + s.ItemSpacing.x + build_run_button_width + s.WindowPadding.x;
	ImGui::SetCursorPosX(viewport_size.x - build_area_width);
	
	ImGui::PushItemWidth(192);
	if(ImGui::BeginCombo("##builds", "dl")) {
		ImGui::EndCombo();
	}
	ImGui::PopItemWidth();
	
	ImGui::SameLine();
	if(ImGui::Button("Build & Run")) {
		
	}
	
	ImGui::End();
}

static void create_dock_layout() {
	ImGuiID dockspace_id = ImGui::GetID("dock_space");
	
	ImGui::DockBuilderRemoveNode(dockspace_id);
	ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_DockSpace);
	ImGui::DockBuilderSetNodeSize(dockspace_id, ImVec2(1.f, 1.f));

	ImGuiID mod_list, details;
	ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Left, 0.5f, &mod_list, &details);
	
	ImGui::DockBuilderDockWindow("Mod List", mod_list);
	ImGui::DockBuilderDockWindow("Details", details);

	ImGui::DockBuilderFinish(dockspace_id);
}

static void begin_docking(f32 buttons_window_height) {
	ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoDocking;
	
	ImGuiViewport* viewport = ImGui::GetMainViewport();
	ImVec2 size = viewport->Size;
	size.y -= buttons_window_height;
	ImGui::SetNextWindowPos(viewport->Pos);
	ImGui::SetNextWindowSize(size);
	ImGui::SetNextWindowViewport(viewport->ID);
	ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
	ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
	window_flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
	window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus;
	
	static bool p_open;
	ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
	ImGui::Begin("dock_space", &p_open, window_flags);
	ImGui::PopStyleVar();
	
	ImGui::PopStyleVar(2);

	ImGuiID dockspace_id = ImGui::GetID("dock_space");
	ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), ImGuiDockNodeFlags_NoTabBar);
}

static ImFont* load_font_from_gui_wad(SectorRange range) {
	std::vector<u8> compressed_font = g_guiwad.read_multiple<u8>(range.offset.bytes(), range.size.bytes());
	std::vector<u8>& decompressed_font = g_launcher.font_buffers.emplace_back();
	decompress_wad(decompressed_font, compressed_font);
	
	ImFontConfig font_cfg;
	font_cfg.FontDataOwnedByAtlas = false;
	
	ImGuiIO& io = ImGui::GetIO();
	return io.Fonts->AddFontFromMemoryTTF(decompressed_font.data(), decompressed_font.size(), 22, &font_cfg);
}

static Texture load_image_from_launcher_wad(SectorRange range) {
	std::vector<u8> compressed_image = g_launcher.wad.read_multiple<u8>(range.offset.bytes(), range.size.bytes());
	
	std::vector<u8> image;
	decompress_wad(image, compressed_image);
	
	s32 width = *(s32*) &image[0];
	s32 height = *(s32*) &image[4];
	std::vector<u8> data(&image[16], &image[16 + width * height * 4]);
	return Texture::create_rgba(width, height, std::move(data));
}
