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

#include "editor_gui.h"
#include "gui/commands.h"

#include <nfd.h>
#include <gui/gui.h>
#include <gui/config.h>
#include <gui/build_settings.h>
#include <gui/command_output.h>
#include <editor/app.h>
#include <editor/gui/view_3d.h>
#include <editor/gui/inspector.h>
#include <editor/gui/asset_selector.h>

struct Layout {
	const char* name;
	void (*menu_bar_extras)();
	void (*tool_bar)();
	std::vector<const char*> visible_windows;
	bool hovered = false;
};

static void menu_bar();
static void level_editor_menu_bar();
static void occlusion_things(Level* level);
static void tool_bar();
static void begin_dock_space();
static void dockable_windows();
static void dockable_window(const char* window, void (*func)());
static void end_dock_space();
static void create_dock_layout();

static bool layout_button(Layout& layout, size_t i);

static Layout layouts[] = {
	//{"Asset Browser", nullptr, nullptr, {}},
	{"Level Editor", level_editor_menu_bar, tool_bar, {}}
};
static size_t selected_layout = 0;
static ImRect available_rect;

void editor_gui() {
	available_rect = ImRect(ImVec2(0, 0), ImGui::GetMainViewport()->Size);
	
	menu_bar();
	if(layouts[selected_layout].tool_bar) {
		layouts[selected_layout].tool_bar();
	}
	
	begin_dock_space();
	dockable_windows();
	
	static bool is_first_frame = true;
	if(is_first_frame) {
		create_dock_layout();
		is_first_frame = false;
	}
	
	end_dock_space();
}

static void menu_bar() {
	if(ImGui::BeginMainMenuBar()) {
		bool open_error_popup = false;
		static std::string error_message;
		
		if(ImGui::BeginMenu("File")) {
			if(ImGui::MenuItem("Save")) {
				if(BaseEditor* editor = g_app->get_editor()) {
					try {
						editor->save(fs::path());
					} catch(SaveError& e) {
						if(e.retry) {
							nfdchar_t* path;
							nfdresult_t result = NFD_SaveDialog("asset", nullptr, &path);
							if(result == NFD_OKAY) {
								try {
									editor->save(path);
								} catch(SaveError& e) {
									error_message = e.message;
									open_error_popup = true;
								}
								free(path);
							}
						} else {
							error_message = e.message;
							open_error_popup = true;
						}
					}
				} else {
					error_message = "No editor open.";
					open_error_popup = true;
				}
			}
			ImGui::EndMenu();
		}
		
		if(ImGui::BeginMenu("Edit")) {
			if(ImGui::MenuItem("Undo")) {
				if(BaseEditor* editor = g_app->get_editor()) {
					try {
						editor->undo();
					} catch(RuntimeError& e) {
						error_message = e.message;
						open_error_popup = true;
					}
				} else {
					error_message = "No editor open.";
					open_error_popup = true;
				}
			}
			if(ImGui::MenuItem("Redo")) {
				if(BaseEditor* editor = g_app->get_editor()) {
					try {
						editor->redo();
					} catch(RuntimeError& e) {
						error_message = e.message;
						open_error_popup = true;
					}
				} else {
					error_message = "No editor open.";
					open_error_popup = true;
				}
			}
			ImGui::EndMenu();
		}
		
		if(open_error_popup) {
			ImGui::OpenPopup("Error");
		}
		
		ImGui::SetNextWindowSize(ImVec2(300, 200));
		if(ImGui::BeginPopupModal("Error")) {
			ImGui::TextWrapped("%s", error_message.c_str());
			if(ImGui::Button("Okay")) {
				error_message.clear();
				ImGui::CloseCurrentPopup();
			}
			ImGui::EndPopup();
		}
		
		if(ImGui::BeginMenu("View")) {
			app& a = *g_app;
			if(ImGui::MenuItem("Reset Camera")) {
				reset_camera(&a);
			}
			if(ImGui::BeginMenu("Visibility")) {
				ImGui::Checkbox("Tfrags", &a.render_settings.draw_tfrags);
				ImGui::Checkbox("Moby Instances", &a.render_settings.draw_moby_instances);
				ImGui::Checkbox("Moby Groups", &a.render_settings.draw_moby_groups);
				ImGui::Checkbox("Tie Instances", &a.render_settings.draw_tie_instances);
				ImGui::Checkbox("Tie Groups", &a.render_settings.draw_tie_groups);
				ImGui::Checkbox("Shrub Instances", &a.render_settings.draw_shrub_instances);
				ImGui::Checkbox("Shrub Groups", &a.render_settings.draw_shrub_groups);
				ImGui::Checkbox("Point Lights", &a.render_settings.draw_point_lights);
				ImGui::Checkbox("Env Sample Points", &a.render_settings.draw_env_sample_points);
				ImGui::Checkbox("Env Transitions", &a.render_settings.draw_env_transitions);
				ImGui::Checkbox("Cuboids", &a.render_settings.draw_cuboids);
				ImGui::Checkbox("Spheres", &a.render_settings.draw_spheres);
				ImGui::Checkbox("Cylinders", &a.render_settings.draw_cylinders);
				ImGui::Checkbox("Pills", &a.render_settings.draw_pills);
				ImGui::Checkbox("Cameras", &a.render_settings.draw_cameras);
				ImGui::Checkbox("Sound Instances", &a.render_settings.draw_sound_instances);
				ImGui::Checkbox("Paths", &a.render_settings.draw_paths);
				ImGui::Checkbox("Grind Paths", &a.render_settings.draw_grind_paths);
				ImGui::Checkbox("Areas", &a.render_settings.draw_areas);
				ImGui::Checkbox("Collision", &a.render_settings.draw_collision);
				ImGui::Separator();
				ImGui::Checkbox("Selected Instance Normals", &a.render_settings.draw_selected_instance_normals);
				ImGui::EndMenu();
			}
			ImGui::EndMenu();
		}
		
		for(size_t i = 0; i < ARRAY_SIZE(layouts); i++) {
			Layout& layout = layouts[i];
			if(layout_button(layout, i)) {
				selected_layout = i;
			}
		}
		
		static gui::PackerParams params;
		std::vector<std::string>& game_builds = g_app->game_bank->game_info.builds;
		std::vector<std::string>& mod_builds = g_app->mod_bank->game_info.builds;
		
		ImGui::SetNextItemWidth(200);
		gui::build_settings(params, &game_builds, mod_builds, false); 
		
		static CommandThread pack_command;
		static std::string iso_path;
		if(ImGui::Button("Build & Run##the_button")) {
			params.game_path = g_app->game_path;
			params.mod_paths = {g_app->mod_path};
			iso_path = gui::run_packer(params, pack_command);
			ImGui::OpenPopup("Build & Run##the_popup");
		}
		
		gui::command_output_screen("Build & Run##the_popup", pack_command, []() {}, []() {
			gui::EmulatorParams params;
			params.iso_path = iso_path;
			gui::run_emulator(params, false);
		});
		
		if(layouts[selected_layout].menu_bar_extras) {
			layouts[selected_layout].menu_bar_extras();
		}
		
		auto& pos = g_app->render_settings.camera_position;
		auto& rot = g_app->render_settings.camera_rotation;
		ImGui::Text("Cam (toggle with Z): X=%.2f Y=%.2f Z=%.2f Pitch=%.2f Yaw=%.2f",
			pos.x, pos.y, pos.z, rot.x, rot.y);
		
		ImGui::Text("Octant: %d %d %d", (s32) (pos.x / 4), (s32) (pos.y / 4), (s32) (pos.z / 4));
		
		available_rect.Min.y += ImGui::GetWindowSize().y - 1;
		ImGui::Text("Frame Time: %.2fms", g_app->delta_time * 1000.f);
		
		ImGui::EndMainMenuBar();
	}
}

static void level_editor_menu_bar() {
	const char* preview_value = "(select level)";
	Level* level = g_app->get_level();
	if(level) {
		preview_value = "(level)";
	}
	
	static AssetSelector level_selector;
	level_selector.required_type = LevelAsset::ASSET_TYPE;
	ImGui::SetNextItemWidth(200);
	if(Asset* asset = asset_selector("##level_selector", preview_value, level_selector, g_app->asset_forest)) {
		g_app->load_level(asset->as<LevelAsset>());
	}
	
	occlusion_things(level);
}

static void occlusion_things(Level* level) {
	static CommandThread occl_command;
	static gui::RebuildOcclusionParams occl_params;
	
	const char* occlusion_problem = nullptr;
	if(!occlusion_problem && !level) {
		occlusion_problem = "No level loaded.";
	}
	if(!occlusion_problem && !level->level_wad().has_occlusion())  {
		occlusion_problem = "Missing occlusion asset.";
	}
	if(!occlusion_problem && !level->level_wad().get_occlusion().has_octants())  {
		occlusion_problem = "Occlusion asset has missing octants attribute.";
	}
	if(!occlusion_problem && !level->level_wad().get_occlusion().has_grid())  {
		occlusion_problem = "Occlusion asset has missing grid attribute.";
	}
	if(!occlusion_problem && !level->level_wad().get_occlusion().has_mappings())  {
		occlusion_problem = "Occlusion asset has mappings attribute.";
	}
	
	if(ImGui::Button("Rebuild Occlusion##the_button") && !occlusion_problem) {
		// Setup the file structure so that the new occlusion file can be
		// written out in place of the old one.
		if(&level->level_wad().get_occlusion().bank() != g_app->mod_bank && level->level().parent()) {
			s32 level_id = level->level_wad().id();
			std::string path = generate_asset_path<LevelAsset>("levels", "level", level_id, *level->level().parent());
			OcclusionAsset& old_occl = level->level_wad().get_occlusion();
			
			// Create a new .asset file for the occlusion data.
			AssetFile& occlusion_file = g_app->mod_bank->asset_file(path);
			AssetLink link = level->level_wad().get_occlusion().absolute_link();
			OcclusionAsset& new_occl = occlusion_file.asset_from_link(OcclusionAsset::ASSET_TYPE, link).as<OcclusionAsset>();
			
			// Copy the old grid and mappings to the mod asset bank as placeholders.
			std::string octants = old_occl.file().read_text_file(old_occl.octants().path);
			std::unique_ptr<InputStream> grid_src = old_occl.file().open_binary_file_for_reading(old_occl.grid());
			std::unique_ptr<InputStream> mappings_src = old_occl.file().open_binary_file_for_reading(old_occl.mappings());
			FileReference octants_ref = new_occl.file().write_text_file("occlusion_octants.csv", octants.c_str());
			auto [grid_dest, grid_ref] = new_occl.file().open_binary_file_for_writing("occlusion_grid.bin");
			auto [mappings_dest, mappings_ref] = new_occl.file().open_binary_file_for_writing("occlusion_mappings.bin");
			Stream::copy(*grid_dest, *grid_src, grid_src->size());
			Stream::copy(*mappings_dest, *mappings_src, grid_src->size());
			new_occl.set_octants(octants_ref);
			new_occl.set_grid(grid_ref);
			new_occl.set_mappings(mappings_ref);
			
			// Write the .asset file.
			occlusion_file.write();
		}
		
		occl_params.game_path = g_app->game_path;
		occl_params.mod_path = g_app->mod_path;
		occl_params.level_wad_asset = level->level_wad().absolute_link().to_string();
		gui::run_occlusion_rebuild(occl_params, occl_command);
		
		ImGui::OpenPopup("Rebuild Occlusion##the_popup");
	} else if(ImGui::IsItemHovered() && occlusion_problem) {
		ImGui::BeginTooltip();
		ImGui::Text("%s", occlusion_problem);
		ImGui::EndTooltip();
	}
	
	gui::command_output_screen("Rebuild Occlusion##the_popup", occl_command, []() {});
}

static void tool_bar() {
	ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.f);
	ImGuiViewport* view = ImGui::GetMainViewport();
	ImGui::SetNextWindowPos(available_rect.Min - ImVec2(1, 0));
	
	ImGui::SetNextWindowSize(ImVec2(56 * g_config.ui.scale, view->Size.y));
	ImGui::Begin("Tools", nullptr,
		ImGuiWindowFlags_NoDecoration |
		ImGuiWindowFlags_NoMove);
	ImGui::PopStyleVar();
	
	for(std::size_t i = 0; i < g_app->tools.size(); i++) {
		bool active = i == g_app->active_tool_index;
		if(!active) {
			ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
		}

		bool clicked = ImGui::ImageButton(
			(void*) (intptr_t) g_app->tools[i]->icon(), ImVec2(32*g_config.ui.scale, 32*g_config.ui.scale),
						ImVec2(0, 0), ImVec2(1, 1), -1);
		if(!active) {
			ImGui::PopStyleColor();
		}
		if(clicked) {
			g_app->active_tool_index = i;
		}
	}
	
	available_rect.Min.x += ImGui::GetWindowSize().x;
	
	ImGui::End();
}

static void begin_dock_space() {
	ImGuiWindowFlags window_flags =  ImGuiWindowFlags_NoDocking;
	ImGui::SetNextWindowPos(available_rect.Min);
	ImGui::SetNextWindowSize(available_rect.Max - available_rect.Min);
	ImGui::SetNextWindowViewport(ImGui::GetWindowViewport()->ID);
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
	ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), ImGuiDockNodeFlags_None);
}

static void dockable_windows() {
	ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
	dockable_window("3D View", view_3d);
	ImGui::PopStyleVar();
	dockable_window("Inspector", inspector);
}

static void dockable_window(const char* window, void (*func)()) {
	ImGui::Begin(window);
	func();
	ImGui::End();
}

static void end_dock_space() {
	ImGui::End();
}

static void create_dock_layout() {
	ImGuiID dockspace_id = ImGui::GetID("dock_space");
	
	ImGui::DockBuilderRemoveNode(dockspace_id);
	ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_DockSpace);
	ImGui::DockBuilderSetNodeSize(dockspace_id, ImVec2(1.f, 1.f));

	ImGuiID left_centre, right;
	ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Left, 8.f / 10.f, &left_centre, &right);
	
	ImGui::DockBuilderDockWindow("3D View", left_centre);
	ImGui::DockBuilderDockWindow("Inspector", right);

	ImGui::DockBuilderFinish(dockspace_id);
}

static bool layout_button(Layout& layout, size_t i) {
	bool selected = i == selected_layout;
	ImGuiID id = ImGui::GetID(layout.name);
	ImGuiCol colour =
		selected       ? ImGuiCol_TabActive  :
		layout.hovered ? ImGuiCol_TabHovered :
		ImGuiCol_Tab;
	ImDrawList* dl = ImGui::GetWindowDrawList();
	ImVec2 pos = ImGui::GetCursorPos();
	ImVec2 size = ImGui::TabItemCalcSize(layout.name, false);
	ImRect bb(pos, pos + size);
	ImGui::ItemAdd(bb, id);
	ImGui::TabItemBackground(dl, bb, ImGuiTabItemFlags_None, ImGui::GetColorU32(colour));
	ImGui::TabItemLabelAndCloseButton(dl, bb, ImGuiTabItemFlags_None, ImGui::GetStyle().FramePadding, layout.name, ImGui::GetID(layout.name), 0, true, nullptr, nullptr);
	bool pressed = ImGui::ButtonBehavior(bb, id, &layout.hovered, nullptr, ImGuiButtonFlags_PressedOnClickRelease);
	ImGui::SetCursorPos(pos + ImVec2(size.x + 4*g_config.ui.scale, 0));
	return pressed;
}
