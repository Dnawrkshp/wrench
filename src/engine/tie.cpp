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

#include "tie.h"

#include <algorithm>

static GcUyaDlTieClassHeader read_tie_header(Buffer src, Game game);
static TiePacket read_tie_packet(Buffer src, const TieClass& tie, int lod, const TiePacketHeader& header);
static void write_tie_packet(OutBuffer dest, const TiePacket& packet);
static ColladaScene recover_tie(const TieClass& tie);

TieClass read_tie_class(Buffer src, Game game) {
	TieClass tie;
	
	GcUyaDlTieClassHeader header = read_tie_header(src, game);
	tie.scale = header.scale;

	tie.ad_gifs = src.read_multiple<TieAdGifs>(header.ad_gif_ofs, header.texture_count, "ad gifs").copy();
	if(game == Game::DL) {
		tie.normals =
			src.read_multiple<TieNormal>(header.vert_normals + 0x10, header.vert_normal_count, "normals").copy();
	} else {
		tie.normals = src.read_multiple<TieNormal>(header.vert_normals, header.vert_normal_count, "normals").copy();
	}

	for (s32 i = 0; i < 3; i++) {
		Buffer rgba_remap_buffer = src.subbuf(header.vert_normals + header.rgba_remap_ofs[i]);
		auto tie_rgba_remaps = rgba_remap_buffer.read_multiple<u8>(0, 0x20, "rgba remap");
		auto rgba_remap_ofs = 0x20;
		for(s32 j = 0; j < tie_rgba_remaps.size(); ++j) {
			int remap_size = tie_rgba_remaps[j] * 0x10;
			if(remap_size == 0)
				break;

			auto remap_header = rgba_remap_buffer.read<TieRgbaRemapSegmentHeader>(rgba_remap_ofs, "rgba remap"); rgba_remap_ofs += 0x10;
			auto rgba_remap0 = rgba_remap_buffer.read_multiple<TieRgbaRemap>(rgba_remap_ofs, remap_header.remap_size[0] / sizeof(TieRgbaRemap), "rgba remap").copy(); rgba_remap_ofs += remap_header.remap_size[0];
			auto fat_rgba_remap1 = rgba_remap_buffer.read_multiple<TieFatRgbaRemap>(rgba_remap_ofs, remap_header.remap_size[1] / sizeof(TieFatRgbaRemap), "rgba remap").copy(); rgba_remap_ofs += remap_header.remap_size[1];
			auto fat_rgba_remap2 = rgba_remap_buffer.read_multiple<TieFatRgbaRemap>(rgba_remap_ofs, remap_header.remap_size[2] / sizeof(TieFatRgbaRemap), "rgba remap").copy(); rgba_remap_ofs += remap_header.remap_size[2];
			auto fat_rgba_remap3 = rgba_remap_buffer.read_multiple<TieFatRgbaRemap>(rgba_remap_ofs, remap_header.remap_size[3] / sizeof(TieFatRgbaRemap), "rgba remap").copy(); rgba_remap_ofs += remap_header.remap_size[3];
			auto fat_rgba_remap4 = rgba_remap_buffer.read_multiple<TieFatRgbaRemap>(rgba_remap_ofs, remap_header.remap_size[4] / sizeof(TieFatRgbaRemap), "rgba remap").copy(); rgba_remap_ofs += remap_header.remap_size[4];
			//auto rgba_remap5 = rgba_remap_buffer.read_multiple<TieRgbaRemap>(rgba_remap_ofs, remap_header.remap_size[5] / sizeof(TieRgbaRemap), "rgba remap").copy(); rgba_remap_ofs += remap_header.remap_size[5];
			rgba_remap_ofs += remap_header.remap_size[5];

			tie.rgba_remap[i].insert(tie.rgba_remap[i].end(), rgba_remap0.begin(), rgba_remap0.end());
			tie.fat_rgba_remap[i].insert(tie.fat_rgba_remap[i].end(), fat_rgba_remap1.begin(), fat_rgba_remap1.end());
			tie.fat_rgba_remap[i].insert(tie.fat_rgba_remap[i].end(), fat_rgba_remap2.begin(), fat_rgba_remap2.end());
			tie.fat_rgba_remap[i].insert(tie.fat_rgba_remap[i].end(), fat_rgba_remap3.begin(), fat_rgba_remap3.end());
			tie.fat_rgba_remap[i].insert(tie.fat_rgba_remap[i].end(), fat_rgba_remap4.begin(), fat_rgba_remap4.end());
			//tie.rgba_remap[i].insert(tie.rgba_remap[i].end(), rgba_remap5.begin(), rgba_remap5.end());

			rgba_remap_ofs = align32(rgba_remap_ofs, 0x10);
		}
	}
	
	for(s32 i = 0; i < 3; i++) {
		TieLod& lod = tie.lods[i];
		lod.packets.reserve(header.packet_count[i]);
		
		Buffer lod_buffer = src.subbuf(header.packets[i]);
		auto packet_table = lod_buffer.read_multiple<TiePacketHeader>(0, header.packet_count[i], "packet header");
		for(s32 j = 0; j < header.packet_count[i]; j++) {
			lod.packets.emplace_back(read_tie_packet(src.subbuf(header.packets[i] + packet_table[j].data), tie, i, packet_table[j]));
		}
	}
	
	return tie;
}

void write_tie_class(OutBuffer dest, const TieClass& tie) {
	
}

static GcUyaDlTieClassHeader read_tie_header(Buffer src, Game game) {
	GcUyaDlTieClassHeader header = {};
	
	if(game == Game::RAC) {
		RacTieClassHeader rac_header = src.read<RacTieClassHeader>(0, "header");
		memcpy(header.packets, rac_header.packets, sizeof(header.packets));
		memcpy(header.packet_count, rac_header.packet_count, sizeof(header.packet_count));
		header.texture_count = rac_header.texture_count;
		header.near_dist = rac_header.near_dist;
		header.mid_dist = rac_header.mid_dist;
		header.far_dist = rac_header.far_dist;
		header.ad_gif_ofs = rac_header.ad_gif_ofs;
		header.scale = rac_header.scale;
		header.bsphere = rac_header.bsphere;
	} else {
		header = src.read<GcUyaDlTieClassHeader>(0, "header");
	}
	
	return header;
}

static void populate_tie_normal(const TieClass& tie, int lod, int vertex_idx, TieVertex& vertex, bool is_fat)
{
	size_t j;

	for(s32 i = 0; i < 3; i++) {
		const std::vector<TieRgbaRemap>& rgba_remap = tie.rgba_remap[i];
		const std::vector<TieFatRgbaRemap>& fat_rgba_remap = tie.fat_rgba_remap[i];

		//
		if(!is_fat) {
			for(j = 0; j < rgba_remap.size(); j++) {
				if(rgba_remap[j].src_ofs == (vertex_idx * 4)) {
					auto normal_idx = (rgba_remap[j].dest_ofs & 0x7fc) / 4;
					if(normal_idx >= tie.normals.size())
						continue;

					const TieNormal& normal = tie.normals[normal_idx];
					vertex.azimuth = normal.azimuth;
					vertex.elevation = normal.elevation;
					return;
				}
			}

			for(j = 0; j < fat_rgba_remap.size(); j++) {
				if((fat_rgba_remap[j].remap & 0x7fc) == (vertex_idx * 4)) {
					auto normal_idx = ((fat_rgba_remap[j].remap >> 12) & 0x7fc) / 4;
					if(normal_idx >= tie.normals.size())
						continue;
					const TieNormal& normal = tie.normals[normal_idx];
					vertex.azimuth = normal.azimuth;
					vertex.elevation = normal.elevation;
					return;
				}
			}
		} else {
			for(j = 0; j < fat_rgba_remap.size(); j++) {
				if((fat_rgba_remap[j].remap & 0x7fc) == (vertex_idx * 4)) {
					auto normal_idx = ((fat_rgba_remap[j].remap >> 12) & 0x7fc) / 4;
					if(normal_idx >= tie.normals.size())
						continue;

					const TieNormal& normal = tie.normals[normal_idx];
					vertex.azimuth = normal.azimuth;
					vertex.elevation = normal.elevation;
					return;
				}
			}

			for(j = 0; j < rgba_remap.size(); j++) {
				if(rgba_remap[j].src_ofs == (vertex_idx * 4)) {
					auto normal_idx = (rgba_remap[j].dest_ofs & 0x7fc) / 4;
					if(normal_idx >= tie.normals.size())
						continue;

					const TieNormal& normal = tie.normals[normal_idx];
					vertex.azimuth = normal.azimuth;
					vertex.elevation = normal.elevation;
					return;
				}
			}
		}
	}

	int a = 0;
}

static TiePacket read_tie_packet(Buffer src, const TieClass& tie, int lod, const TiePacketHeader& header) {
	static size_t vertex_idx = 0;
	static int last_lod = 0;

	// reset rolling index with each lod
	if (last_lod != lod) {
		vertex_idx = 0;
		last_lod = lod;
	}

	TiePacket packet;

	auto ad_gif_dest_offsets = src.read_multiple<s32>(0x0, 4, "ad gif destination offsets");
	auto ad_gif_src_offsets = src.read_multiple<s32>(0x10, 4, "ad gif source offsets");
	TieUnpackHeader unpack_header = src.read<TieUnpackHeader>(0x20, "header");
	
	s32 strip_ofs = 0x2c;
	auto strips = src.read_multiple<TieStrip>(strip_ofs, unpack_header.strip_count, "strips");

	Buffer vertex_buffer = src.subbuf(header.vert_ofs * 0x10, header.vert_size * 0x10);
	s32 dinky_vertex_count = (unpack_header.dinky_vertices_size_plus_four - 4) / 2;
	std::vector<TieDinkyVertex> dinky_vertices = vertex_buffer.read_multiple<TieDinkyVertex>(0, dinky_vertex_count, "dinky vertices").copy();
	auto fat_vertices = vertex_buffer.read_all<TieFatVertex>(dinky_vertex_count * 0x10);
	
	// Combine the dinky and fat vertices.
	std::vector<TieVertex> raw_vertices;
	for(const TieDinkyVertex& src : dinky_vertices) {
		TieVertex& dest_1 = raw_vertices.emplace_back();

		dest_1.x = src.x;
		dest_1.y = src.y;
		dest_1.z = src.z;
		dest_1.s = src.s;
		dest_1.t = src.t;
		dest_1.q = src.q;
		dest_1.gs_packet_write_ofs = src.gs_packet_write_ofs;
		dest_1.gs_packet_write_ofs_2 = src.gs_packet_write_ofs_2;

		populate_tie_normal(tie, lod, vertex_idx, dest_1, false);

		if(src.gs_packet_write_ofs_2 != 0 && src.gs_packet_write_ofs_2 != src.gs_packet_write_ofs) {
			TieVertex& dest_2 = raw_vertices.emplace_back();

			dest_2 = dest_1;
			dest_2.gs_packet_write_ofs = src.gs_packet_write_ofs_2;
		}

		vertex_idx++;
	}
	for(const TieFatVertex& src : fat_vertices) {
		TieVertex& dest_1 = raw_vertices.emplace_back();

		dest_1.x = src.x;
		dest_1.y = src.y;
		dest_1.z = src.z;
		dest_1.gs_packet_write_ofs = src.gs_packet_write_ofs;
		dest_1.s = src.s;
		dest_1.t = src.t;
		dest_1.q = src.q;
		dest_1.gs_packet_write_ofs_2 = src.gs_packet_write_ofs_2;
		dest_1.azimuth = 0;
		dest_1.elevation = 0;

		populate_tie_normal(tie, lod, vertex_idx, dest_1, true);

		if(src.gs_packet_write_ofs_2 != 0 && src.gs_packet_write_ofs_2 != src.gs_packet_write_ofs) {
			TieVertex& dest_2 = raw_vertices.emplace_back();
			dest_2 = dest_1;
			dest_2.gs_packet_write_ofs = src.gs_packet_write_ofs_2;
		}

		vertex_idx++;
	}
	
	// The vertices in the file are not sorted by their GS packet address,
	// probably to help with buffering. For the purposes of reading ties, we
	// want to read them in the order they appear in the GS packet.
 	std::sort(BEGIN_END(raw_vertices), [](TieVertex& lhs, TieVertex& rhs)
		{ return lhs.gs_packet_write_ofs < rhs.gs_packet_write_ofs; });
	
	// Each packet must have a minimum of 4 regular vertices, so there may be
	// duplicates to pad out small packets. These can be safely removed.
	std::vector<TieVertex> vertices;
	if(raw_vertices.size() >= 1) {
		vertices.emplace_back(raw_vertices[0]);
	}
	for(size_t i = 1; i < raw_vertices.size(); i++) {
		if(raw_vertices[i].gs_packet_write_ofs != raw_vertices[i - 1].gs_packet_write_ofs) {
			vertices.emplace_back(raw_vertices[i]);
		}
	}

	s32 material_index = -1;
	TiePrimitive* prim = nullptr;
	
	s32 next_strip = 0;
	s32 next_vertex = 0;
	s32 next_ad_gif = 0;
	s32 next_offset = 0;
	
	// The first AD GIF is always at the beginning of the GIF packet.
	{
		material_index = ad_gif_src_offsets[next_ad_gif++] / 0x50;
		next_offset += 6;
	}
	
	// Interpret the data in the order it would appear in the GS packet.
	while(next_strip < strips.size() || next_vertex < vertices.size()) {
		// Data used to generate GIF tags for each of the strips.
		if(next_strip < strips.size() && strips[next_strip].gif_tag_offset == next_offset) {

			prim = &packet.primitives.emplace_back();
			prim->material_index = material_index;
			prim->face = strips[next_strip].pad_3 != 0;

			next_strip++;
			next_offset += 1;
			
			continue;
		}
		
		// Regular vertices.
		if(next_vertex < vertices.size() && vertices[next_vertex].gs_packet_write_ofs == next_offset) {
			verify(prim != nullptr, "Tie has bad GS packet data.");
			prim->vertices.emplace_back(vertices[next_vertex]);
			
			next_vertex++;
			next_offset += 3;
			
			continue;
		}
		
		// AD GIF data to change the texture.
		if(next_ad_gif < ad_gif_src_offsets.size() && ad_gif_dest_offsets[next_ad_gif - 1] == next_offset) {
			material_index = ad_gif_src_offsets[next_ad_gif] / 0x50;
			
			next_ad_gif++;
			next_offset += 6;
			
			continue;
		}
		
		verify_not_reached("Bad GS packet, expected next offset 0x%x.", next_offset);
	}
	
	return packet;
}

static void write_tie_packet(OutBuffer dest, const TiePacket& packet) {
	
}

ColladaScene recover_tie_class(const TieClass& tie) {
	ColladaScene scene;
	
	for(s32 i = 0; i < (s32) tie.ad_gifs.size(); i++) {
		ColladaMaterial& material = scene.materials.emplace_back();
		material.name = stringf("%d", i);
		material.surface = MaterialSurface(i);
		material.clamp_s = tie.ad_gifs.at(i).d4_clamp_1.data_lo;
		material.clamp_t = tie.ad_gifs.at(i).d4_clamp_1.data_hi;
	}
	
	for(s32 i = 0; i < (s32) tie.ad_gifs.size(); i++) {
		scene.texture_paths.emplace_back(stringf("%d.png", i));
	}
	
	Mesh& mesh = scene.meshes.emplace_back();
	mesh.name = "mesh";
	mesh.flags |= MESH_HAS_TEX_COORDS;
	
	for(const TiePacket& packet : tie.lods[0].packets) {
		for(const TiePrimitive& primitive : packet.primitives) {
			s32 base_vertex = (s32) mesh.vertices.size();
			
			SubMesh& submesh = mesh.submeshes.emplace_back();
			submesh.material = primitive.material_index;
			
			for(s32 i = 0; i < (s32) primitive.vertices.size(); i++) {
				Vertex& dest = mesh.vertices.emplace_back();
				const TieVertex& src = primitive.vertices[i];

				// The normals are stored in spherical coordinates, then there's a
				// cosine/sine lookup table at the top of the scratchpad.
				f32 normal_azimuth_radians = src.azimuth * (WRENCH_PI / 128.f);
				f32 normal_elevation_radians = src.elevation * (WRENCH_PI / 128.f);
				f32 cos_azimuth = cosf(normal_azimuth_radians);
				f32 sin_azimuth = sinf(normal_azimuth_radians);
				f32 cos_elevation = cosf(normal_elevation_radians);
				f32 sin_elevation = sinf(normal_elevation_radians);

				// This bit is done on VU0.
				f32 nx = cos_azimuth * cos_elevation;
				f32 ny = sin_azimuth * cos_elevation;
				f32 nz = sin_elevation;

				dest.pos.x = src.x * (tie.scale / 1024.f);
				dest.pos.y = src.y * (tie.scale / 1024.f);
				dest.pos.z = src.z * (tie.scale / 1024.f);
				dest.tex_coord.s = vu_fixed12_to_float(src.s);
				dest.tex_coord.t = vu_fixed12_to_float(src.t);
				dest.normal = glm::vec3(nx, ny, nz);
				
				if(i >= 2) {
					//if(i % 2 == primitive.face) {
						submesh.faces.emplace_back(base_vertex + i - 2, base_vertex + i - 1, base_vertex + i);
					//} else {
					//	submesh.faces.emplace_back(base_vertex + i - 0, base_vertex + i - 1, base_vertex + i - 2);
					//}
				}
			}
		}
	}

	for(size_t i = 0; i < scene.meshes.size(); i++) {
		fix_winding_orders_of_triangles_based_on_normals(scene.meshes[i]);
	}

	return scene;
}
