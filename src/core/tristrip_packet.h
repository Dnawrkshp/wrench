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

#ifndef CORE_TRISTRIP_PACKET_H
#define CORE_TRISTRIP_PACKET_H

#include <core/mesh.h>
#include <core/mesh_graph.h>
#include <core/material.h>
#include <core/tristrip.h>

// This models the limited maximum size of a given packet. For each constraint,
// the number of the given objects in a packet will be multiplied by their
// respective costs and these results will be summed. If the sum is greater than
// the max cost, the packet is too big so this can be used to reject changes to
// a packet e.g. by limiting the length of a strip.
struct TriStripConstraints {
	u8 num_constraints = 0;
	s32 constant_cost[4];
	s32 strip_cost[4];
	s32 vertex_cost[4];
	s32 index_cost[4];
	s32 material_cost[4];
	s32 max_cost[4];
};

struct TriStripConfig {
	TriStripConstraints constraints;
	bool support_instancing;
};

struct TriStripRunningTotals {
	s32 strip_count = 0;
	s32 vertex_count = 0;
	s32 index_count = 0;
	s32 material_count = 0;
};

struct GeometryPacket {
	s32 primitive_begin = 0;
	s32 primitive_count = 0;
};

struct GeometryPackets {
	std::vector<GeometryPacket> packets;
	std::vector<GeometryPrimitive> primitives;
	std::vector<s32> indices;
};

GeometryPackets generate_tristrip_packets(const GeometryPrimitives& input, const std::vector<Material>& materials, const std::vector<EffectiveMaterial>& effectives, const TriStripConfig& config);

// Gets fed tristrips (as well as triangle lists) and incrementally splits them
// up into packets based on the constraints passed to it at construction time.
class TriStripPacketGenerator {
	const std::vector<Material>& _materials;
	const std::vector<EffectiveMaterial>& _effectives;
	const TriStripConstraints& _constraints;
	bool _support_instancing;
	
	TriStripRunningTotals _totals;
	GeometryPacket* _packet;
	s32 _current_effective_material = -1;
	
	GeometryPackets _output;
	
public:
	TriStripPacketGenerator(const std::vector<Material>& materials, const std::vector<EffectiveMaterial>& effectives, const TriStripConstraints& constraints, bool support_instancing);
	void add_list(const s32* indices, s32 index_count, s32 effective_material);
	void add_strip(const s32* indices, s32 index_count, s32 effective_material);
	GeometryPackets get_output();
	
private:
	void new_packet();
	bool try_add_strip();
	bool try_add_unique_vertex();
	bool try_add_duplicate_vertex();
	bool try_add_material();
};

#endif
