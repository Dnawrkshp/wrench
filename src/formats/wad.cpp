/*
	wrench - A set of modding tools for the Ratchet & Clank PS2 games.
	Copyright (C) 2019-2020 chaoticgd

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

#include "wad.h"

#include <thread>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <algorithm>

// Enable/disable debug output for the decompression function.
#define WAD_DEBUG(cmd) //cmd
#define WAD_DEBUG_STATS(cmd) //cmd // Print statistics about packet counts/min and max values.
// If this code breaks, dump the correct output and point to that here.
//#define WAD_DEBUG_EXPECTED_PATH "<file path goes here>"

bool validate_wad(char* magic) {
	return std::memcmp(magic, "WAD", 3) == 0;
}

void decompress_wad(array_stream& dest, array_stream& src) {
	decompress_wad_n(dest, src, 0);
}

// We don't want to use stream::copy_n since it uses virtual functions.
void copy_bytes(array_stream& dest, array_stream& src, std::size_t bytes) {
	for(std::size_t i = 0; i < bytes; i++) {
		dest.write8(src.read8());
	}
}

void decompress_wad_n(array_stream& dest, array_stream& src, std::size_t bytes_to_decompress) {

	WAD_DEBUG(
		#ifdef WAD_DEBUG_EXPECTED_PATH
			file_stream expected(WAD_DEBUG_EXPECTED_PATH);
			std::optional<file_stream*> expected_ptr(&expected);
		#else
			std::optional<file_stream*> expected_ptr;
		#endif
	)
	
	WAD_DEBUG_STATS(
		uint32_t stat_literals_count = 0;
		uint32_t stat_literals_max_size = 0;
		
		size_t big_match_min_lookback_diff = UINT32_MAX;
		size_t big_match_max_lookback_diff = 0;
		int big_match_min_bytes = INT32_MAX;
		int big_match_max_bytes = INT32_MIN;
		size_t little_match_min_lookback_diff = UINT32_MAX;
		size_t little_match_max_lookback_diff = 0;
		int little_match_min_bytes = INT32_MAX;
		int little_match_max_bytes = INT32_MIN;
	)
	
	auto header = src.read<wad_header>(0);
	if(!validate_wad(header.magic)) {
		throw stream_format_error("Invalid WAD header.");
	}

	while(
		src.pos < header.total_size &&
		(bytes_to_decompress == 0 || dest.pos < bytes_to_decompress)) {

		WAD_DEBUG(
			dest.print_diff(expected_ptr);
			std::cout << "{dest.pos -> " << dest.pos << ", src.pos -> " << src.pos << "}\n\n";
		)

		WAD_DEBUG(
			static int count = 0;
			std::cout << "*** PACKET " << count++ << " ***\n";
		)

		uint8_t flag_byte = src.read8();
		WAD_DEBUG(std::cout << "flag_byte = " << std::hex << (flag_byte & 0xff) << "\n";)

		std::size_t lookback_offset = -1;
		int bytes_to_copy = 0;

		if(flag_byte < 0x10) { // Literal packet (0x0-0xf).
			uint32_t num_bytes;
			if(flag_byte != 0) {
				num_bytes = flag_byte + 3;
			} else {
				num_bytes = src.read8() + 18;
			}
			WAD_DEBUG(std::cout << " => copy 0x" << (int) num_bytes << " (decision) bytes from compressed stream at 0x" << src.pos << " to 0x" << dest.pos << "\n";)
			copy_bytes(dest, src, num_bytes);
			
			WAD_DEBUG_STATS(
				stat_literals_count++;
				stat_literals_max_size = std::max(num_bytes, stat_literals_max_size);
			)
			
			if(src.pos < src.buffer.size() && src.peek8() < 0x10) {
				// The game disallows this so lets complain.
				throw std::runtime_error("WAD decompression failed: Two literals in a row? Implausible!");
			}
			
			continue;
		} else if(flag_byte < 0x20) { // (0x10-0x1f)
			WAD_DEBUG(std::cout << " -- packet type C\n";)

			bytes_to_copy = flag_byte & 7;
			if(bytes_to_copy == 0) {
				bytes_to_copy = src.read8() + 7;
			}
			
			uint8_t b0 = src.read8();
			uint8_t b1 = src.read8();
			
			lookback_offset = dest.pos + ((flag_byte & 8) * -0x800 - ((b0 >> 2) + b1 * 0x40));
			if(lookback_offset != dest.pos) {
				bytes_to_copy += 2;
				lookback_offset -= 0x4000;
			} else if(bytes_to_copy != 1) {
				WAD_DEBUG(std::cout << " -- padding detected\n";)
				while(src.pos % 0x1000 != 0x10) {
					src.pos++;
				}
				continue;
			}
		} else if(flag_byte < 0x40) { // Big match packet (0x20-0x3f).
			WAD_DEBUG(std::cout << " -- packet type B\n";)

			bytes_to_copy = flag_byte & 0x1f;
			if(bytes_to_copy == 0) {
				bytes_to_copy = src.read8() + 0x1f;
			}
			bytes_to_copy += 2;

			uint8_t b1 = src.read8();
			uint8_t b2 = src.read8();
			lookback_offset = dest.pos - ((b1 >> 2) + b2 * 0x40) - 1;
			
			WAD_DEBUG_STATS(
				size_t lookback_diff = dest.pos - lookback_offset;
				big_match_min_lookback_diff = std::min(lookback_diff, big_match_min_lookback_diff);
				big_match_max_lookback_diff = std::max(lookback_diff, big_match_max_lookback_diff);
				big_match_min_bytes = std::min(bytes_to_copy, big_match_min_bytes);
				big_match_max_bytes = std::max(bytes_to_copy, big_match_max_bytes);
			)
		} else { // Little match packet (0x40-0xff).
			WAD_DEBUG(std::cout << " -- packet type A\n";)

			uint8_t b1 = src.read8();
			WAD_DEBUG(std::cout << " -- pos_major = " << (int) b1 << ", pos_minor = " << (int) ((flag_byte >> 2) & 7) << "\n";)
			lookback_offset = dest.pos - b1 * 8 - ((flag_byte >> 2) & 7) - 1;
			bytes_to_copy = (flag_byte >> 5) + 1;
			
			WAD_DEBUG_STATS(
				size_t lookback_diff = dest.pos - lookback_offset;
				little_match_min_lookback_diff = std::min(lookback_diff, little_match_min_lookback_diff);
				little_match_max_lookback_diff = std::max(lookback_diff, little_match_max_lookback_diff);
				little_match_min_bytes = std::min(bytes_to_copy, little_match_min_bytes);
				little_match_max_bytes = std::max(bytes_to_copy, little_match_max_bytes);
			)
		}

		if(bytes_to_copy != 1) {
			WAD_DEBUG(std::cout << " => copy 0x" << (int) bytes_to_copy << " bytes from uncompressed stream at 0x" << lookback_offset << "\n";)
			for(int i = 0; i < bytes_to_copy; i++) {
				dest.write8(dest.peek8(lookback_offset + i));
			}
		}
		
		uint32_t snd_pos = src.peek8(src.pos - 2) & 3;
		if(snd_pos != 0) {
			WAD_DEBUG(std::cout << " => copy 0x" << snd_pos << " (snd_pos) bytes from compressed stream at 0x" << src.pos << " to 0x" << dest.pos << "\n";)
			copy_bytes(dest, src, snd_pos);
			continue;
		}
	}

	WAD_DEBUG(std::cout << "Stopped reading at " << src.pos << "\n";)
	
	WAD_DEBUG_STATS(
		std::cout << "\n*** stats ***\n";
		std::cout << "literals count = " << stat_literals_count << "\n";
		std::cout << "literals max size = " << stat_literals_max_size << "\n";
		
		std::cout << "big_match_min_lookback_diff = " << big_match_min_lookback_diff << "\n";
		std::cout << "big_match_max_lookback_diff = " << big_match_max_lookback_diff << "\n";
		std::cout << "big_match_min_bytes = " << big_match_min_bytes << "\n";
		std::cout << "big_match_max_bytes = " << big_match_max_bytes << "\n";
		
		
		std::cout << "little_match_min_lookback_diff = " << little_match_min_lookback_diff << "\n";
		std::cout << "little_match_max_lookback_diff = " << little_match_max_lookback_diff << "\n";
		std::cout << "little_match_min_bytes = " << little_match_min_bytes << "\n";
		std::cout << "little_match_max_bytes = " << little_match_max_bytes << "\n";
	)
}

// Used for calculating the bounds of the sliding window.
std::size_t sub_clamped(std::size_t lhs, std::size_t rhs) {
	if(rhs > lhs) {
		return 0;
	}
	return lhs - rhs;
}

#define WAD_COMPRESS_DEBUG(cmd) //cmd
// NOTE: You must zero the size field of the expected file while debugging (0x3 through 0x6).
// WAD_COMPRESS_DEBUG_EXPECTED_PATH isn't actually very useful, as this encoder
// is NOT identical to the one used by Insomniac. It was useful when starting
// to write the algorithm, but probably won't be for debugging it.
//#define WAD_COMPRESS_DEBUG_EXPECTED_PATH "dumps/mobyseg_compressed_null"

void encode_wad_packet(
		array_stream& dest,
		const uint8_t* src,
		size_t& src_pos,
		size_t src_end,
		uint32_t& last_flag);

size_t get_wad_packet_size(uint8_t* src, size_t bytes_left);

static const int DO_NOT_INJECT_FLAG = 0x100;

static const size_t TYPE_A_MAX_LOOKBACK = 2044;
static const size_t MAX_LITERAL_SIZE = 255;//255 + 18;
static const size_t MAX_MATCH_SIZE = 0x100;

void compress_wad(array_stream& dest, array_stream& src, int thread_count) {
	WAD_COMPRESS_DEBUG(
		#ifdef WAD_COMPRESS_DEBUG_EXPECTED_PATH
			file_stream expected(WAD_COMPRESS_DEBUG_EXPECTED_PATH);
			std::optional<stream*> expected_ptr(&expected);
		#else
			std::optional<stream*> expected_ptr;
		#endif
	)
	
	std::vector<array_stream> intermediates(thread_count);
	
	// Compress the data into a stream of packets.
	if(thread_count == 1) {
		uint32_t last_flag = DO_NOT_INJECT_FLAG;
		
		src.seek(0);
		for(size_t i = 0; src.pos < src.buffer.size(); i++) {
			WAD_COMPRESS_DEBUG(
				std::cout << "{dest.pos -> " << dest.pos << ", src.pos -> " << src.pos << "}\n\n";
			)
			
			WAD_COMPRESS_DEBUG(
				std::cout << "*** PACKET " << i << " ***\n";
			)
			
			encode_wad_packet(intermediates[0], (uint8_t*) src.buffer.data(), src.pos, src.buffer.size(), last_flag);
			intermediates[0].pos = intermediates[0].buffer.size();
		}
	} else {
		size_t min_block_size = 0x100 * thread_count;
		size_t total_size = src.buffer.size();
		total_size += min_block_size - (total_size % min_block_size);
		size_t block_size = total_size / thread_count;
		
		std::vector<std::thread> threads(thread_count);
		for(int i = 0; i < thread_count; i++) {
			const uint8_t* src_ptr = (uint8_t*) src.buffer.data();
			size_t src_pos_val = block_size * i;
			size_t src_end = std::min(src.buffer.size(), (block_size * (i + 1)));
			threads[i] = std::thread([&intermediates, i, src_ptr, src_pos_val, src_end]() {
				size_t src_pos = src_pos_val; // Make non-const.
				uint32_t last_flag = DO_NOT_INJECT_FLAG;
				array_stream thread_dest;
				while(src_pos < src_end) {
					encode_wad_packet(thread_dest, src_ptr, src_pos, src_end, last_flag);
					thread_dest.pos = thread_dest.buffer.size();
				}
				intermediates[i].buffer = std::move(thread_dest.buffer);
			});
		}
		for(int i = 0; i < thread_count; i++) {
			threads[i].join();
		}
	}
	
	dest.seek(0);
	const char* header = "\x57\x41\x44\x00\x00\x00\x00\x57\x52\x45\x4e\x43\x48\x30\x31\x30";
	dest.write_n(header, 0x10);
	
	// Append the compressed data and insert padding where required.
	for(int i = 0; i < thread_count; i++) {
		array_stream& intermediate = intermediates[i];
		intermediate.pos = 0;
		while(intermediate.pos < intermediate.size()) {
			size_t packet_size = get_wad_packet_size(
				(uint8_t*) intermediate.buffer.data() + intermediate.pos,
				intermediate.buffer.size() - intermediate.pos);
			// The different blocks each thread generates may begin/end with
			// literal packets. Two consecutive literal packets aren't allowed,
			// so we add a dummy packet in between. We need to do this while
			// respecting the 0x2000 buffer size (see comment below), so we do
			// it here.
			bool insert_dummy = i != 0 && intermediate.pos == 0;
			size_t insert_size = packet_size;
			if(insert_dummy) {
				insert_size += 3;
			}
			// dest.pos is offset 0x10 bytes by the header:
			//  0x0000 WAD. .... .... ....
			//	0x0010 [data]
			//   ...
			//	0x2000 [data]
			//	0x2010 [start of new block]
			if(((dest.pos + 0x1ff0) % 0x2000) + insert_size > 0x2000 - 3) {
				// Every 0x2000 bytes or so there must be a pad packet or the
				// game crashes with a teq exception. This is because the game
				// copies the compressed data into the EE core's scratchpad, which
				// is 0x4000 bytes in size.
				dest.write8(0x12);
				dest.write8(0x0);
				dest.write8(0x0);
				while(dest.pos % 0x2000 != 0x10) {
					dest.write8(0xee);
				}
			}
			if(insert_dummy) {
				dest.write8(0x11);
				dest.write8(0);
				dest.write8(0);
			}
			dest.buffer.insert(dest.buffer.end(),
				intermediate.buffer.begin() + intermediate.pos,
				intermediate.buffer.begin() + intermediate.pos + packet_size);
			dest.pos += packet_size;
			intermediate.pos += packet_size;
		}
	}
	
	uint32_t total_size = dest.pos;
	dest.seek(3);
	dest.write<uint32_t>(total_size);
}

const std::vector<char> DUMMY_PACKET = { 0x11, 0, 0 };

void encode_wad_packet(
		array_stream& dest,
		const uint8_t* src, // Beginning of the whole buffer. Minimum lookback position.
		size_t& src_pos, // Position in the buffer to compress.
		size_t src_end, // End of the part of the buffer to compress.
		uint32_t& last_flag) {
	// Just emit literals at the end so we don't have to worry about overrunning
	// the input buffer.
	if(src_pos + MAX_MATCH_SIZE >= src_end) {
		WAD_COMPRESS_DEBUG(std::cout << "\n*** SPECIAL EOF PACKET ***\n");
		if(last_flag < 0x10 || last_flag == DO_NOT_INJECT_FLAG) {
			last_flag = 0x11;
			dest.buffer.insert(dest.buffer.end(), DUMMY_PACKET.begin(), DUMMY_PACKET.end());
			return;
		}
		size_t literal_size = std::min(src_end - src_pos, MAX_LITERAL_SIZE);
		if(literal_size <= 3) {
			dest.buffer[dest.pos - 2] |= literal_size; // Inject the literal size into the last packet.
		} else if(literal_size <= 18) {
			// We can encode the size in the flag byte.
			dest.buffer.push_back(literal_size - 3);
		} else {
			// We have to push it as a seperate byte (leave the flag as zero).
			dest.buffer.push_back(0);
			dest.buffer.push_back(literal_size - 18);
		}
		last_flag = DO_NOT_INJECT_FLAG;
		auto iter = src + src_pos;
		dest.buffer.insert(dest.buffer.end(), iter, iter + literal_size);
		src_pos += literal_size;
		return;
	}
	
	// Determine where the next repeating pattern is.
	size_t literal_size = MAX_LITERAL_SIZE;
	size_t match_offset = 0;
	size_t match_size = 0;
	const auto find_match = [&]() {
		literal_size = MAX_LITERAL_SIZE;
		match_offset = 0;
		match_size = 0;
		
		for(size_t i = 0; i < MAX_LITERAL_SIZE; i++) {
			size_t high = src_pos + i;
			size_t low = sub_clamped(high, TYPE_A_MAX_LOOKBACK);
			for(size_t j = low; j < high; j++) {
				// Count number of equal bytes.
				size_t target = src_pos + i;
				size_t k = 0;
				for(; k < MAX_MATCH_SIZE; k++) {
					auto l_val = src[target + k];
					auto r_val = src[j + k];
					if(l_val != r_val) {
						break;
					}
				}
				
				if(k >= 3 && k > match_size) {
					match_offset = j;
					match_size = k;
				}
			}
			if(match_size >= 3) {
				literal_size = i;
				break;
			}
		}
	};
	find_match();
	
	if(literal_size == 0) { // Match packet.
		std::size_t delta = src_pos - match_offset - 1;

		// Max bytes_to_copy for a packet of type A is 0x8.
		if(match_size <= 0x8) { // A type
			assert(match_size >= 3);
		
			WAD_COMPRESS_DEBUG(
				std::cout << "match at 0x" << std::hex << match_offset << " of size 0x" << match_size << "\n";
			)

			uint8_t pos_major = delta / 8;
			uint8_t pos_minor = delta % 8;

			WAD_COMPRESS_DEBUG(
				std::cout << "pos_major = " << (int) pos_major << ", pos_minor = " << (int) pos_minor << "\n";
			)

			dest.buffer.push_back(((match_size - 1) << 5) | (pos_minor << 2));
			dest.buffer.push_back(pos_major);
		} else { // B type
			WAD_COMPRESS_DEBUG(std::cout << "B type detected!\n";)
			
			if(match_size > (0b11111 + 2)) {
				dest.buffer.push_back(1 << 5); // flag
				dest.buffer.push_back(match_size - (0b11111 + 2));
			} else {
				dest.buffer.push_back((1 << 5) | (match_size - 2)); // flag
			}

			uint8_t pos_minor = delta % 0x40;
			uint8_t pos_major = delta / 0x40;
			
			dest.buffer.push_back(pos_minor << 2);
			dest.buffer.push_back(pos_major);
		}

		WAD_COMPRESS_DEBUG(
			std::cout << " => copy 0x" << std::hex << (int) match_size
				<< " bytes from uncompressed stream at 0x" << match_offset << " (source = " << src_pos << ")\n";
		)
		src_pos += match_size;
	} else { // Literal packet.
		if(last_flag < 0x10) { // Two literals in a row? Implausible!
			last_flag = 0x11;
			dest.buffer.insert(dest.buffer.end(), DUMMY_PACKET.begin(), DUMMY_PACKET.end());
			return;
		}
	
		if(literal_size <= 3) {
			// If the last flag is a literal, or there's already a small literal
			// injected into the last packet, we need to push a new dummy packet
			// that we can stuff a literal into on the next iteration.
			if(last_flag == DO_NOT_INJECT_FLAG) {
				last_flag = 0x11;
				dest.buffer.insert(dest.buffer.end(), DUMMY_PACKET.begin(), DUMMY_PACKET.end());
				return;
			}
			
			WAD_COMPRESS_DEBUG(
				std::cout << " => copy 0x" << std::hex << (int) literal_size
					<< " (snd_pos) bytes from compressed stream at 0x" << dest_pos + packet.size() << " (target = " << src_pos << ")\n";
			)
			dest.buffer[dest.pos - 2] |= literal_size;
			dest.buffer.insert(dest.buffer.end(), src + src_pos, src + src_pos + literal_size);
			src_pos += literal_size;
			last_flag = DO_NOT_INJECT_FLAG;
			return;
		} else if(literal_size <= 18) {
			// We can encode the size in the flag byte.
			dest.buffer.push_back(literal_size - 3); // flag
		} else {
			// We have to push it as a seperate byte.
			dest.buffer.push_back(0); // flag
			dest.buffer.push_back(literal_size - 18);
		}
		
		WAD_COMPRESS_DEBUG(
			std::cout << " => copy 0x" << std::hex << (int) literal_size
				<< " bytes from compressed stream at 0x" << dest_pos + packet.size() << " (target = " << src_pos << ")\n";
		)
		dest.buffer.insert(dest.buffer.end(), src + src_pos, src + src_pos + literal_size);
		src_pos += literal_size;
	}
	
	last_flag = dest.buffer[dest.pos];
	WAD_COMPRESS_DEBUG(std::cout << "flag_byte = " << std::hex << (packet[0] & 0xff) << "\n");
}

size_t get_wad_packet_size(uint8_t* src, size_t bytes_left) {
	size_t size_of_packet = 1; // flag
	uint8_t flag_byte = src[0];
	if(flag_byte < 0x10) { // Literal packet (0x0-0xf).
		if(flag_byte != 0) {
			size_of_packet += flag_byte + 3; // smalllit
		} else {
			size_of_packet += 1 + (src[1] + 18); // size + biglit
		}
		if(size_of_packet < bytes_left && src[size_of_packet] < 0x10) {
			throw std::runtime_error("Compression failed: Intermediate buffer corrupted (double literal)!\n");
		}
		return size_of_packet; // We can't put a tiny literal inside another literal.
	} else if(flag_byte < 0x20) { // (0x10-0x1f)
		uint8_t bytes_to_copy = flag_byte & 7;
		if(bytes_to_copy == 0) {
			size_of_packet++; // bytes_to_copy
		}
		size_of_packet += 1 + 1; // b0 + b1
	} else if(flag_byte < 0x40) { // Big match packet (0x20-0x3f).
		uint8_t bytes_to_copy = flag_byte & 0x1f;
		if(bytes_to_copy == 0) {
			size_of_packet++;
		}
		size_of_packet += 1 + 1; // b1 + b2
	} else { // Little match packet (0x40-0xff).
		size_of_packet++; // pos_major
	}
	size_of_packet += src[size_of_packet - 2] & 3; // Add on the tiny literal.
	return size_of_packet;
}
