/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2014
 * National Astronomical Observatory of Japan
 * 2-21-1, Osawa, Mitaka, Tokyo, 181-8588, Japan.
 * 
 * This file is part of Sakura.
 * 
 * Sakura is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the 
 * Free Software Foundation, either version 3 of the License, or (at your 
 * option) any later version.
 * 
 * Sakura is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License 
 * along with Sakura.  If not, see <http://www.gnu.org/licenses/>.
 * @SAKURA_LICENSE_HEADER_END@
 */
#ifndef _LIBSAKURA_LOCATOR_TCC_
#define _LIBSAKURA_LOCATOR_TCC_

namespace {

template<class XDataType>
size_t LocateData(size_t start_position, size_t end_position, size_t num_base,
		XDataType const base_position[], XDataType located_position) {
	assert(num_base > 0);
	assert(start_position <= end_position);
	assert(end_position < num_base);
	assert(base_position != nullptr);

	// If length of the array is just 1, return 0
	if (num_base == 1)
		return 0;

	// base_position must be sorted in ascending order
	if (located_position <= base_position[0]) {
		// out of range
		return 0;
	} else if (located_position > base_position[num_base - 1]) {
		// out of range
		return num_base;
	} else if (located_position < base_position[start_position]) {
		// x_located is not in the range (start_position, end_position)
		// call this function to search other location
		return LocateData(0, start_position, num_base, base_position,
				located_position);
	} else if (located_position > base_position[end_position]) {
		// located_position is not in the range (start_position, end_position)
		// call this function to search other location
		return LocateData(end_position, num_base - 1, num_base, base_position,
				located_position);
	} else {
		// do bisection
		size_t left_index = start_position;
		size_t right_index = end_position;
		while (right_index > left_index + 1) {
			size_t middle_index = (right_index + left_index) / 2;
			if (located_position > base_position[middle_index]) {
				left_index = middle_index;
			} else {
				right_index = middle_index;
			}
		}
		return right_index;
	}
}

template<class XDataType>
size_t Locate(size_t num_base, size_t num_located,
		XDataType const base_position[], XDataType const located_position[],
		size_t location_list[]) {
	// input arrays must be sorted in ascending order
	size_t num_location_list = 0;
	if (num_base == 1) {
		if (located_position[num_located - 1] <= base_position[0]) {
			num_location_list = 1;
			location_list[0] = 0;
		} else if (located_position[0] >= base_position[0]) {
			num_location_list = 1;
			location_list[0] = 1;
		} else {
			num_location_list = 2;
			location_list[0] = 0;
			location_list[1] = 1;
		}
	} else if (located_position[num_located - 1] <= base_position[0]) {
		// all located_position's are on the left (lower) side of base_position
		num_location_list = 1;
		location_list[0] = 0;
	} else if (located_position[0] >= base_position[num_base - 1]) {
		// all located_position's are on the right (upper) side of base_position
		num_location_list = 1;
		location_list[0] = num_base;
	} else {
		size_t start_position = 0;
		size_t end_position = num_base - 1;
		size_t previous_location = num_base + 1;
		for (size_t i = 0; i < num_located && start_position <= end_position; ++i) {
			size_t location = LocateData(start_position, end_position, num_base,
					base_position, located_position[i]);
			if (location != previous_location) {
				location_list[num_location_list] = location;
				num_location_list += 1;
				start_position = location;
				previous_location = location;
			}
		}
	}
	return num_location_list;
}

}

#endif /* _LIBSAKURA_LOCATOR_TCC_ */
