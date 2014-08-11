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
#ifndef _LIBSAKURA_INTERPOLATION_LINEAR_TCC_
#define _LIBSAKURA_INTERPOLATION_LINEAR_TCC_

namespace {

// TODO: documentation must be added
template<class XDataType, class YDataType>
class LinearXInterpolatorImpl {
public:
	void PrepareForInterpolation(uint8_t polynomial_order, size_t num_base,
			size_t num_array, XDataType const base_position[],
			YDataType const base_data[]) {
	}
	void Interpolate1D(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t num_location,
			size_t const location[], size_t offset) {
		for (size_t k = 1; k < num_location; ++k) {
			size_t left_index = offset + k - 1;
			for (size_t j = 0; j < num_array; ++j) {
				size_t offset_index_left = j * num_base + left_index;
				XDataType dydx =
						static_cast<XDataType>(base_data[offset_index_left + 1]
								- base_data[offset_index_left])
								/ (base_position[left_index + 1]
										- base_position[left_index]);
				YDataType base_term = base_data[offset_index_left];
				YDataType *y_work = &interpolated_data[j * num_interpolated];
				for (size_t i = location[k - 1]; i < location[k]; ++i) {
					y_work[i] = base_term
							+ static_cast<YDataType>(dydx
									* (interpolated_position[i]
											- base_position[left_index]));
				}
			}
		}
	}
};

template<class XDataType, class YDataType>
class LinearYInterpolatorImpl {
public:
	void PrepareForInterpolation(uint8_t polynomial_order, size_t num_base,
			size_t num_array, XDataType const base_position[],
			YDataType const base_data[]) {
	}
	void Interpolate1D(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t num_location,
			size_t const location[], size_t offset) {
		for (size_t k = 1; k < num_location; ++k) {
			size_t left_index = offset + k - 1;
			for (size_t i = location[k - 1]; i < location[k]; ++i) {
				YDataType fraction =
						static_cast<YDataType>((interpolated_position[i]
								- base_position[left_index])
								/ (base_position[left_index + 1]
										- base_position[left_index]));
				YDataType *work = &interpolated_data[num_array * i];
				YDataType const *left_value = &base_data[num_array * left_index];
				YDataType const *right_value = &base_data[num_array
						* (left_index + 1)];
				for (size_t j = 0; j < num_array; ++j) {
					work[j] = left_value[j]
							+ fraction * (right_value[j] - left_value[j]);
				}
			}
		}
	}
};

}

#endif /* _LIBSAKURA_INTERPOLATION_LINEAR_TCC_ */
