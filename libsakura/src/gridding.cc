/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2016
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
// Author: Kohji Nakamura <k.nakamura@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution

#include <cassert>
#include <cmath>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/localdef.h"
#include "libsakura/packed_type.h"

namespace {
#include "libsakura/packed_operation.h"

typedef int32_t integer;

template<typename T>
T Square(T n) {
	return n * n;
}

template<typename INT, typename REAL>
INT Round(REAL v) {
	return INT(v >= 0. ? v + REAL(0.5) : v - REAL(0.5));
}

// To overload conj(complex), the name starts with lower case
inline constexpr float conj(float v) {
	return v;
}

inline constexpr size_t At2(size_t B, integer a, integer b) {
	return a * B + b;
}
inline constexpr size_t At3(size_t B, size_t C, integer a, integer b, integer c) {
	return a * (B * C) + At2(C, b, c);
}
inline constexpr size_t At4(size_t B, size_t C, size_t D, integer a, integer b, integer c,
		integer d) {
	return a * (B * C * D) + At3(C, D, b, c, d);
}

inline void GridPosition(double const xy[/*2*/], integer sampling,
		integer loc[/*2*/], integer off[/*2*/]) {
	for (integer idim = 0; idim < 2; ++idim) {
		double pos = xy[idim];
		loc[idim] = Round<integer, double>(pos);
		off[idim] = Round<integer, double>((loc[idim] - pos) * sampling);
	}
}

inline bool OnGrid(double const xy[/*2*/], integer width, integer height,
		integer loc[/*2*/], integer support) {
	auto &x = loc[0];
	auto &y = loc[1];

	bool on_grid = (x - support >= 0) && (x + support < width)
			&& (y - support >= 0) && (y + support < height);

	// rough check against the case where integer can't represent the location correctly.
	on_grid = on_grid && (xy[0] - support >= -1)
			&& (xy[0] + support - 1 < width) && (xy[1] - support >= -1)
			&& (xy[1] + support - 1 < height);
	return on_grid;
}

struct WeightOnly {
	static inline constexpr float func(float weight, float value, float accumulator) {
		return weight + accumulator;
	}
};

struct WeightedValue {
	static inline constexpr float func(float weight, float value, float accumulator) {
		return weight * conj(value) + accumulator;
	}
};

struct VWeightOnly {
	static inline LIBSAKURA_SYMBOL(SimdPacketNative) func(
	LIBSAKURA_SYMBOL(SimdPacketNative) weight,
	LIBSAKURA_SYMBOL(SimdPacketNative) value,
	LIBSAKURA_SYMBOL(SimdPacketNative) accumulator) {
		return LIBSAKURA_SYMBOL(SimdMath)<LIBSAKURA_SYMBOL(SimdArchNative),
				float>::Add(weight, accumulator);
	}
};

struct VWeightedValue {
	static inline LIBSAKURA_SYMBOL(SimdPacketNative) func(
	LIBSAKURA_SYMBOL(SimdPacketNative) weight,
	LIBSAKURA_SYMBOL(SimdPacketNative) value,
	LIBSAKURA_SYMBOL(SimdPacketNative) accumulator) {
		return LIBSAKURA_SYMBOL(SimdMath)<LIBSAKURA_SYMBOL(SimdArchNative),
				float>::MulAdd(weight, value, accumulator);
	}
};

template<typename WeightFunc>
struct ScalarImpl {
	static inline void ChannelLoop(float convolution_factor,
			integer const num_channels,
			uint32_t const channel_map_arg[/*num_channels*/],
			float const weight_arg[/*num_channels*/],
			bool const mask[/*num_channels*/],
			float const value[/*num_channels*/], integer num_channels_for_grid,
			double weight_sum[/*num_channels_for_grid*/],
			float grid[/*num_channels_for_grid*/],
			float weight_of_grid[/*num_channels_for_grid*/]) {
		auto weight = AssumeAligned(weight_arg);
		auto channel_map = AssumeAligned(channel_map_arg);
		for (integer ichan = 0; ichan < num_channels; ++ichan) {
			if (mask[ichan]) {
				integer out_chan = channel_map[ichan];
				float the_weight = weight[ichan] * convolution_factor;
				grid[out_chan] = WeightFunc::func(the_weight, value[ichan], grid[out_chan]);
				weight_of_grid[out_chan] += the_weight;
				weight_sum[out_chan] += the_weight;
			}
		}
	}
};

#if defined(__AVX__)
template<typename WeightFunc, typename WeightFunc4Scalar>
struct VectorizedImpl {
	static inline void ChannelLoop(float convolution_factor,
			integer const num_channels,
			uint32_t const channel_map[/*num_channels*/],
			float const weight[/*num_channels*/],
			bool const mask[/*num_channels*/],
			float const value[/*num_channels*/], integer num_channels_for_grid,
			double weight_sum[/*num_channels_for_grid*/],
			float grid[/*num_channels_for_grid*/],
			float weight_of_grid[/*num_channels_for_grid*/]) {
		LIBSAKURA_SYMBOL(SimdPacketNative) pconvolution_factor;
		pconvolution_factor.set1(convolution_factor);

		auto vmask =
#if defined(__AVX2__)
				reinterpret_cast<__m64 const *>(mask);
#else
		reinterpret_cast<float const *>(mask);
#endif
		LIBSAKURA_SYMBOL(SimdPacketNative) const *vvalue =
				(LIBSAKURA_SYMBOL(SimdPacketNative) const *) value;

		LIBSAKURA_SYMBOL(SimdPacketNative) const *vweight =
				(LIBSAKURA_SYMBOL(SimdPacketNative) const *) weight;

		LIBSAKURA_SYMBOL(SimdPacketNative) *vgrid =
				(LIBSAKURA_SYMBOL(SimdPacketNative) *) grid;
		LIBSAKURA_SYMBOL(SimdPacketNative) *vweight_of_grid =
				(LIBSAKURA_SYMBOL(SimdPacketNative) *) weight_of_grid;

		LIBSAKURA_SYMBOL(SimdPacketNative) *vweight_sum =
				(LIBSAKURA_SYMBOL(SimdPacketNative) *) weight_sum;

#if defined(__AVX2__)
		__m128 const zero128i = _mm_setzero_ps();
#endif
		assert(
				num_channels % LIBSAKURA_SYMBOL(SimdPacketNative)::kNumFloat == 0);
		integer chan_limit = num_channels
				/ LIBSAKURA_SYMBOL(SimdPacketNative)::kNumFloat;
		for (integer ichan = 0; ichan < chan_limit; ++ichan) {
			LIBSAKURA_SYMBOL(SimdPacketNative) pweight;
#if defined(__AVX2__)
			pweight.raw_float = _mm256_cvtepi32_ps(
					_mm256_cvtepu8_epi32(
							_mm_castps_si128(_mm_loadl_pi(zero128i, vmask))));
			++vmask;
#else
			pweight.raw_float =
					_mm256_cvtepi32_ps(
							_mm256_castps_si256(
									_mm256_insertf128_ps(_mm256_castps128_ps256(
													_mm_castsi128_ps(_mm_cvtepi8_epi32(_mm_castps_si128(_mm_load1_ps(&vmask[0]))))),
											_mm_castsi128_ps(_mm_cvtepi8_epi32(_mm_castps_si128(_mm_load1_ps(&vmask[1])))),
											1)));
			vmask += 2;
#endif
			pweight = LIBSAKURA_SYMBOL(SimdMath)<
			LIBSAKURA_SYMBOL(SimdArchNative), float>::Mul(pweight,
					pconvolution_factor);
			pweight = LIBSAKURA_SYMBOL(SimdMath)<
			LIBSAKURA_SYMBOL(SimdArchNative), float>::Mul(pweight, *vweight);
			++vweight;
			*vgrid = WeightFunc::func(pweight, *vvalue, *vgrid);
			++vvalue;
			++vgrid;

			*vweight_of_grid = LIBSAKURA_SYMBOL(SimdMath)<
			LIBSAKURA_SYMBOL(SimdArchNative), float>::Add(*vweight_of_grid,
					pweight);
			++vweight_of_grid;

			using Math = LIBSAKURA_SYMBOL(SimdMath)<
			LIBSAKURA_SYMBOL(SimdArchNative), double>;
			using Converter = LIBSAKURA_SYMBOL(SimdConvert)<
			LIBSAKURA_SYMBOL(SimdArchNative)>;
			vweight_sum[0] = Math::Add(vweight_sum[0],
					Converter::FloatToDouble(pweight.v_prior.v[0]));
			vweight_sum[1] = Math::Add(vweight_sum[1],
					Converter::FloatToDouble(pweight.v_prior.v[1]));
			vweight_sum += 2;
		}
	}
};
#else
template<typename WeightFunc, typename WeightFunc4Scalar>
struct VectorizedImpl {
	static inline void ChannelLoop(
			float convolution_factor,
			integer const num_channels,
			uint32_t const channel_map[/*num_channels*/],
			float const weight[/*num_channels*/],
			bool const mask[/*num_channels*/],
			float const value[/*num_channels*/],
			integer num_channels_for_grid,
			double weight_sum[/*num_channels_for_grid*/],
			float grid[/*num_channels_for_grid*/],
			float weight_of_grid[/*num_channels_for_grid*/]
	) {
		assert(false);
	}
};

#endif

template<typename Impl>
inline void Grid(integer locx, integer locy, integer const doubled_support,
		integer const sampling, size_t const integral_radius[],
		integer const num_polarization,
		uint32_t const polarization_map_arg[/*num_polarization*/],
		integer const num_channels,
		uint32_t const channel_map[/*num_channels*/],
		bool const mask_arg/*[num_polarization]*/[/*num_channels*/],
		float const value_arg/*[num_polarization]*/[/*num_channels*/],
		float const weight[/*num_channels*/],
		integer num_convolution_table/*= sqrt(2)*support*sampling + extra*/,
		float const convolution_table_arg[/*num_convolution_table*/],
		integer num_polarization_for_grid, integer num_channels_for_grid,
		integer const width, integer const height,
		double weight_sum_arg/*[num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float weight_of_grid_arg/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float grid_arg/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/]) {

	auto polarization_map = AssumeAligned(polarization_map_arg);
	auto mask = AssumeAligned(mask_arg);
	auto value = AssumeAligned(value_arg);
	auto convolution_table = AssumeAligned(convolution_table_arg);
	auto weight_sum = AssumeAligned(weight_sum_arg);
	auto weight_of_grid = AssumeAligned(weight_of_grid_arg);
	auto grid = AssumeAligned(grid_arg);

	size_t ir = 0;

	for (integer iy = 0; iy < doubled_support; ++iy) {
		integer ay = locy + iy;
		for (integer ix = 0; ix < doubled_support; ++ix) {
			integer ax = locx + ix;
			float convolution_factor = convolution_table[integral_radius[ir]];
			float *grid_pos_local = &grid[At4(width, num_polarization_for_grid,
					num_channels_for_grid, ay, ax, 0, 0)];
			float *weight_of_grid_pos_local = &weight_of_grid[At4(width,
					num_polarization_for_grid, num_channels_for_grid, ay, ax, 0,
					0)];

			for (integer ipol = 0; ipol < num_polarization; ++ipol) {
				integer const out_pol = polarization_map[ipol];
				bool const *mask_local = &mask[At2(num_channels, ipol, 0)];
				float const *value_local = &value[At2(num_channels, ipol, 0)];

				double *weight_sum_local = &weight_sum[At2(
						num_channels_for_grid, out_pol, 0)];
				float *grid_local = &grid_pos_local[At2(num_channels_for_grid,
						out_pol, 0)];
				float *weight_of_grid_local = &weight_of_grid_pos_local[At2(
						num_channels_for_grid, out_pol, 0)];

				Impl::ChannelLoop(convolution_factor, num_channels, channel_map,
						weight, mask_local, value_local, num_channels_for_grid,
						weight_sum_local, grid_local, weight_of_grid_local);
			} // ipol
			++ir;
		} // ix
	} // iy
}

template<typename OptimizedImpl>
inline void InternalGrid(size_t num_spectra, size_t start_spectrum,
		size_t end_spectrum,
		bool const spectrum_mask[/*num_spectra*/],
		double const x_arg[/*num_spectra*/],
		double const y_arg[/*num_spectra*/], integer support, integer sampling,
		integer num_polarization,
		uint32_t const polarization_map[/*num_polarization*/],
		integer num_channels, uint32_t const channel_map[/*num_channels*/],
		bool const mask/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const value/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const weight/*[num_spectra]*/[/*num_channels*/],
		integer num_convolution_table/*= sqrt(2)*support*sampling + extra*/,
		float const convolution_table[/*num_convolution_table*/],
		integer num_polarization_for_grid, integer num_channels_for_grid,
		integer width, integer height,
		double weight_sum/*[num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float weight_of_grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/]) {

	auto x = AssumeAligned(x_arg);
	auto y = AssumeAligned(y_arg);
	integer const doubled_support = 2 * support + 1;
	for (size_t spectrum = start_spectrum; spectrum < end_spectrum;
			++spectrum) {
		if (spectrum_mask[spectrum]) {
			double xy[2] = { x[spectrum], y[spectrum] };
			integer point[2], offset[2];
			GridPosition(xy, sampling, point, offset);
			if (OnGrid(xy, width, height, point, support)) {
				SIMD_ALIGN
				size_t integral_radius[Square(doubled_support)];
				{
					float initial_relative_location_x = -(support + 1)
							* sampling + offset[0];
					float relative_location_y = -(support + 1) * sampling
							+ offset[1];
					size_t ir = 0;
					for (integer iy = 0; iy < doubled_support; ++iy) {
						relative_location_y += sampling;
						float relative_location_x = initial_relative_location_x;
						for (integer ix = 0; ix < doubled_support; ++ix) {
							relative_location_x += sampling;
							integral_radius[ir] = static_cast<size_t>(sqrt(
									Square(relative_location_x)
											+ Square(relative_location_y)));
							assert(integral_radius[ir] < num_convolution_table);
							++ir;
						}
					}
				}

				Grid<OptimizedImpl>(point[0] - support, point[1] - support,
						doubled_support, sampling, integral_radius,
						num_polarization, polarization_map, num_channels,
						channel_map,
						&mask[At3(num_polarization, num_channels, spectrum, 0,
								0)],
						&value[At3(num_polarization, num_channels, spectrum, 0,
								0)], &weight[At2(num_channels, spectrum, 0)],
						num_convolution_table, convolution_table,
						num_polarization_for_grid, num_channels_for_grid, width,
						height, weight_sum, weight_of_grid, grid);
			}
		}
	}
}

inline bool IsVectorOperationApplicable(size_t num_channels,
		uint32_t const channel_map_arg[/*num_channels*/]) {
#if !defined(__AVX__) || defined(ARCH_SCALAR)
	return false;
#endif
	constexpr size_t kElementsInPacket = LIBSAKURA_ALIGNMENT / sizeof(float);
	if (num_channels % kElementsInPacket != 0
			|| num_channels < kElementsInPacket) {
		return false;
	}

	auto channel_map = AssumeAligned(channel_map_arg);
	for (size_t i = 0; i < num_channels; ++i) {
		if (static_cast<size_t>(channel_map[i]) != i) {
			return false;
		}
	}
	return true;
}

void GridConvolvingCasted(size_t num_spectra, size_t start_spectrum,
		size_t end_spectrum,
		bool const spectrum_mask[/*num_spectra*/],
		double const x[/*num_spectra*/], double const y[/*num_spectra*/],
		integer support, integer sampling, integer num_polarization,
		uint32_t const polarization_map[/*num_polarization*/],
		integer num_channels, uint32_t const channel_map[/*num_channels*/],
		bool const mask/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const value/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const weight/*[num_spectra]*/[/*num_channels*/],
		bool weight_only,
		integer num_convolution_table/*= ceil(sqrt(2.)*(support+1)*sampling)*/,
		float const convolution_table[/*num_convolution_table*/],
		integer num_polarization_for_grid, integer num_channels_for_grid,
		integer width, integer height,
		double weight_sum/*[num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float weight_of_grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/]) {
	if (weight_only) {
		if (IsVectorOperationApplicable(num_channels, channel_map)) {
			InternalGrid<VectorizedImpl<VWeightOnly, WeightOnly> >(num_spectra,
					start_spectrum, end_spectrum, spectrum_mask, x, y, support,
					sampling, num_polarization, polarization_map, num_channels,
					channel_map, mask, value, weight, num_convolution_table,
					convolution_table, num_polarization_for_grid,
					num_channels_for_grid, width, height, weight_sum,
					weight_of_grid, grid);
		} else {
			InternalGrid<ScalarImpl<WeightOnly> >(num_spectra, start_spectrum,
					end_spectrum, spectrum_mask, x, y, support, sampling,
					num_polarization, polarization_map, num_channels,
					channel_map, mask, value, weight, num_convolution_table,
					convolution_table, num_polarization_for_grid,
					num_channels_for_grid, width, height, weight_sum,
					weight_of_grid, grid);
		}
	} else {
		if (IsVectorOperationApplicable(num_channels, channel_map)) {
			InternalGrid<VectorizedImpl<VWeightedValue, WeightedValue> >(num_spectra,
					start_spectrum, end_spectrum, spectrum_mask, x, y, support,
					sampling, num_polarization, polarization_map, num_channels,
					channel_map, mask, value, weight, num_convolution_table,
					convolution_table, num_polarization_for_grid,
					num_channels_for_grid, width, height, weight_sum,
					weight_of_grid, grid);
		} else {
			InternalGrid<ScalarImpl<WeightedValue> >(num_spectra,
					start_spectrum, end_spectrum, spectrum_mask, x, y, support,
					sampling, num_polarization, polarization_map, num_channels,
					channel_map, mask, value, weight, num_convolution_table,
					convolution_table, num_polarization_for_grid,
					num_channels_for_grid, width, height, weight_sum,
					weight_of_grid, grid);
		}
	}
}

inline void GridConvolving(size_t num_spectra,
		size_t start_spectrum, size_t end_spectrum,
		bool const spectrum_mask[/*num_spectra*/],
		double const x[/*num_spectra*/], double const y[/*num_spectra*/],
		size_t support, size_t sampling, size_t num_polarization,
		uint32_t const polarization_map[/*num_polarization*/],
		size_t num_channels, uint32_t const channel_map[/*num_channels*/],
		bool const mask/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const value/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const weight/*[num_spectra]*/[/*num_channels*/],
		bool weight_only,
		size_t num_convolution_table/*= ceil(sqrt(2.)*(support+1)*sampling)*/,
		float const convolution_table[/*num_convolution_table*/],
		size_t num_polarization_for_grid, size_t num_channels_for_grid,
		size_t width, size_t height,
		double weight_sum/*[num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float weight_of_grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/]) {
	GridConvolvingCasted(num_spectra, start_spectrum, end_spectrum,
			spectrum_mask, x, y, (integer) support, (integer) sampling,
			(integer) num_polarization, polarization_map,
			(integer) num_channels, channel_map, mask, value, weight,
			weight_only, (integer) num_convolution_table, convolution_table,
			(integer) num_polarization_for_grid,
			(integer) num_channels_for_grid, (integer) width, (integer) height,
			weight_sum, weight_of_grid, grid);
}

}

#define CHECK_ARGS(x) do { \
	if (!(x)) { \
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument); \
	} \
} while (false)

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GridConvolvingFloat)(
		size_t num_spectra, size_t start_spectrum, size_t end_spectrum,
		bool const spectrum_mask[/*num_spectra*/],
		double const x[/*num_spectra*/], double const y[/*num_spectra*/],
		size_t support, size_t sampling, size_t num_polarization,
		uint32_t const polarization_map[/*num_polarization*/],
		size_t num_channels, uint32_t const channel_map[/*num_channels*/],
		bool const mask/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const value/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const weight/*[num_spectra]*/[/*num_channels*/], bool weight_only,
		size_t num_convolution_table/*= ceil(sqrt(2.)*(support+1)*sampling)*/,
		float const convolution_table[/*num_convolution_table*/],
		size_t num_polarization_for_grid, size_t num_channels_for_grid,
		size_t width, size_t height,
		double weight_sum/*[num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float weight_of_grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/]) noexcept {
	CHECK_ARGS(spectrum_mask != nullptr);
	CHECK_ARGS(x != nullptr);
	CHECK_ARGS(y != nullptr);
	CHECK_ARGS(polarization_map != nullptr);
	CHECK_ARGS(channel_map != nullptr);
	CHECK_ARGS(mask != nullptr);
	CHECK_ARGS(weight_only || value != nullptr);
	CHECK_ARGS(weight != nullptr);
	CHECK_ARGS(convolution_table != nullptr);
	CHECK_ARGS(weight_sum != nullptr);
	CHECK_ARGS(weight_of_grid != nullptr);
	CHECK_ARGS(grid != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(spectrum_mask));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(x));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(y));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(polarization_map));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(channel_map));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(mask));
	CHECK_ARGS(weight_only || LIBSAKURA_SYMBOL(IsAligned)(value));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(weight));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(convolution_table));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(weight_sum));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(grid));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(weight_of_grid));
	CHECK_ARGS(
			0 <= start_spectrum && start_spectrum <= end_spectrum
					&& end_spectrum <= num_spectra);
	CHECK_ARGS(0 < support && support <= 256);
	CHECK_ARGS(0 < sampling && sampling <= INT32_MAX);
	CHECK_ARGS(static_cast<size_t>(support) * sampling <= INT32_MAX / 32);
	CHECK_ARGS(0 < num_polarization && num_polarization <= INT32_MAX);
	CHECK_ARGS(0 < num_channels && num_channels <= INT32_MAX);
	CHECK_ARGS(
			static_cast<int>(ceil(sqrt(2.) * (support + 1) * sampling)) <= num_convolution_table && num_convolution_table <= INT32_MAX / 32);
	CHECK_ARGS(
			0 < num_polarization_for_grid && num_polarization_for_grid <= INT32_MAX);
	CHECK_ARGS(0 < num_channels_for_grid && num_channels_for_grid <= INT32_MAX);
	CHECK_ARGS(0 < width && width <= INT32_MAX);
	CHECK_ARGS(0 < height && height <= INT32_MAX);
	CHECK_ARGS(static_cast<size_t>(width) * height <= INT32_MAX);
#if !defined(NDEBUG)
	for (size_t i = 0; i < num_polarization; ++i) {
		CHECK_ARGS(
				0 <= polarization_map[i]
				&& polarization_map[i] < num_polarization_for_grid);
	}
	for (size_t i = 0; i < num_channels; ++i) {
		CHECK_ARGS(
				0 <= channel_map[i] && channel_map[i] < num_channels_for_grid);
	}
#endif /* !defined(NDEBUG) */

	try {
		GridConvolving(num_spectra, start_spectrum, end_spectrum,
				spectrum_mask, x, y, support, sampling, num_polarization,
				polarization_map, num_channels, channel_map, mask, value,
				weight, weight_only, num_convolution_table, convolution_table,
				num_polarization_for_grid, num_channels_for_grid, width, height,
				weight_sum, weight_of_grid, grid);
	} catch (...) {
		assert(false); // No exception should be raised for the current implementation.
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}
