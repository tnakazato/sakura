// Author: Kohji Nakamura <k.nakamura@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution

#include <iostream>
#include <cstdio>
#include <cassert>

#include "libsakura/optimized_implementation_factory_impl.h"

using ::LIBSAKURA_PREFIX::Gridding;
using ::std::cout;
using ::std::endl;

namespace {
#include "libsakura/localdef.h"
#include "libsakura/packed_operation.h"

typedef Gridding::integer integer;

template<typename T>
T Square(T n) {
	return n * n;
}

template<typename INT, typename REAL>
INT Round(REAL v) {
	return INT(v >= 0. ? v + REAL(0.5) : v - REAL(0.5));
}

// To overload conj(complex), the name starts with lower case
inline float conj(float v) {
	return v;
}

inline integer At2(integer B, integer a, integer b) {
	return a * B + b;
}
inline integer At3(integer B, integer C, integer a, integer b, integer c) {
	return a * (B * C) + At2(C, b, c);
}
inline integer At4(integer B, integer C, integer D, integer a, integer b,
		integer c, integer d) {
	return a * (B * C * D) + At3(C, D, b, c, d);
}

inline void GridPosition(double const xy[/*2*/], integer sampling, integer loc[/*2*/],
		integer off[/*2*/]) {
	for (integer idim = 0; idim < 2; ++idim) {
		float pos = xy[idim];
		loc[idim] = Round<integer, float>(pos);
		off[idim] = Round<integer, float>((loc[idim] - pos) * sampling);
	}
}

inline bool OnGrid(integer width, integer height, integer loc[/*2*/], integer support) {
	integer &x = loc[0];
	integer &y = loc[1];

	return (x - support >= 0) && (x + support < width) && (y - support >= 0)
			&& (y + support < height);
}

namespace ForSpeed {
struct WeightOnly {
	static inline float func(float weight,
			float value) {
		return weight;
	}
};

struct WeightedValue {
	static inline float func(float weight,
			float value) {
		return weight * conj(value);
	}
};

struct VWeightOnly {
	static inline LIBSAKURA_SYMBOL(SimdPacketNative) func(LIBSAKURA_SYMBOL(SimdPacketNative) weight,
			LIBSAKURA_SYMBOL(SimdPacketNative) value) {
		return weight;
	}
};

struct VWeightedValue {
	static inline LIBSAKURA_SYMBOL(SimdPacketNative) func(LIBSAKURA_SYMBOL(SimdPacketNative) weight,
			LIBSAKURA_SYMBOL(SimdPacketNative) value) {
		return LIBSAKURA_SYMBOL(SimdMath)<LIBSAKURA_SYMBOL(SimdArchNative), float>::Mul(weight, value);
	}
};

#if defined(__AVX__)
template<typename WeightFunc>
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
		LIBSAKURA_SYMBOL(SimdPacketNative) pconvolution_factor;
		pconvolution_factor.set1(convolution_factor);
		LIBSAKURA_SYMBOL(SimdPacketMMX) const*vmask =
				(LIBSAKURA_SYMBOL(SimdPacketMMX) const*)mask;

		LIBSAKURA_SYMBOL(SimdPacketNative) const *vvalue =
				(LIBSAKURA_SYMBOL(SimdPacketNative) const *)value;

		LIBSAKURA_SYMBOL(SimdPacketNative) const *vweight =
				(LIBSAKURA_SYMBOL(SimdPacketNative) const *)weight;
		for (integer ichan = 0;
				ichan < num_channels
								/ LIBSAKURA_SYMBOL(SimdPacketNative)::kNumFloat;
				++ichan) {

			LIBSAKURA_SYMBOL(SimdPacketNative) pweight;
			pweight.v_prior[0] = LIBSAKURA_SYMBOL(SimdConvert)::byteToFloat(vmask[0]);
			pweight.v_prior[1] = LIBSAKURA_SYMBOL(SimdConvert)::byteToFloat(vmask[1]);
			vmask+=2;
			pweight = LIBSAKURA_SYMBOL(SimdMath)<
					LIBSAKURA_SYMBOL(SimdArchNative), float>::Mul(pweight,
					pconvolution_factor);
			pweight = LIBSAKURA_SYMBOL(SimdMath)<
					LIBSAKURA_SYMBOL(SimdArchNative), float>::Mul(pweight,
					*vweight);
			vweight++;
			LIBSAKURA_SYMBOL(SimdPacketNative) pvalue = WeightFunc::func(pweight, *vvalue);
			vvalue++;
#if 0
			for (int i = 0; i < LIBSAKURA_SYMBOL(SimdPacketNative)::kNumFloat;
					i++) {
				uint32_t out_chan = channel_map[i];
				grid[out_chan] += pvalue.v_float.v[i];
				weight_of_grid[out_chan] += pweight.v_float.v[i];
				weight_sum[out_chan] += pweight.v_float.v[i];
			}
#else
#endif
			channel_map += LIBSAKURA_SYMBOL(SimdPacketNative)::kNumFloat;
		}
	}
};
#else
template<typename WeightFunc>
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
	}
	static inline void ChannelLoop_(
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
		/* READONLY
		 */
		for (integer ichan = 0; ichan < num_channels; ++ichan) {
			integer out_chan = channel_map[ichan];
			float conv_factor = mask[ichan] ? convolution_factor : 0.;
			float the_weight = weight[ichan] * conv_factor;
			float nvalue = WeightFunc::func(the_weight, value[ichan]);
			grid[out_chan] += nvalue;
			weight_of_grid[out_chan] += the_weight;
			weight_sum[out_chan] += the_weight;
		}
	}
};

#endif

template<typename WeightFunc>
struct ScalarImpl {
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
		/* READONLY
		 */
		for (integer ichan = 0; ichan < num_channels; ++ichan) {
			integer out_chan = channel_map[ichan];
			float conv_factor = mask[ichan] ? convolution_factor : 0.;
			float the_weight = weight[ichan] * conv_factor;
			float nvalue = WeightFunc::func(the_weight, value[ichan]);
			grid[out_chan] += nvalue;
			weight_of_grid[out_chan] += the_weight;
			weight_sum[out_chan] += the_weight;
		}
	}
};

template <typename Impl>
inline void Grid(
		integer locx, integer locy,
		integer const doubled_support, integer const sampling,
		integer const integral_radius[],
		integer const num_polarization,
		uint32_t const polarization_map[/*num_polarization*/],
		integer const num_channels,
		uint32_t const channel_map[/*num_channels*/],
		bool const mask/*[num_polarization]*/[/*num_channels*/],
		float const value/*[num_polarization]*/[/*num_channels*/],
		float const weight[/*num_channels*/],
		integer num_convolution_table/*= sqrt(2)*support*sampling + extra*/,
		float const convolution_table[/*num_convolution_table*/],
		integer num_polarization_for_grid, integer num_channels_for_grid,
		integer const width, integer const height,
		double weight_sum/*[num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float weight_of_grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/]) {

	integer const pixels = width * height;
	integer ir = 0;

	for (integer iy = 0; iy < doubled_support; ++iy) {
		integer ay = locy + iy;
		for (integer ix = 0; ix < doubled_support; ++ix) {
			integer ax = locx + ix;
			float convolution_factor = convolution_table[integral_radius[ir]];
			integer const point_idx = At2(width, ay, ax);
			float *grid_pos_local = &grid[At4(width, num_polarization_for_grid, num_channels_for_grid, ay, ax, 0, 0)];
			float *weight_of_grid_pos_local = &weight_of_grid[At4(width, num_polarization_for_grid, num_channels_for_grid, ay, ax, 0, 0)];

			for (integer ipol = 0; ipol < num_polarization; ++ipol) {
				integer const out_pol = polarization_map[ipol];
				bool const *mask_local = &mask[At2(num_channels, ipol, 0)];
				float const *value_local = &value[At2(num_channels, ipol, 0)];

				double *weight_sum_local = &weight_sum[At2(num_channels_for_grid, out_pol, 0)];
				float *grid_local = &grid_pos_local[At2(num_channels_for_grid, out_pol, 0)];
				float *weight_of_grid_local = &weight_of_grid_pos_local[At2(num_channels_for_grid, out_pol, 0)];

				Impl::ChannelLoop(
						convolution_factor,
						num_channels,
						channel_map, weight, mask_local, value_local,
						num_channels_for_grid,
						weight_sum_local, grid_local, weight_of_grid_local);
			} // ipol
			ir++;
		} // ix
	} // iy
}

template<typename T, integer nvispol, integer npol>
inline void doGrid(
		Gridding::value_t const values/*[nrow]*/[/*nvischan*/][nvispol],
		integer nvischan,
		Gridding::flag_t const flag/*[nrow]*/[/*nvischan*/][nvispol],
		float const weight/*[nrow]*/[/*nvischan*/],
		Gridding::value_t grid/*[ny][nx]*/[/*nchan*/][npol],
		float wgrid/*[ny][nx]*/[/*nchan*/][npol], integer nx, integer ny,
		integer nchan, integer support, float const convTable[],
		integer const chanmap[/*nvischan*/], integer const polmap[/*nvispol*/],
		double sumwt[/*nchan*/][npol], integer locx, integer locy,
		integer const irad[/*square(Wsupport)*/], integer const Wsupport,
		integer irow) {
	integer ir = 0;
	// do iy=-support,support
	for (integer iy = 0; iy < Wsupport; ++iy) {
		integer ay = locy + iy;
		// do ix=-support,support
		for (integer ix = 0; ix < Wsupport; ++ix) {
			integer ax = locx + ix;
			float wt = convTable[irad[ir]];
			for (integer ichan = 0; ichan < nvischan; ++ichan) {
				integer achan = chanmap[ichan];
				float weight_ = weight[At2(nvischan, irow, ichan)];
				bool nop = bool(weight_ > 0.0);
				integer idx = At3(nx, nchan, ay, ax, achan);
				for (integer ipol = 0; ipol < nvispol; ++ipol) {
					integer apol = polmap[ipol];
					float wt_ =
							wt
									* (nop
											&& (!flag[At2(nvischan, irow, ichan)][ipol]));
					Gridding::value_t nvalue = T::func(weight_, values,
							nvischan, irow, ichan, ipol);
#if 1
					grid[idx][apol] += nvalue * wt_;
					wgrid[idx][apol] += weight_ * wt_;
#else
					integer idx = At4(npol, ny, nx, achan, apol, ay, ax);
					((Gridding::value_t*)grid)[idx] += nvalue * wt_;
					((float *)wgrid)[idx] += weight_ * wt_;
#endif
					sumwt[achan][apol] += weight_ * wt_;
				} // ipol
			} // ichan
			ir++;
		} // ix
	} // iy
}

#if defined( __AVX__)
#include <immintrin.h>

template<typename T, typename U>
void dump(T x) {
	union {
		T v;
		U f[4];
	} tmp;
	tmp.v = x;
	for (int i = 0; i < 4; i++) {
		cout << tmp.f[i] << ",";
	}
	cout << endl;
}

inline __m128 shuffleps(__m128 x, __m128i mask) {
	return _mm_castsi128_ps(_mm_shuffle_epi8(_mm_castps_si128(x), mask));
}

template<integer nvispol>
struct WeightOnlySIMD {
	static inline __m128 func(__m128 weight,
			Gridding::value_t const values/*[nrow]*/[/*nvischan*/][nvispol],
			integer nvischan, integer irow, integer ichan) {
		return weight;
	}
};

template<integer nvispol>
struct WeightedValueSIMD {
	static inline __m128 func(__m128 weight,
			Gridding::value_t const values/*[nrow]*/[/*nvischan*/][nvispol],
			integer nvischan, integer irow, integer ichan) {
		return _mm_mul_ps(weight,
				_mm_load_ps(values[At2(nvischan, irow, ichan)]));
	}
};

template<typename T>
inline void doGridSIMD(
		Gridding::value_t const values/*[nrow]*/[/*nvischan*/][4],
		integer nvischan,
		Gridding::flag_t const flag/*[nrow]*/[/*nvischan*/][4],
		float const weight/*[nrow]*/[/*nvischan*/],
		Gridding::value_t grid/*[ny][nx]*/[/*nchan*/][4],
		float wgrid/*[ny][nx]*/[/*nchan*/][4], integer nx, integer ny,
		integer nchan, integer support, float const convTable[],
		integer const chanmap[/*nvischan*/], integer const polmap[/*nvispol*/],
		double sumwt[/*nchan*/][4], integer locx, integer locy,
		integer const irad[/*square(Wsupport)*/], integer const Wsupport,
		integer irow) {
	size_t const npol = 4;
	size_t const nvispol = 4;
	__m128  const zero = _mm_setzero_ps();
	__m128i  const one = _mm_set1_epi32(1);
	__m128i  const destMap[4] = { _mm_set_epi32(-1, -1, -1, 0x03020100),
			_mm_set_epi32(-1, -1, 0x03020100, -1), _mm_set_epi32(-1, 0x03020100,
					-1, -1), _mm_set_epi32(0x03020100, -1, -1, -1) };

	__m128i  const src0 = _mm_set1_epi32(0x03020100);
	__m128i  const src1 = _mm_set1_epi32(0x07060504);
	__m128i  const src2 = _mm_set1_epi32(0x0B0A0908);
	__m128i  const src3 = _mm_set1_epi32(0x0F0E0D0C);

	__m128  const shufMask0_ = _mm_castsi128_ps(
			_mm_shuffle_epi8(src0, destMap[polmap[0]]));
	__m128  const shufMask1_ = _mm_castsi128_ps(
			_mm_shuffle_epi8(src1, destMap[polmap[1]]));
	__m128  const shufMask2_ = _mm_castsi128_ps(
			_mm_shuffle_epi8(src2, destMap[polmap[2]]));
	__m128  const shufMask3_ = _mm_castsi128_ps(
			_mm_shuffle_epi8(src3, destMap[polmap[3]]));

	__m128i  const shufMask0 = _mm_castps_si128(
			_mm_or_ps(_mm_cmpeq_ps(shufMask0_, zero), shufMask0_));
	__m128i  const shufMask1 = _mm_castps_si128(
			_mm_or_ps(_mm_cmpeq_ps(shufMask1_, zero), shufMask1_));
	__m128i  const shufMask2 = _mm_castps_si128(
			_mm_or_ps(_mm_cmpeq_ps(shufMask2_, zero), shufMask2_));
	__m128i  const shufMask3 = _mm_castps_si128(
			_mm_or_ps(_mm_cmpeq_ps(shufMask3_, zero), shufMask3_));

	integer ir = 0;
	// do iy=-support,support
	for (integer iy = 0; iy < Wsupport; ++iy) {
		integer ay = locy + iy;
		// do ix=-support,support
		for (integer ix = 0; ix < Wsupport; ++ix) {
			integer ax = locx + ix;
			__m128 wt = _mm_set1_ps(convTable[irad[ir]]);
			for (integer ichan = 0; ichan < nvischan; ++ichan) {
				integer achan = chanmap[ichan];
				integer idx = At3(nx, nchan, ay, ax, achan);
				float weight_ = weight[At2(nvischan, irow, ichan)];
				__m128 weight = _mm_set1_ps(weight_);
				__m128 mask = _mm_cmpgt_ps(weight, zero);
				__m128i flag_ =
						_mm_cvtepu8_epi32(
								_mm_castps_si128(
										_mm_load_ss(
												(float*) flag[At2(nvischan,
														irow, ichan)])));
				mask = _mm_and_ps(mask,
						_mm_cmpeq_ps(_mm_castsi128_ps(flag_), zero));

				__m128 wt_ = _mm_and_ps(wt, mask);

				__m128 nvalue = T::func(weight, values, nvischan, irow, ichan);

				{
					//wgrid[idx] += weight * wt_;
					__m128 tmp = _mm_mul_ps(weight, wt_);
					__m128 sum = shuffleps(tmp, shufMask0);
					sum = _mm_add_ps(sum, shuffleps(tmp, shufMask1));
					sum = _mm_add_ps(sum, shuffleps(tmp, shufMask2));
					sum = _mm_add_ps(sum, shuffleps(tmp, shufMask3));
					_mm_store_ps(wgrid[idx],
							_mm_add_ps(sum, _mm_load_ps(wgrid[idx])));

					//sumwt[achan] += weight * wt_;
					_mm256_store_pd(sumwt[achan],
							_mm256_add_pd(_mm256_cvtps_pd(sum),
									_mm256_load_pd(sumwt[achan])));
				}
				{
					//grid[idx] += nvalue * wt_;
					__m128 tmp = _mm_mul_ps(nvalue, wt_);
					__m128 sum = shuffleps(tmp, shufMask0);
					sum = _mm_add_ps(sum, shuffleps(tmp, shufMask1));
					sum = _mm_add_ps(sum, shuffleps(tmp, shufMask2));
					sum = _mm_add_ps(sum, shuffleps(tmp, shufMask3));
					_mm_store_ps(grid[idx],
							_mm_add_ps(sum, _mm_load_ps(grid[idx])));
				}
			} // ichan
			ir++;
		} // ix
	} // iy
}

#endif

#define FUNCTYPE(W,nvispol,npol)			\
    GridFunc_ ## W ## _ ## nvispol ## _ ## npol ## _t
#define FUNCTYPEDEF(W,nvispol,npol)		\
    typedef  void (*FUNCTYPE(W, nvispol, npol))				\
    (Gridding::value_t const values/*[nrow]*/[/*nvischan*/][nvispol],	\
     integer nvischan,							\
     Gridding::flag_t const flag/*[nrow]*/[/*nvischan*/][nvispol],	\
     float const weight/*[nrow]*/[/*nvischan*/],			\
     Gridding::value_t grid/*[ny][nx]*/[/*nchan*/][npol],		\
     float wgrid/*[ny][nx]*/[/*nchan*/][npol],				\
     integer nx, integer ny, integer nchan,				\
     integer support,							\
     float const convTable[],						\
     integer const chanmap[/*nvischan*/],				\
     integer const polmap[/*nvispol*/],					\
     double sumwt[/*nchan*/][npol],					\
     integer locx, integer locy,					\
     integer const irad[/*square(Wsupport)*/],				\
     integer const Wsupport,						\
     integer irow);
#define FUNC(W,P,Q) FUNCTYPEDEF(W, P, Q)
#define FUNCS(W,P) FUNC(W, P, 1) FUNC(W, P, 2) FUNC(W, P, 3) FUNC(W, P, 4)
#define FUNCSS(W) FUNCS(W,1) FUNCS(W,2) FUNCS(W,3) FUNCS(W,4)
FUNCSS(WeightedValue)
FUNCSS(WeightOnly)
#undef FUNCTYPEDEF
#undef FUNC

#define FUNC(W,nvispol,npol)  FUNCTYPE(W, nvispol, npol)		\
    gridFunc_ ## W ## _ ## nvispol ## _ ## npol = doGrid<W<nvispol>, nvispol, npol>;
/* To instantiate template functions */

#undef FUNC
#undef FUNCS
#undef FUNCSS

typedef void (*GridFunc_t)(
		Gridding::value_t const values/*[nrow][nvischan]*/[/*nvispol*/],
		integer nvischan,
		Gridding::flag_t const flag/*[nrow][nvischan]*/[/*nvispol*/],
		float const weight/*[nrow]*/[/*nvischan*/],
		Gridding::value_t grid/*[ny][nx][nchan]*/[/*npol*/],
		float wgrid/*[ny][nx][nchan]*/[/*npol*/], integer nx, integer ny,
		integer nchan, integer support, float const convTable[],
		integer const chanmap[/*nvischan*/], integer const polmap[/*nvispol*/],
		double sumwt/*[nchan]*/[/*npol*/], integer locx, integer locy,
		integer const irad[/*square(Wsupport)*/], integer const Wsupport,
		integer irow);

GridFunc_t gridFuncs[2][/*nvispol*/4][/*npol*/4] = {
	};
#undef FUNCSS
#undef FUNCS
#undef FUNC
#undef FUNCTYPE

struct Adapter {
	typedef GridFunc_t FuncType;

	static inline FuncType chooseFunc(bool do_weight, integer num_channels) {

	}

	static inline void callFunc(FuncType func,
			float const values/*[nrow][nvischan]*/[/*nvispol*/],
			integer nvispol, integer nvischan,
			Gridding::flag_t const flag/*[nrow][nvischan]*/[/*nvispol*/],
			float const weight/*[nrow]*/[/*nvischan*/],
			float grid[], float wgrid[], integer nx, integer ny,
			integer npol, integer nchan, integer support,
			float const conv_table[], integer const chanmap[/*nvischan*/],
			integer const polmap[/*nvispol*/],
			double sumwt/*[nchan]*/[/*npol*/], integer locx, integer locy,
			integer const irad[/*square(Wsupport)*/], integer const Wsupport,
			integer irow) {
		func(values, nvischan, flag, weight, grid, wgrid, nx, ny, nchan,
				support, conv_table, chanmap, polmap, sumwt, locx, locy, irad,
				Wsupport, irow);
	}
};
} // ForSpeed

template<typename OptimizedImpl>
inline void InternalGrid(
		integer num_spectra,
		integer start_spectrum, integer end_spectrum,
		bool const spectrum_mask[/*num_spectra*/],
		double const x[/*num_spectra*/],
		double const y[/*num_spectra*/],
		integer support, integer sampling,
		integer num_polarization,
		uint32_t const polarization_map[/*num_polarization*/],
		integer num_channels,
		uint32_t const channel_map[/*num_channels*/],
		bool const mask/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const value/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const weight/*[num_spectra]*/[/*num_channels*/],
		integer num_convolution_table/*= sqrt(2)*support*sampling + extra*/,
		float const convolution_table[/*num_convolution_table*/],
		integer num_polarization_for_grid, integer num_channels_for_grid,
		integer width, integer height,
		double weight_sum/*[num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float weight_of_grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/]
		) {
	integer const doubled_support = 2 * support + 1;
	for (int spectrum = start_spectrum; spectrum < end_spectrum; ++spectrum) {
		if (spectrum_mask[spectrum]) {
			double xy[2] = { x[spectrum], y[spectrum]};
			integer point[2], offset[2];
			GridPosition(xy, sampling, point, offset);
			if (OnGrid(width, height, point, support)) {
				integer integral_radius[Square(doubled_support)];
				{
					float initial_relative_location_x = -(support + 1) * sampling + offset[0];
					float relative_location_y = -(support + 1) * sampling + offset[1];
					integer ir = 0;
					for (integer iy = 0; iy < doubled_support; ++iy) {
						relative_location_y += sampling;
						float relative_location_x = initial_relative_location_x;
						for (integer ix = 0; ix < doubled_support; ++ix) {
							relative_location_x += sampling;
							integral_radius[ir] = static_cast<integer>(sqrt(
									Square(relative_location_x) + Square(relative_location_y)));
							assert(integral_radius[ir] < num_convolution_table);
							++ir;
						}
					}
				}

				ForSpeed::Grid<OptimizedImpl>(
						point[0] - support, point[1] - support,
						doubled_support, sampling,
						integral_radius,
						num_polarization,
						polarization_map,
						num_channels,
						channel_map,
						&mask[At3(num_polarization, num_channels, spectrum, 0, 0)],
						&value[At3(num_polarization, num_channels, spectrum, 0, 0)],
						&weight[At2(num_channels, spectrum, 0)],
						num_convolution_table,
						convolution_table,
						num_polarization_for_grid, num_channels_for_grid,
						width, height,
						weight_sum,
						weight_of_grid,
						grid);
			}
		}
	}
}

inline bool IsVectorOperationApplicable(int num_channels) {
	size_t elements_in_packet = LIBSAKURA_SYMBOL(GetAlignment)() / sizeof(float);
	if (num_channels % elements_in_packet == 0
			&& num_channels >= elements_in_packet * 2) {
		return true;
	}
	return false;
}

void GridConvolvingCasted(integer num_spectra,
		integer start_spectrum, integer end_spectrum,
		bool const spectrum_mask[/*num_spectra*/],
		double const x[/*num_spectra*/],
		double const y[/*num_spectra*/],
		integer support, integer sampling,
		integer num_polarization,
		uint32_t const polarization_map[/*num_polarization*/],
		integer num_channels,
		uint32_t const channel_map[/*num_channels*/],
		bool const mask/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const value/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const weight/*[num_spectra]*/[/*num_channels*/],
		bool do_weight,
		integer num_convolution_table/*= ceil(sqrt(2.)*(support+1)*sampling)*/,
		float const convolution_table[/*num_convolution_table*/],
		integer num_polarization_for_grid, integer num_channels_for_grid,
		integer width, integer height,
		double weight_sum/*[num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float weight_of_grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/]
		) {
	if (do_weight) {
		if (IsVectorOperationApplicable(num_channels)) {
			InternalGrid<ForSpeed::VectorizedImpl<ForSpeed::VWeightOnly> >(num_spectra,
					start_spectrum, end_spectrum,
					spectrum_mask,
					x, y,
					support, sampling,
					num_polarization,
					polarization_map,
					num_channels,
					channel_map,
					mask,
					value,
					weight,
					num_convolution_table,
					convolution_table,
					num_polarization_for_grid, num_channels_for_grid,
					width, height,
					weight_sum,
					weight_of_grid,
					grid);
		} else {
			InternalGrid<ForSpeed::ScalarImpl<ForSpeed::WeightOnly> >(num_spectra,
					start_spectrum, end_spectrum,
					spectrum_mask,
					x, y,
					support, sampling,
					num_polarization,
					polarization_map,
					num_channels,
					channel_map,
					mask,
					value,
					weight,
					num_convolution_table,
					convolution_table,
					num_polarization_for_grid, num_channels_for_grid,
					width, height,
					weight_sum,
					weight_of_grid,
					grid);
		}
	} else {
		if (IsVectorOperationApplicable(num_channels)) {
			InternalGrid<ForSpeed::VectorizedImpl<ForSpeed::VWeightedValue> >(num_spectra,
					start_spectrum, end_spectrum,
					spectrum_mask,
					x, y,
					support, sampling,
					num_polarization,
					polarization_map,
					num_channels,
					channel_map,
					mask,
					value,
					weight,
					num_convolution_table,
					convolution_table,
					num_polarization_for_grid, num_channels_for_grid,
					width, height,
					weight_sum,
					weight_of_grid,
					grid);
		} else {
			InternalGrid<ForSpeed::ScalarImpl<ForSpeed::WeightedValue> >(num_spectra,
					start_spectrum, end_spectrum,
					spectrum_mask,
					x, y,
					support, sampling,
					num_polarization,
					polarization_map,
					num_channels,
					channel_map,
					mask,
					value,
					weight,
					num_convolution_table,
					convolution_table,
					num_polarization_for_grid, num_channels_for_grid,
					width, height,
					weight_sum,
					weight_of_grid,
					grid);
		}
	}
}
}

namespace LIBSAKURA_PREFIX {

void ADDSUFFIX(Gridding, ARCH_SUFFIX)::GridConvolving(size_t num_spectra,
		size_t start_spectrum, size_t end_spectrum,
		bool const spectrum_mask[/*num_spectra*/],
		double const x[/*num_spectra*/],
		double const y[/*num_spectra*/],
		size_t support, size_t sampling,
		size_t num_polarization,
		uint32_t const polarization_map[/*num_polarization*/],
		size_t num_channels,
		uint32_t const channel_map[/*num_channels*/],
		bool const mask/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const value/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const weight/*[num_spectra]*/[/*num_channels*/],
		bool do_weight,
		size_t num_convolution_table/*= ceil(sqrt(2.)*(support+1)*sampling)*/,
		float const convolution_table[/*num_convolution_table*/],
		size_t num_polarization_for_grid, size_t num_channels_for_grid,
		size_t width, size_t height,
		double weight_sum/*[num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float weight_of_grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/]
		) const {
	GridConvolvingCasted((integer)num_spectra,
			(integer)start_spectrum, (integer)end_spectrum,
			spectrum_mask,
			x, y,
			(integer)support, (integer)sampling,
			(integer)num_polarization,
			polarization_map,
			(integer)num_channels,
			channel_map,
			mask,
			value,
			weight,
			do_weight,
			(integer)num_convolution_table,
			convolution_table,
			(integer)num_polarization_for_grid, (integer)num_channels_for_grid,
			(integer)width, (integer)height,
			weight_sum,
			weight_of_grid,
			grid
			);
}

#if 0
void ADDSUFFIX(Gridding, ARCH_SUFFIX)::Transform(integer ny, integer nx,
		integer nchan, integer npol,
		Gridding::value_t grid_from/*[ny][nx][nchan]*/[/*npol*/],
		float wgrid_from/*[ny][nx][nchan]*/[/*npol*/],
		Gridding::value_t gridTo/*[nchan][npol][ny]*/[/*nx*/],
		float wgridTo/*[nchan][npol][ny]*/[/*nx*/]) const {
	for (integer iy = 0; iy < ny; ++iy) {
		for (integer ix = 0; ix < nx; ++ix) {
			for (integer ichan = 0; ichan < nchan; ++ichan) {
				for (integer ipol = 0; ipol < npol; ++ipol) {
					gridTo[At4(npol, ny, nx, ichan, ipol, iy, ix)] =
							grid_from[At4(nx, nchan, npol, iy, ix, ichan, ipol)];
					wgridTo[At4(npol, ny, nx, ichan, ipol, iy, ix)] =
							wgrid_from[At4(nx, nchan, npol, iy, ix, ichan, ipol)];
				}
			}
		}
	}
}
#endif
} // LIBSAKURA_PREFIX
