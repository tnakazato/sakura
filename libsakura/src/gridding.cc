// Author: Kohji Nakamura <k.nakamura@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution

#include <iostream>
#include <cstdio>
#include <cassert>
#include <libsakura/OptimizedImplementationFactoryImpl.h>
#include <libsakura/localdef.h>

using namespace std;
using namespace libsakura_PREFIX;

namespace {
typedef sakura::Gridding::integer integer;

template<typename T>
T square(T n) {
	return n * n;
}

template<typename INT, typename REAL>
INT nint(REAL v) {
	return INT(v >= 0. ? v + REAL(0.5) : v - REAL(0.5));
}

// overload conj(complex)
inline float conj(float v) {
	return v;
}

inline integer at2(integer B, integer a, integer b) {
	return a * B + b;
}
inline integer at3(integer B, integer C, integer a, integer b, integer c) {
	return a * (B * C) + at2(C, b, c);
}
inline integer at4(integer B, integer C, integer D, integer a, integer b,
		integer c, integer d) {
	return a * (B * C * D) + at3(C, D, b, c, d);
}

void gridPos(double const xy[/*2*/], integer sampling, integer loc[/*2*/],
		integer off[/*2*/]) {
	for (integer idim = 0; idim < 2; ++idim) {
		float pos = xy[idim];
		loc[idim] = nint<integer, float>(pos);
		off[idim] = nint<integer, float>((loc[idim] - pos) * sampling);
	}
}

bool onGrid(integer nx, integer ny, integer loc[/*2*/], integer support) {
	integer &x = loc[0];
	integer &y = loc[1];

	return (x - support >= 0) && (x + support < nx) && (y - support >= 0)
			&& (y + support < ny);
}

namespace ForSpace {
struct WeightOnly {
	static inline Gridding::value_t func(float weight,
			Gridding::value_t const values/*[nrow][nvischan]*/[/*nvispol*/],
			integer nvischan, integer nvispol, integer irow, integer ichan,
			integer ipol) {
		return Gridding::value_t(weight);
	}
};

struct WeightedValue {
	static inline Gridding::value_t func(float weight,
			Gridding::value_t const values/*[nrow][nvischan]*/[/*nvispol*/],
			integer nvischan, integer nvispol, integer irow, integer ichan,
			integer ipol) {
		return weight * conj(values[at3(nvischan, nvispol, irow, ichan, ipol)]);
	}
};

template<typename T>
inline void doGrid(
		Gridding::value_t const values/*[nrow][nvischan]*/[/*nvispol*/],
		integer nvispol, integer nvischan,
		Gridding::flag_t const flag/*[nrow][nvischan]*/[/*nvispol*/],
		float const weight/*[nrow]*/[/*nvischan*/],
		Gridding::value_t grid/*[nchan][npol][ny]*/[/*nx*/],
		float wgrid/*[nchan][npol][ny]*/[/*nx*/], integer nx, integer ny,
		integer npol, integer nchan, integer support, float const convTable[],
		integer const chanmap[/*nvischan*/], integer const polmap[/*nvispol*/],
		double sumwt/*[nchan]*/[/*npol*/], integer locx, integer locy,
		integer const irad[/*square(Wsupport)*/], integer const Wsupport,
		integer irow) {
	for (integer ichan = 0; ichan < nvischan; ++ichan) {
		integer achan = chanmap[ichan];
		float weight_ = weight[at2(nvischan, irow, ichan)];
		if (weight_ > 0.0) {
			for (integer ipol = 0; ipol < nvispol; ++ipol) {
				integer apol = polmap[ipol];
				if (!flag[at3(nvischan, nvispol, irow, ichan, ipol)]) {
					Gridding::value_t nvalue = T::func(weight_, values,
							nvischan, nvispol, irow, ichan, ipol);
					float norm = 0.0;
					integer ir = 0;
					// do iy=-support,support
					for (integer iy = 0; iy < Wsupport; ++iy) {
						integer ay = locy + iy;
						// do ix=-support,support
						for (integer ix = 0; ix < Wsupport; ++ix) {
							integer ax = locx + ix;
							float wt = convTable[irad[ir]];
							integer idx = at4(npol, ny, nx, achan, apol, ay,
									ax);
							grid[idx] += nvalue * wt;
							wgrid[idx] += weight_ * wt;
							norm += wt;
							ir++;
						} // ix
					} // iy
					sumwt[at2(npol, achan, apol)] += weight_ * norm;
				} // if
			} // ipol
		} // if
	} // ichan
}

typedef void
(*GridFunc_t)(Gridding::value_t const values/*[nrow][nvischan]*/[/*nvispol*/],
		integer nvispol, integer nvischan,
		Gridding::flag_t const flag/*[nrow][nvischan]*/[/*nvispol*/],
		float const weight/*[nrow]*/[/*nvischan*/],
		Gridding::value_t grid/*[nchan][npol][ny]*/[/*nx*/],
		float wgrid/*[nchan][npol][ny]*/[/*nx*/], integer nx, integer ny,
		integer npol, integer nchan, integer support, float const convTable[],
		integer const chanmap[/*nvischan*/], integer const polmap[/*nvispol*/],
		double sumwt/*[nchan]*/[/*npol*/], integer locx, integer locy,
		integer const irad[/*square(Wsupport)*/], integer const Wsupport,
		integer irow);

GridFunc_t gridFuncs[2] = { doGrid<WeightedValue>, doGrid<WeightOnly> };

struct Adapter {
	typedef GridFunc_t FuncType;

	static inline void check(integer nvispol, integer npol,
			integer const polmap[/*nvispol*/], integer nvischan, integer nchan,
			integer const chanmap[/*nvischan*/]) {
		assert(0 < npol);
		assert(0 < nvispol);
	}

	static inline FuncType chooseFunc(bool dowt, integer nvispol,
			integer npol) {
		return gridFuncs[dowt];
	}

	static inline void callFunc(FuncType func,
			Gridding::value_t const values/*[nrow][nvischan]*/[/*nvispol*/],
			integer nvispol, integer nvischan,
			Gridding::flag_t const flag/*[nrow][nvischan]*/[/*nvispol*/],
			float const weight/*[nrow]*/[/*nvischan*/],
			Gridding::value_t grid[], float wgrid[], integer nx, integer ny,
			integer npol, integer nchan, integer support,
			float const convTable[], integer const chanmap[/*nvischan*/],
			integer const polmap[/*nvispol*/],
			double sumwt/*[nchan]*/[/*npol*/], integer locx, integer locy,
			integer const irad[/*square(Wsupport)*/], integer const Wsupport,
			integer irow) {
		func(values, nvispol, nvischan, flag, weight, grid, wgrid, nx, ny, npol,
				nchan, support, convTable, chanmap, polmap, sumwt, locx, locy,
				irad, Wsupport, irow);
	}
};
} // ForSpace

namespace ForSpeed {
template<integer nvispol>
struct WeightOnly {
	static inline Gridding::value_t func(float weight,
			Gridding::value_t const values/*[nrow]*/[/*nvischan*/][nvispol],
			integer nvischan, integer irow, integer ichan, integer ipol) {
		return Gridding::value_t(weight);
	}
};

template<integer nvispol>
struct WeightedValue {
	static inline Gridding::value_t func(float weight,
			Gridding::value_t const values/*[nrow]*/[/*nvischan*/][nvispol],
			integer nvischan, integer irow, integer ichan, integer ipol) {
		return weight * conj(values[at2(nvischan, irow, ichan)][ipol]);
	}
};

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
				float weight_ = weight[at2(nvischan, irow, ichan)];
				bool nop = bool(weight_ > 0.0);
				integer idx = at3(nx, nchan, ay, ax, achan);
				for (integer ipol = 0; ipol < nvispol; ++ipol) {
					integer apol = polmap[ipol];
					float wt_ =
							wt
									* (nop
											&& (!flag[at2(nvischan, irow, ichan)][ipol]));
					Gridding::value_t nvalue = T::func(weight_, values,
							nvischan, irow, ichan, ipol);
#if 1
					grid[idx][apol] += nvalue * wt_;
					wgrid[idx][apol] += weight_ * wt_;
#else
					integer idx = at4(npol, ny, nx, achan, apol, ay, ax);
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

template <typename T, typename U>
void dump(T x) {
	union {
		T v;
		U f[4];
	}tmp;
	tmp.v = x;
	for (int i = 0; i < 4; i++) {
		cout << tmp.f[i] << ",";
	}
	cout << endl;
}

inline __m128 shuffleps(__m128 x, __m128i mask) {
	return _mm_castsi128_ps(_mm_shuffle_epi8(_mm_castps_si128(x), mask));
}

template <integer nvispol>
struct WeightOnlySIMD {
	static inline __m128
	func(__m128 weight,
			Gridding::value_t const values/*[nrow]*/[/*nvischan*/][nvispol],
			integer nvischan,
			integer irow, integer ichan) {
		return weight;
	}
};

template <integer nvispol>
struct WeightedValueSIMD {
	static inline __m128
	func(__m128 weight,
			Gridding::value_t const values/*[nrow]*/[/*nvischan*/][nvispol],
			integer nvischan,
			integer irow, integer ichan) {
		return _mm_mul_ps(weight,
				_mm_load_ps(values[at2(nvischan, irow, ichan)]));
	}
};

template <typename T>
inline void
doGridSIMD(Gridding::value_t const values/*[nrow]*/[/*nvischan*/][4],
		integer nvischan,
		Gridding::flag_t const flag/*[nrow]*/[/*nvischan*/][4],
		float const weight/*[nrow]*/[/*nvischan*/],
		Gridding::value_t grid/*[ny][nx]*/[/*nchan*/][4],
		float wgrid/*[ny][nx]*/[/*nchan*/][4],
		integer nx, integer ny, integer nchan,
		integer support,
		float const convTable[],
		integer const chanmap[/*nvischan*/],
		integer const polmap[/*nvispol*/],
		double sumwt[/*nchan*/][4],
		integer locx, integer locy,
		integer const irad[/*square(Wsupport)*/],
		integer const Wsupport,
		integer irow) {
	size_t const npol = 4;
	size_t const nvispol = 4;
	__m128 const zero = _mm_setzero_ps();
	__m128i const one = _mm_set1_epi32(1);
	__m128i const destMap[4] = {
		_mm_set_epi32( -1, -1, -1, 0x03020100 ),
		_mm_set_epi32( -1, -1, 0x03020100, -1 ),
		_mm_set_epi32( -1, 0x03020100, -1, -1 ),
		_mm_set_epi32( 0x03020100, -1, -1, -1 )
	};

	__m128i const src0 = _mm_set1_epi32(0x03020100);
	__m128i const src1 = _mm_set1_epi32(0x07060504);
	__m128i const src2 = _mm_set1_epi32(0x0B0A0908);
	__m128i const src3 = _mm_set1_epi32(0x0F0E0D0C);

	__m128 const shufMask0_ = _mm_castsi128_ps(_mm_shuffle_epi8(src0, destMap[polmap[0]]));
	__m128 const shufMask1_ = _mm_castsi128_ps(_mm_shuffle_epi8(src1, destMap[polmap[1]]));
	__m128 const shufMask2_ = _mm_castsi128_ps(_mm_shuffle_epi8(src2, destMap[polmap[2]]));
	__m128 const shufMask3_ = _mm_castsi128_ps(_mm_shuffle_epi8(src3, destMap[polmap[3]]));

	__m128i const shufMask0 = _mm_castps_si128(_mm_or_ps(_mm_cmpeq_ps(shufMask0_, zero), shufMask0_));
	__m128i const shufMask1 = _mm_castps_si128(_mm_or_ps(_mm_cmpeq_ps(shufMask1_, zero), shufMask1_));
	__m128i const shufMask2 = _mm_castps_si128(_mm_or_ps(_mm_cmpeq_ps(shufMask2_, zero), shufMask2_));
	__m128i const shufMask3 = _mm_castps_si128(_mm_or_ps(_mm_cmpeq_ps(shufMask3_, zero), shufMask3_));

	integer ir = 0;
	// do iy=-support,support
	for (integer iy=0; iy < Wsupport; ++iy) {
		integer ay = locy + iy;
		// do ix=-support,support
		for (integer ix=0; ix < Wsupport; ++ix) {
			integer ax = locx + ix;
			__m128 wt = _mm_set1_ps(convTable[irad[ir]]);
			for (integer ichan=0; ichan < nvischan; ++ichan) {
				integer achan = chanmap[ichan];
				integer idx = at3(nx, nchan, ay, ax, achan);
				float weight_ = weight[at2(nvischan, irow, ichan)];
				__m128 weight = _mm_set1_ps(weight_);
				__m128 mask = _mm_cmpgt_ps(weight, zero);
				__m128i flag_ =
				_mm_cvtepu8_epi32(_mm_castps_si128(
								_mm_load_ss((float*)flag[at2(nvischan, irow, ichan)])));
				mask = _mm_and_ps(mask,
						_mm_cmpeq_ps(_mm_castsi128_ps(flag_), zero));

				__m128 wt_ = _mm_and_ps(wt, mask);

				__m128 nvalue =
				T::func(weight, values, nvischan, irow, ichan);

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

template <>
inline void
doGrid<WeightOnly<4>, 4, 4>(Gridding::value_t const values/*[nrow]*/[/*nvischan*/][4],
		integer nvischan,
		Gridding::flag_t const flag/*[nrow]*/[/*nvischan*/][4],
		float const weight/*[nrow]*/[/*nvischan*/],
		Gridding::value_t grid/*[ny][nx]*/[/*nchan*/][4],
		float wgrid/*[ny][nx]*/[/*nchan*/][4],
		integer nx, integer ny, integer nchan,
		integer support,
		float const convTable[],
		integer const chanmap[/*nvischan*/],
		integer const polmap[/*nvispol*/],
		double sumwt[/*nchan*/][4],
		integer locx, integer locy,
		integer const irad[/*square(Wsupport)*/],
		integer const Wsupport,
		integer irow) {
	doGridSIMD<WeightOnlySIMD<4> >(
			values,
			nvischan,
			flag, weight,
			grid, wgrid,
			nx, ny,
			nchan,
			support, convTable,
			chanmap, polmap,
			sumwt,
			locx, locy,
			irad,
			Wsupport,
			irow);
}

template <>
inline void
doGrid<WeightedValue<4>, 4, 4>(Gridding::value_t const values/*[nrow]*/[/*nvischan*/][4],
		integer nvischan,
		Gridding::flag_t const flag/*[nrow]*/[/*nvischan*/][4],
		float const weight/*[nrow]*/[/*nvischan*/],
		Gridding::value_t grid/*[ny][nx]*/[/*nchan*/][4],
		float wgrid/*[ny][nx]*/[/*nchan*/][4],
		integer nx, integer ny, integer nchan,
		integer support,
		float const convTable[],
		integer const chanmap[/*nvischan*/],
		integer const polmap[/*nvispol*/],
		double sumwt[/*nchan*/][4],
		integer locx, integer locy,
		integer const irad[/*square(Wsupport)*/],
		integer const Wsupport,
		integer irow) {
	doGridSIMD<WeightedValueSIMD<4> >(
			values,
			nvischan,
			flag, weight,
			grid, wgrid,
			nx, ny,
			nchan,
			support, convTable,
			chanmap, polmap,
			sumwt,
			locx, locy,
			irad,
			Wsupport,
			irow);
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
FUNCSS(WeightedValue)
FUNCSS(WeightOnly)

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
#define FUNC(W,P,Q) reinterpret_cast<GridFunc_t>	\
      (FUNCTYPE(W,P,Q)(doGrid<W<P>, P, Q>))
#define FUNCS(W,P) { FUNC(W, P, 1), FUNC(W, P, 2),	\
		     FUNC(W, P, 3), FUNC(W, P, 4) }
#define FUNCSS(W) { FUNCS(W,1), FUNCS(W,2), FUNCS(W,3), FUNCS(W,4) }
		FUNCSS(WeightedValue),
		FUNCSS(WeightOnly)
	};
#undef FUNCSS
#undef FUNCS
#undef FUNC
#undef FUNCTYPE

struct Adapter {
	typedef GridFunc_t FuncType;

	static inline void check(integer nvispol, integer npol,
			integer const polmap[/*nvispol*/], integer nvischan, integer nchan,
			integer const chanmap[/*nvischan*/]) {
		assert(0 < npol && npol <= 4);
		assert(0 < nvispol && nvispol <= 4);
		for (integer i = 0; i < nvispol; ++i) {
			assert(0 <= polmap[i] && polmap[i] < npol);
		}
		for (integer i = 0; i < nvischan; ++i) {
			assert(0 <= chanmap[i] && chanmap[i] < nchan);
		}
	}

	static inline FuncType chooseFunc(bool dowt, integer nvispol,
			integer npol) {
		return gridFuncs[dowt][nvispol - 1][npol - 1];
	}

	static inline void callFunc(FuncType func,
			Gridding::value_t const values/*[nrow][nvischan]*/[/*nvispol*/],
			integer nvispol, integer nvischan,
			Gridding::flag_t const flag/*[nrow][nvischan]*/[/*nvispol*/],
			float const weight/*[nrow]*/[/*nvischan*/],
			Gridding::value_t grid[], float wgrid[], integer nx, integer ny,
			integer npol, integer nchan, integer support,
			float const convTable[], integer const chanmap[/*nvischan*/],
			integer const polmap[/*nvispol*/],
			double sumwt/*[nchan]*/[/*npol*/], integer locx, integer locy,
			integer const irad[/*square(Wsupport)*/], integer const Wsupport,
			integer irow) {
		func(values, nvischan, flag, weight, grid, wgrid, nx, ny, nchan,
				support, convTable, chanmap, polmap, sumwt, locx, locy, irad,
				Wsupport, irow);
	}
};
} // ForSpeed

template<typename T>
inline void internalGridsd(double const xy[/*nrow*/][2],
		Gridding::value_t const values/*[nrow][nvischan]*/[/*nvispol*/],
		integer nvispol, integer nvischan, bool dowt,
		Gridding::flag_t const flag/*[nrow][nvischan]*/[/*nvispol*/],
		integer const rflag[/*nrow*/],
		float const weight/*[nrow]*/[/*nvischan*/], integer nrow, integer irow,
		Gridding::value_t grid[], float wgrid[], integer nx, integer ny,
		integer npol, integer nchan, integer support, integer sampling,
		float const convTable[], integer const chanmap[/*nvischan*/],
		integer const polmap[/*nvispol*/], double sumwt/*[nchan]*/[/*npol*/]) {
	assert(sizeof(Gridding::flag_t) == 1);
	T::check(nvispol, npol, polmap, nvischan, nchan, chanmap);
	typename T::FuncType gridFunc = T::chooseFunc(dowt, nvispol, npol);

	integer rbeg = 0, rend = nrow;
	if (irow >= 0) {
		rbeg = irow;
		rend = irow + 1;
	}

	integer const Wsupport = 2 * support + 1;
	for (irow = rbeg; irow < rend; ++irow) {
		if (!rflag[irow]) {
			integer loc[2], off[2];
			gridPos(xy[irow], sampling, loc, off);
			if (onGrid(nx, ny, loc, support)) {
				integer irad[square(Wsupport)];
				{
					float rlocyInitial = -(support + 1) * sampling + off[0];
					float rlocy = -(support + 1) * sampling + off[1];
					integer ir = 0;
					for (integer iy = 0; iy < Wsupport; ++iy) {
						rlocy += sampling;
						float rlocx = rlocyInitial;
						for (integer ix = 0; ix < Wsupport; ++ix) {
							rlocx += sampling;
							irad[ir] = static_cast<integer>(sqrt(
									square(rlocx) + square(rlocy)));
							++ir;
						}
					}
				}

				T::callFunc(gridFunc, values, nvispol, nvischan, flag, weight,
						grid, wgrid, nx, ny, npol, nchan, support, convTable,
						chanmap, polmap, sumwt, loc[0] - support,
						loc[1] - support, irad, Wsupport, irow);
			}
		}
	}
}
}

namespace libsakura_PREFIX {
void ADDSUFFIX(Gridding, ARCH_SUFFIX)::gridsd(double const xy[/*nrow*/][2],
		Gridding::value_t const values/*[nrow][nvischan]*/[/*nvispol*/],
		integer nvispol, integer nvischan, bool dowt,
		Gridding::flag_t const flag/*[nrow][nvischan]*/[/*nvispol*/],
		integer const rflag[/*nrow*/],
		float const weight/*[nrow]*/[/*nvischan*/], integer nrow, integer irow,
		Gridding::value_t grid/*[nchan][npol][ny]*/[/*nx*/],
		float wgrid/*[nchan][npol][ny]*/[/*nx*/], integer nx, integer ny,
		integer npol, integer nchan, integer support, integer sampling,
		float const convTable[], integer const chanmap[/*nvischan*/],
		integer const polmap[/*nvispol*/],
		double sumwt/*[nchan]*/[/*npol*/]) const {
	internalGridsd<ForSpace::Adapter>(xy, values, nvispol, nvischan, dowt, flag,
			rflag, weight, nrow, irow, grid, wgrid, nx, ny, npol, nchan,
			support, sampling, convTable, chanmap, polmap, sumwt);
}
void ADDSUFFIX(Gridding, ARCH_SUFFIX)::gridsdForSpeed(
		double const xy[/*nrow*/][2],
		Gridding::value_t const values/*[nrow][nvischan]*/[/*nvispol*/],
		integer nvispol, integer nvischan, bool dowt,
		Gridding::flag_t const flag/*[nrow][nvischan]*/[/*nvispol*/],
		integer const rflag[/*nrow*/],
		float const weight/*[nrow]*/[/*nvischan*/], integer nrow, integer irow,
		Gridding::value_t grid/*[ny][nx][nchan]*/[/*npol*/],
		float wgrid/*[ny][nx][nchan]*/[/*npol*/], integer nx, integer ny,
		integer npol, integer nchan, integer support, integer sampling,
		float const convTable[], integer const chanmap[/*nvischan*/],
		integer const polmap[/*nvispol*/],
		double sumwt/*[nchan]*/[/*npol*/]) const {
	internalGridsd<ForSpeed::Adapter>(xy, values, nvispol, nvischan, dowt, flag,
			rflag, weight, nrow, irow, grid, wgrid, nx, ny, npol, nchan,
			support, sampling, convTable, chanmap, polmap, sumwt);
}

void ADDSUFFIX(Gridding, ARCH_SUFFIX)::transform(integer ny, integer nx,
		integer nchan, integer npol,
		Gridding::value_t gridFrom/*[ny][nx][nchan]*/[/*npol*/],
		float wgridFrom/*[ny][nx][nchan]*/[/*npol*/],
		Gridding::value_t gridTo/*[nchan][npol][ny]*/[/*nx*/],
		float wgridTo/*[nchan][npol][ny]*/[/*nx*/]) const {
	for (integer iy = 0; iy < ny; ++iy) {
		for (integer ix = 0; ix < nx; ++ix) {
			for (integer ichan = 0; ichan < nchan; ++ichan) {
				for (integer ipol = 0; ipol < npol; ++ipol) {
					gridTo[at4(npol, ny, nx, ichan, ipol, iy, ix)] =
							gridFrom[at4(nx, nchan, npol, iy, ix, ichan, ipol)];
					wgridTo[at4(npol, ny, nx, ichan, ipol, iy, ix)] =
							wgridFrom[at4(nx, nchan, npol, iy, ix, ichan, ipol)];
				}
			}
		}
	}
}
}
