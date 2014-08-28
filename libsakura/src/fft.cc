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
#include <cassert>
#include <cstdint>

#if defined(__SSE2__) && !defined(ARCH_SCALAR)
#include "emmintrin.h"
#endif

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/localdef.h"

namespace {
template<typename T>
struct LastDimFlip {
	static void flip(size_t len, size_t dstPos, T const src[], T dst[]) {
		auto srcAligned = AssumeAligned(src, sizeof(T));
		auto dstAligned = AssumeAligned(dst, sizeof(T));
		size_t i;
		size_t end = len - dstPos;
		for (i = 0; i < end; ++i, ++dstPos) {
			dstAligned[dstPos] = srcAligned[i];
		}
		dstPos = 0;
		for (; i < len; ++i, ++dstPos) {
			dstAligned[dstPos] = srcAligned[i];
		}
	}
};

template<typename T>
struct LastDimNoFlip {
	static void flip(size_t len, size_t dstPos, T const src[], T dst[]) {
		auto srcAligned = AssumeAligned(src, sizeof(T));
		auto dstAligned = AssumeAligned(dst, sizeof(T));
		for (size_t i = 0; i < len; ++i) {
			dstAligned[i] = srcAligned[i];
		}
	}
};

#if defined(__SSE2__) && !defined(ARCH_SCALAR)
template<>
struct LastDimFlip<LIBSAKURA_PREFIX::FFT::Type16> {
	static void flip(size_t len, size_t dstPos,
	LIBSAKURA_PREFIX::FFT::Type16 const src[],
	LIBSAKURA_PREFIX::FFT::Type16 dst[]) {
		STATIC_ASSERT(sizeof(__m128d) == sizeof(*src));
		size_t i;
		size_t end = len - dstPos;
		for (i = 0; i < end; ++i, ++dstPos) {
			_mm_store_pd(reinterpret_cast<double*>(&dst[dstPos]),
			_mm_load_pd(reinterpret_cast<double const*>(&src[i])));
		}
		dstPos = 0;
		for (; i < len; ++i, ++dstPos) {
			_mm_store_pd(reinterpret_cast<double*>(&dst[dstPos]),
			_mm_load_pd(reinterpret_cast<double const*>(&src[i])));
		}
	}
};

template<>
struct LastDimNoFlip<LIBSAKURA_PREFIX::FFT::Type16> {
	static void flip(size_t len, size_t dstPos,
			LIBSAKURA_PREFIX::FFT::Type16 const src[],
			LIBSAKURA_PREFIX::FFT::Type16 dst[]) {
		STATIC_ASSERT(sizeof(__m128d) == sizeof(*src));
		for (size_t i = 0; i < len; ++i) {
			_mm_store_pd(reinterpret_cast<double*>(&dst[i]),
			_mm_load_pd(reinterpret_cast<double const*>(&src[i])));
		}
	}
};

#endif

template<typename T, typename LastDim, size_t Extra>
void flip_(size_t lowerPlaneElements, size_t dim, size_t const len[],
		T const src[], T dst[]) {
	if (dim <= 0)
		return;
	size_t curLen = len[dim - 1];
	size_t dstPos = (curLen + Extra) / 2;
	if (dim == 1) {
		LastDim::flip(curLen, dstPos, src, dst);
	} else {
		size_t lowerElements = lowerPlaneElements / len[dim - 2];
		for (size_t i = 0; i < curLen; ++i) {
			flip_<T, LastDim, Extra>(lowerElements, dim - 1, &len[0],
					&src[i * lowerPlaneElements],
					&dst[dstPos * lowerPlaneElements]);
			dstPos = (dstPos + 1) % curLen;
		}
	}
}

template<typename T>
void flip(bool reverse, bool innerMostUntouched, size_t dim, size_t const len[],
		T const src[], T dst[]) {
	STATIC_ASSERT(
			sizeof(T) == 1 || sizeof(T) == 2 || sizeof(T) == 4 || sizeof(T) == 8
					|| sizeof(T) == 16 || sizeof(T) == 32);
	size_t lowerPlaneElements = 1;
	for (size_t i = 0; i + 1 < dim; ++i) {
		lowerPlaneElements *= len[i];
	}
	if (reverse) {
		if (innerMostUntouched) {
			flip_<T, LastDimNoFlip<T>, 1>(lowerPlaneElements, dim, len, src,
					dst);
		} else {
			flip_<T, LastDimFlip<T>, 1>(lowerPlaneElements, dim, len, src, dst);
		}
	} else {
		if (innerMostUntouched) {
			flip_<T, LastDimNoFlip<T>, 0>(lowerPlaneElements, dim, len, src,
					dst);
		} else {
			flip_<T, LastDimFlip<T>, 0>(lowerPlaneElements, dim, len, src, dst);
		}
	}
}

} // namespace

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(FFT, ARCH_SUFFIX)::Flip4(bool reverse, bool innerMostUntouched,
		size_t dim, size_t const len[], Type4 const src[], Type4 dst[]) const {
	flip<Type4>(reverse, innerMostUntouched, dim, len, src, dst);
}

void ADDSUFFIX(FFT, ARCH_SUFFIX)::Flip8(bool reverse, bool innerMostUntouched,
		size_t dim, size_t const len[], Type8 const src[], Type8 dst[]) const {
	flip<Type8>(reverse, innerMostUntouched, dim, len, src, dst);
}

void ADDSUFFIX(FFT, ARCH_SUFFIX)::Flip16(bool reverse, bool innerMostUntouched,
		size_t dim, size_t const len[], Type16 const src[],
		Type16 dst[]) const {
	flip<Type16>(reverse, innerMostUntouched, dim, len, src, dst);
}

} // namespace LIBSAKURA_PREFIX
