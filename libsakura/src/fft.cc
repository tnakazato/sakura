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
#include <cassert>
#include <cstdlib>
#include <cstdint>

#if defined(__SSE2__) && !defined(ARCH_SCALAR)
#include "emmintrin.h"
#endif

#include "libsakura/sakura.h"
#include "libsakura/localdef.h"

namespace {
typedef union {
	float x;
} Type4;

typedef union {
	double x;
} Type8;

typedef union {
	struct {
		double x, y;
	} x;
} Type16;

#define RESTRICT __restrict

template<typename T>
struct LastDimFlip {
	static void Flip(size_t len, size_t dst_pos, T const *RESTRICT src,
			T *RESTRICT dst) {
		auto src_aligned = AssumeAligned<decltype(src), sizeof(T)>(src);
		auto dst_aligned = AssumeAligned<decltype(dst), sizeof(T)>(dst);
		assert(len >= dst_pos);
		size_t end = len - dst_pos;
		size_t i;
		for (i = 0; i < end; ++i, ++dst_pos) {
			dst_aligned[dst_pos] = src_aligned[i];
		}
		dst_pos = 0;
		for (; i < len; ++i, ++dst_pos) {
			dst_aligned[dst_pos] = src_aligned[i];
		}
	}
};

template<typename T>
struct LastDimNoFlip {
	static void Flip(size_t len, size_t dst_pos, T const *RESTRICT src,
			T *RESTRICT dst) {
		auto src_aligned = AssumeAligned<decltype(src), sizeof(T)>(src);
		auto dst_aligned = AssumeAligned<decltype(dst), sizeof(T)>(dst);
		for (size_t i = 0; i < len; ++i) {
			dst_aligned[i] = src_aligned[i];
		}
	}
};

#if defined(__SSE2__) && !defined(ARCH_SCALAR)
template<>
struct LastDimFlip<Type16> {
	static void Flip(size_t len, size_t dst_pos, Type16 const *RESTRICT src,
			Type16 *RESTRICT dst) {
		STATIC_ASSERT(sizeof(__m128d) == sizeof(*src));
		assert(len >= dst_pos);
		size_t end = len - dst_pos;
		size_t i;
		for (i = 0; i < end; ++i, ++dst_pos) {
			_mm_store_pd(reinterpret_cast<double*>(&dst[dst_pos]),
			_mm_load_pd(reinterpret_cast<double const*>(&src[i])));
		}
		dst_pos = 0;
		for (; i < len; ++i, ++dst_pos) {
			_mm_store_pd(reinterpret_cast<double*>(&dst[dst_pos]),
			_mm_load_pd(reinterpret_cast<double const*>(&src[i])));
		}
	}
};

template<>
struct LastDimNoFlip<Type16> {
	static void Flip(size_t len, size_t dst_pos, Type16 const *RESTRICT src,
			Type16 *RESTRICT dst) {
		STATIC_ASSERT(sizeof(__m128d) == sizeof(*src));
		for (size_t i = 0; i < len; ++i) {
			_mm_store_pd(reinterpret_cast<double*>(&dst[i]),
			_mm_load_pd(reinterpret_cast<double const*>(&src[i])));
		}
	}
};

#endif

template<typename T, typename LastDim, size_t kExtra>
void FlipLowLevel(size_t lower_plane_elements, size_t dim, size_t const len[],
		T const src[], T dst[]) {
	if (dim <= 0)
		return;
	size_t current_len = len[dim - 1];
	size_t dst_pos = (current_len + kExtra) / 2;
	if (dim == 1) {
		LastDim::Flip(current_len, dst_pos, src, dst);
	} else {
		size_t lower_elements = lower_plane_elements / len[dim - 2];
		for (size_t i = 0; i < current_len; ++i) {
			FlipLowLevel<T, LastDim, kExtra>(lower_elements, dim - 1, &len[0],
					&src[i * lower_plane_elements],
					&dst[dst_pos * lower_plane_elements]);
			dst_pos = (dst_pos + 1) % current_len;
		}
	}
}

template<typename T>
void Flip(bool reverse, bool inner_most_untouched, size_t dim,
		size_t const len[], T const src[], T dst[]) {
	STATIC_ASSERT(
			sizeof(T) == 1 || sizeof(T) == 2 || sizeof(T) == 4 || sizeof(T) == 8
					|| sizeof(T) == 16 || sizeof(T) == 32);
	size_t lower_plane_elements = 1;
	for (size_t i = 0; i + 1 < dim; ++i) {
		lower_plane_elements *= len[i];
	}
	if (reverse) {
		if (inner_most_untouched) {
			FlipLowLevel<T, LastDimNoFlip<T>, 1>(lower_plane_elements, dim, len,
					src, dst);
		} else {
			FlipLowLevel<T, LastDimFlip<T>, 1>(lower_plane_elements, dim, len,
					src, dst);
		}
	} else {
		if (inner_most_untouched) {
			FlipLowLevel<T, LastDimNoFlip<T>, 0>(lower_plane_elements, dim, len,
					src, dst);
		} else {
			FlipLowLevel<T, LastDimFlip<T>, 0>(lower_plane_elements, dim, len,
					src, dst);
		}
	}
}

} // namespace

#define CHECK_ARGS(x) do { \
	if (!(x)) { \
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument); \
	} \
} while (false)

namespace {

template<typename T>
LIBSAKURA_SYMBOL(Status) FlipMatrix(
bool reverse,
bool inner_most_untouched, size_t dims, size_t const elements[], T const src[],
		T dst[]) {
	CHECK_ARGS(elements != nullptr);
	CHECK_ARGS(src != nullptr);
	CHECK_ARGS(dst != nullptr);
	CHECK_ARGS(IsAligned(src, sizeof(T)));
	CHECK_ARGS(IsAligned(dst, sizeof(T)));

	try {
		STATIC_ASSERT(sizeof(src[0]) == sizeof(T));
		Flip<T>(reverse, inner_most_untouched, dims, elements,
				reinterpret_cast<T const *>(src), reinterpret_cast<T *>(dst));
	} catch (...) {
		assert(false); // No exception should be raised for the current implementation.
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

} // namespace

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(FlipArrayFloat)(
bool inner_most_untouched, size_t dims, size_t const elements[],
		float const src[], float dst[]) noexcept {
	typedef Type4 T;
	STATIC_ASSERT(sizeof(*src) == sizeof(T) && sizeof(*dst) == sizeof(T));
	return FlipMatrix<T>(false, inner_most_untouched, dims, elements,
			reinterpret_cast<T const *>(src), reinterpret_cast<T *>(dst));
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UnflipArrayFloat)(
bool inner_most_untouched, size_t dims, size_t const elements[],
		float const src[], float dst[]) noexcept {
	typedef Type4 T;
	STATIC_ASSERT(sizeof(*src) == sizeof(T) && sizeof(*dst) == sizeof(T));
	return FlipMatrix<T>(true, inner_most_untouched, dims, elements,
			reinterpret_cast<T const *>(src), reinterpret_cast<T *>(dst));
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(FlipArrayDouble)(
bool inner_most_untouched, size_t dims, size_t const elements[],
		double const src[], double dst[]) noexcept {
	typedef Type8 T;
	STATIC_ASSERT(sizeof(*src) == sizeof(T) && sizeof(*dst) == sizeof(T));
	return FlipMatrix<T>(false, inner_most_untouched, dims, elements,
			reinterpret_cast<T const *>(src), reinterpret_cast<T *>(dst));
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UnflipArrayDouble)(
bool inner_most_untouched, size_t dims, size_t const elements[],
		double const src[], double dst[]) noexcept {
	typedef Type8 T;
	STATIC_ASSERT(sizeof(*src) == sizeof(T) && sizeof(*dst) == sizeof(T));
	return FlipMatrix<T>(true, inner_most_untouched, dims, elements,
			reinterpret_cast<T const *>(src), reinterpret_cast<T *>(dst));
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(FlipArrayDouble2)(
bool inner_most_untouched, size_t dims, size_t const elements[],
		double const src[][2], double dst[][2]) noexcept {
	typedef Type16 T;
	STATIC_ASSERT(sizeof(*src) == sizeof(T) && sizeof(*dst) == sizeof(T));
	return FlipMatrix<T>(false, inner_most_untouched, dims, elements,
			reinterpret_cast<T const *>(src), reinterpret_cast<T *>(dst));
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UnflipArrayDouble2)(
bool inner_most_untouched, size_t dims, size_t const elements[],
		double const src[][2], double dst[][2]) noexcept {
	typedef Type16 T;
	STATIC_ASSERT(sizeof(*src) == sizeof(T) && sizeof(*dst) == sizeof(T));
	return FlipMatrix<T>(true, inner_most_untouched, dims, elements,
			reinterpret_cast<T const *>(src), reinterpret_cast<T *>(dst));
}
