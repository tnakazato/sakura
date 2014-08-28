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
#include <cstdlib>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

#define CHECK_ARGS(x) do { \
	if (!(x)) { \
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument); \
	} \
} while (false)

namespace {

template<typename T>
LIBSAKURA_SYMBOL(Status) FlipMatrix(
		void (LIBSAKURA_PREFIX::FFT::*func)(bool reverse,
		bool innerMostUntouched, size_t dim, size_t const len[], T const src[],
				T dst[]) const,
		bool reverse,
		bool innerMostUntouched, size_t dims, size_t const elements[],
		T const src[], T dst[]) {
	CHECK_ARGS(elements != nullptr);
	CHECK_ARGS(src != nullptr);
	CHECK_ARGS(dst != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(src));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(dst));

	try {
		auto fft =
				::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetFFTImpl();
		STATIC_ASSERT(sizeof(src[0]) == sizeof(T));
		(fft->*func)(reverse, innerMostUntouched, dims, elements,
				reinterpret_cast<T const *>(src), reinterpret_cast<T *>(dst));
	} catch (...) {
		assert(false); // no exception should not be raised for the current implementation.
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

} // namespace

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(FlipMatrixFloat)(
bool innerMostUntouched, size_t dims, size_t const elements[],
		float const src[], float dst[]) {
	typedef LIBSAKURA_PREFIX::FFT::Type4 T;
	STATIC_ASSERT(sizeof(*src) == sizeof(T) && sizeof(*dst) == sizeof(T));
	return FlipMatrix<T>(&LIBSAKURA_PREFIX::FFT::Flip4, false,
			innerMostUntouched, dims, elements,
			reinterpret_cast<T const *>(src), reinterpret_cast<T *>(dst));
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UnflipMatrixFloat)(
bool innerMostUntouched, size_t dims, size_t const elements[],
		float const src[], float dst[]) {
	typedef LIBSAKURA_PREFIX::FFT::Type4 T;
	STATIC_ASSERT(sizeof(*src) == sizeof(T) && sizeof(*dst) == sizeof(T));
	return FlipMatrix<T>(&LIBSAKURA_PREFIX::FFT::Flip4, true,
			innerMostUntouched, dims, elements,
			reinterpret_cast<T const *>(src), reinterpret_cast<T *>(dst));
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(FlipMatrixDouble)(
bool innerMostUntouched, size_t dims, size_t const elements[],
		double const src[], double dst[]) {
	typedef LIBSAKURA_PREFIX::FFT::Type8 T;
	STATIC_ASSERT(sizeof(*src) == sizeof(T) && sizeof(*dst) == sizeof(T));
	return FlipMatrix<T>(&LIBSAKURA_PREFIX::FFT::Flip8, false,
			innerMostUntouched, dims, elements,
			reinterpret_cast<T const *>(src), reinterpret_cast<T *>(dst));
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UnflipMatrixDouble)(
bool innerMostUntouched, size_t dims, size_t const elements[],
		double const src[], double dst[]) {
	typedef LIBSAKURA_PREFIX::FFT::Type8 T;
	STATIC_ASSERT(sizeof(*src) == sizeof(T) && sizeof(*dst) == sizeof(T));
	return FlipMatrix<T>(&LIBSAKURA_PREFIX::FFT::Flip8, true,
			innerMostUntouched, dims, elements,
			reinterpret_cast<T const *>(src), reinterpret_cast<T *>(dst));
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(FlipMatrixDouble2)(
bool innerMostUntouched, size_t dims, size_t const elements[],
		double const src[][2], double dst[][2]) {
	typedef LIBSAKURA_PREFIX::FFT::Type16 T;
	STATIC_ASSERT(sizeof(*src) == sizeof(T) && sizeof(*dst) == sizeof(T));
	return FlipMatrix<T>(&LIBSAKURA_PREFIX::FFT::Flip16, false,
			innerMostUntouched, dims, elements,
			reinterpret_cast<T const *>(src), reinterpret_cast<T *>(dst));
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UnflipMatrixDouble2)(
bool innerMostUntouched, size_t dims, size_t const elements[],
		double const src[][2], double dst[][2]) {
	typedef LIBSAKURA_PREFIX::FFT::Type16 T;
	STATIC_ASSERT(sizeof(*src) == sizeof(T) && sizeof(*dst) == sizeof(T));
	return FlipMatrix<T>(&LIBSAKURA_PREFIX::FFT::Flip16, true,
			innerMostUntouched, dims, elements,
			reinterpret_cast<T const *>(src), reinterpret_cast<T *>(dst));
}
