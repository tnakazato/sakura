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

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ComputeStatistics)(
		size_t elements, float const data[], bool const is_valid[],
		LIBSAKURA_SYMBOL(StatisticsResult) *result) {
	CHECK_ARGS(elements <= INT32_MAX);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(is_valid != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(is_valid));
	CHECK_ARGS(result != nullptr);

	try {
		auto stat =
				::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetStatisticsImpl();
		stat->ComputeStatistics(data, is_valid, elements, result);
	} catch (...) {
		assert(false); // no exception should not be raised for the current implementation.
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

namespace {

template<typename T, typename COMPARATOR>
void QuickSort(T data[], size_t elements) {
	// TODO implement quick sort using expression template to avoid overhead of calling COMPARATOR.
	qsort(data, elements, sizeof(T),
			reinterpret_cast<int (*)(void const*,
					void const*)>(COMPARATOR::Compare));}

template<typename T>
class AscendingOrder {
public:
	static int Compare(T const*a, T const*b) {
		if (false) {
			auto tmp = *a - *b;
			int sign = static_cast<int>(std::abs(tmp) / tmp);
			if (static_cast<int>(0. / 0.) == INT_MIN) {
				sign <<= 1;
			} else if (static_cast<int>(0. / 0.) == 0) {
			} else {
				assert(false);
			}
			return sign;
		} else {
			if (*a < *b) {
				return -1;
			} else if (*a > *b) {
				return 1;
			} else {
				return 0;
			}
		}
	}
};

}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SortValidValuesDensely)(
		size_t elements, bool const is_valid[], float data[],
		size_t *new_elements) {
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(is_valid != nullptr);
	CHECK_ARGS(new_elements != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(is_valid));

	try {
		size_t valid_count = 0;
		for (size_t i = 0; i < elements; ++i) {
			if (is_valid[i]) {
				assert(!isnanf(data[i]));
				data[valid_count] = data[i];
				++valid_count;
			}
		}

		QuickSort<float, AscendingOrder<float> >(data, valid_count);
		*new_elements = valid_count;
	} catch (...) {
		assert(false); // no exception should not be raised for the current implementation.
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}
