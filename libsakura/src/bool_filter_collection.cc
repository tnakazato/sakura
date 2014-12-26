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
#include <iostream>
#include <cmath>

#include "libsakura/sakura.h"
#include "libsakura/localdef.h"

// Vectorization by Compiler
namespace {

template<typename DataType, size_t NUM_BOUNDS>
inline void SetTrueInRangesInclusiveVector(size_t num_data,
		DataType const *data, DataType const *lower_bounds,
		DataType const *upper_bounds,
		bool *result) {
	uint8_t *result_alias = reinterpret_cast<uint8_t *>(result);
	DataType const zero(static_cast<DataType>(0));

	for (size_t i = 0; i < num_data; ++i) {
		uint8_t is_in_range = 0;
		for (size_t j = 0; j < NUM_BOUNDS; ++j) {
			is_in_range |= static_cast<uint8_t>((data[i] - lower_bounds[j])
					* (upper_bounds[j] - data[i]) >= zero);
		}
		result_alias[i] = is_in_range;
	}
}

template<typename DataType, size_t NUM_BOUNDS>
inline void SetTrueInRangesInclusiveScalar(size_t num_data,
		DataType const *data, DataType const *lower_bounds,
		DataType const *upper_bounds,
		bool *result) {
	DataType const zero(static_cast<DataType>(0));

	for (size_t i = 0; i < num_data; ++i) {
		bool is_in_range = false;
		for (size_t j = 0; j < NUM_BOUNDS; ++j) {
			if (((data[i] - lower_bounds[j]) * (upper_bounds[j] - data[i])
					>= zero)) {
				is_in_range = true;
				break;
			}
		}
		result[i] = is_in_range;
	}
}

template<typename DataType>
inline void SetTrueInRangesInclusiveGeneric(size_t num_data,
		DataType const *data, size_t num_condition,
		DataType const *lower_bounds, DataType const *upper_bounds,
		bool *result) {
	DataType const zero(static_cast<DataType>(0));
	for (size_t i = 0; i < num_data; ++i) {
		bool is_in_range = false;
		for (size_t j = 0; j < num_condition; ++j) {
			DataType lower_value = lower_bounds[j];
			DataType upper_value = upper_bounds[j];
			if (((data[i] - lower_value) * (upper_value - data[i]) >= zero)) {
				is_in_range = true;
				break;
			}
		}
		result[i] = is_in_range;
	}
}

template<typename DataType, size_t NUM_BOUNDS>
inline void SetTrueInRangesExclusiveScalar(size_t num_data,
		DataType const *data, DataType const *lower_bounds,
		DataType const *upper_bounds,
		bool *result) {
	DataType const zero(static_cast<DataType>(0));

	for (size_t i = 0; i < num_data; ++i) {
		bool is_in_range = false;
		for (size_t j = 0; j < NUM_BOUNDS; ++j) {
			if (((data[i] - lower_bounds[j]) * (upper_bounds[j] - data[i])
					> zero)) {
				is_in_range = true;
				break;
			}
		}
		result[i] = is_in_range;
	}
}

template<typename DataType>
inline void SetTrueInRangesExclusiveGeneric(size_t num_data,
		DataType const *data, size_t num_condition,
		DataType const *lower_bounds, DataType const *upper_bounds,
		bool *result) {
	DataType const zero(static_cast<DataType>(0));
	for (size_t i = 0; i < num_data; ++i) {
		bool is_in_range = false;
		for (size_t j = 0; j < num_condition; ++j) {
			DataType lower_value = lower_bounds[j];
			DataType upper_value = upper_bounds[j];
			if (((data[i] - lower_value) * (upper_value - data[i]) > zero)) {
				is_in_range = true;
				break;
			}
		}
		result[i] = is_in_range;
	}
}

template<typename DataType>
inline void SetTrueGreaterThan(size_t num_data, DataType const *data,
		DataType threshold, bool *result) {
	DataType const zero(static_cast<DataType>(0));
	auto adata = AssumeAligned(data);
	auto aresult = AssumeAligned(result);
	for (size_t i = 0; i < num_data; ++i) {
		aresult[i] = (adata[i] - threshold) > zero;
	}
}
//template<typename DataType>
//inline void SetTrueGreaterThan(size_t num_data, DataType const *data,
//		DataType threshold, bool *result) {
//	DataType const zero(static_cast<DataType>(0));
//	for (size_t i = 0; i < num_data; ++i) {
//		result[i] = (data[i] - threshold) > zero;
//	}
//}

template<typename DataType>
inline void SetTrueGreaterThanOrEquals(size_t num_data, DataType const *data,
		DataType threshold, bool *result) {
	DataType const zero(static_cast<DataType>(0));
	for (size_t i = 0; i < num_data; ++i) {
		result[i] = (data[i] - threshold) >= zero;
	}
}

template<typename DataType>
inline void SetTrueLessThan(size_t num_data, DataType const *data,
		DataType threshold, bool *result) {
	DataType const zero(static_cast<DataType>(0));
	for (size_t i = 0; i < num_data; ++i) {
		result[i] = (threshold - data[i]) > zero;
	}
}

template<typename DataType>
inline void SetTrueLessThanOrEquals(size_t num_data, DataType const *data,
		DataType threshold, bool *result) {
	DataType const zero(static_cast<DataType>(0));
	for (size_t i = 0; i < num_data; ++i) {
		result[i] = (threshold - data[i]) >= zero;
	}
}

template<typename DataType>
inline void SetFalseIfNanOrInf(size_t num_data, DataType const *data,
		bool *result) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));

	// No operation is done when num_data==0.
	for (size_t i = 0; i < num_data; ++i) {
		result[i] = std::isfinite(data[i]);
	}
}

template<typename DataType>
inline void ToBool(size_t num_data, DataType const *data, bool *result) {
//	std::cout << "Invoking ToBoolDefault()" << std::endl;
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	DataType const zero(static_cast<DataType>(0));
	for (size_t i = 0; i < num_data; ++i) {
		result[i] = (data[i] != zero);
	}
}

inline void InvertBool(size_t num_data, bool const *data, bool *result) {
//	std::cout << "Invoking InvertBoolDefault()" << std::endl;
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));

	uint8_t const *data8 = reinterpret_cast<uint8_t const *>(data);
	uint8_t *result8 = reinterpret_cast<uint8_t *>(result);
	uint8_t true8(static_cast<uint8_t>(true));
	STATIC_ASSERT(sizeof(data8[0]) == sizeof(data[0]));
	STATIC_ASSERT(sizeof(result8[0]) == sizeof(result[0]));
	STATIC_ASSERT(sizeof(data[0]) == sizeof(true8));
	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);
	// No operation is done when num_data==0.
	for (size_t i = 0; i < num_data; ++i) {
		result8[i] = (data8[i] ^ true8);
	}
}

template<typename DataType>
void SetTrueInRangesInclusive(
		size_t num_data, DataType const data[/*num_data*/],
		size_t num_condition, DataType const lower_bounds[/*num_condition*/],
		DataType const upper_bounds[/*num_condition*/],
		bool result[/*num_data*/]) {
	typedef void (*SetTrueInRangesInclusiveFunc)(size_t num_data,
			DataType const *data, DataType const *lower_bounds,
			DataType const *upper_bounds,
			bool *result);
	// Use Scalar version for now
	static SetTrueInRangesInclusiveFunc const funcs[] = {
			SetTrueInRangesInclusiveScalar<DataType, 0>,
			SetTrueInRangesInclusiveScalar<DataType, 1>,
			SetTrueInRangesInclusiveScalar<DataType, 2>,
			SetTrueInRangesInclusiveScalar<DataType, 3>,
			SetTrueInRangesInclusiveScalar<DataType, 4>,
			SetTrueInRangesInclusiveScalar<DataType, 5>,
			SetTrueInRangesInclusiveScalar<DataType, 6>,
			SetTrueInRangesInclusiveScalar<DataType, 7>,
			SetTrueInRangesInclusiveScalar<DataType, 8>,
			SetTrueInRangesInclusiveScalar<DataType, 9>,
			SetTrueInRangesInclusiveScalar<DataType, 10>,
			SetTrueInRangesInclusiveScalar<DataType, 11>,
			SetTrueInRangesInclusiveScalar<DataType, 12>,
			SetTrueInRangesInclusiveScalar<DataType, 13>,
			SetTrueInRangesInclusiveScalar<DataType, 14>,
			SetTrueInRangesInclusiveScalar<DataType, 15>,
			SetTrueInRangesInclusiveScalar<DataType, 16> };

	// So far, only unit8_t version is vectorized
	//std::cout << "Invoking SetTrueInRangesInclusiveDefault()" << std::endl;
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	assert(LIBSAKURA_SYMBOL(IsAligned)(upper_bounds));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lower_bounds));
#ifndef NDEBUG
	for (size_t j = 0; j < num_condition; ++j) {
		assert(lower_bounds[j] <= upper_bounds[j]);
	}
#endif
	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);
	// Initialize result with false

	if (num_condition < ELEMENTSOF(funcs)) {
		funcs[num_condition](num_data, data, lower_bounds, upper_bounds,
				result);
	} else {
		SetTrueInRangesInclusiveGeneric(num_data, data, num_condition,
				lower_bounds, upper_bounds, result);
	}
}

template<typename DataType>
void SetTrueInRangesExclusive(
		size_t num_data, DataType const data[/*num_data*/],
		size_t num_condition, DataType const lower_bounds[/*num_condition*/],
		DataType const upper_bounds[/*num_condition*/],
		bool result[/*num_data*/]) {
	typedef void (*SetTrueInRangesExclusiveFunc)(size_t num_data,
			DataType const *data, DataType const *lower_bounds,
			DataType const *upper_bounds,
			bool *result);
	// Use Scalar version for now
	static SetTrueInRangesExclusiveFunc const funcs[] = {
			SetTrueInRangesExclusiveScalar<DataType, 0>,
			SetTrueInRangesExclusiveScalar<DataType, 1>,
			SetTrueInRangesExclusiveScalar<DataType, 2>,
			SetTrueInRangesExclusiveScalar<DataType, 3>,
			SetTrueInRangesExclusiveScalar<DataType, 4>,
			SetTrueInRangesExclusiveScalar<DataType, 5>,
			SetTrueInRangesExclusiveScalar<DataType, 6>,
			SetTrueInRangesExclusiveScalar<DataType, 7>,
			SetTrueInRangesExclusiveScalar<DataType, 8>,
			SetTrueInRangesExclusiveScalar<DataType, 9>,
			SetTrueInRangesExclusiveScalar<DataType, 10>,
			SetTrueInRangesExclusiveScalar<DataType, 11>,
			SetTrueInRangesExclusiveScalar<DataType, 12>,
			SetTrueInRangesExclusiveScalar<DataType, 13>,
			SetTrueInRangesExclusiveScalar<DataType, 14>,
			SetTrueInRangesExclusiveScalar<DataType, 15>,
			SetTrueInRangesExclusiveScalar<DataType, 16> };

	// So far, only unit8_t version is vectorized
	//std::cout << "Invoking SetTrueInRangesInclusiveDefault()" << std::endl;
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	assert(LIBSAKURA_SYMBOL(IsAligned)(upper_bounds));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lower_bounds));
#ifndef NDEBUG
	for (size_t j = 0; j < num_condition; ++j) {
		assert(lower_bounds[j] <= upper_bounds[j]);
	}
#endif
	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);
	// Initialize result with false

	if (num_condition < ELEMENTSOF(funcs)) {
		funcs[num_condition](num_data, data, lower_bounds, upper_bounds,
				result);
	} else {
		SetTrueInRangesExclusiveGeneric(num_data, data, num_condition,
				lower_bounds, upper_bounds, result);
	}
}

} /* anonymous namespace */

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfInRangesInclusiveFloat)(
		size_t num_data, float const data[], size_t num_condition,
		float const lower_bounds[], float const upper_bounds[], bool result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (lower_bounds == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (upper_bounds == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(lower_bounds)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(upper_bounds)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	// lower_bounds should be smaller or equals to corresponding upper_bounds.
	for (size_t i = 0; i < num_condition; ++i) {
		if (lower_bounds[i] > upper_bounds[i])
			return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}

	// Now actual operation
	try {
		SetTrueInRangesInclusive(num_data, data, num_condition,
				lower_bounds, upper_bounds, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfInRangesInclusiveInt)(
		size_t num_data, int const data[], size_t num_condition,
		int const lower_bounds[], int const upper_bounds[], bool result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (lower_bounds == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (upper_bounds == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(lower_bounds)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(upper_bounds)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	// lower_bounds should be smaller or equals to corresponding upper_bounds.
	for (size_t i = 0; i < num_condition; ++i) {
		if (lower_bounds[i] > upper_bounds[i])
			return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}

	// Now actual operation
	try {
		SetTrueInRangesInclusive(num_data, data, num_condition,
				lower_bounds, upper_bounds, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfInRangesExclusiveFloat)(
		size_t num_data, float const data[], size_t num_condition,
		float const lower_bounds[], float const upper_bounds[], bool result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (lower_bounds == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (upper_bounds == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(lower_bounds)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(upper_bounds)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	// lower_bounds should be smaller or equals to corresponding upper_bounds.
	for (size_t i = 0; i < num_condition; ++i) {
		if (lower_bounds[i] > upper_bounds[i])
			return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}

	// Now actual operation
	try {
		SetTrueInRangesExclusive(num_data, data, num_condition,
				lower_bounds, upper_bounds, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfInRangesExclusiveInt)(
		size_t num_data, int const data[], size_t num_condition,
		int const lower_bounds[], int const upper_bounds[], bool result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (lower_bounds == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (upper_bounds == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(lower_bounds)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(upper_bounds)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	// lower_bounds should be smaller or equals to corresponding upper_bounds.
	for (size_t i = 0; i < num_condition; ++i) {
		if (lower_bounds[i] > upper_bounds[i])
			return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}

	// Now actual operation
	try {
		SetTrueInRangesExclusive(num_data, data, num_condition,
				lower_bounds, upper_bounds, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanFloat)(
		size_t num_data, float const data[], float threshold,
		bool result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	try {
		SetTrueGreaterThan(num_data, data, threshold, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	try {
		SetTrueGreaterThan(num_data, data, threshold, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanOrEqualsFloat)(
		size_t num_data, float const data[], float threshold,
		bool result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	try {
		SetTrueGreaterThanOrEquals(num_data, data, threshold, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanOrEqualsInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	try {
		SetTrueGreaterThanOrEquals(num_data, data, threshold, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanFloat)(
		size_t num_data, float const data[], float threshold,
		bool result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	try {
		SetTrueLessThan(num_data, data, threshold, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	try {
		SetTrueLessThan(num_data, data, threshold, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanOrEqualsFloat)(
		size_t num_data, float const data[], float threshold,
		bool result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	try {
		SetTrueLessThanOrEquals(num_data, data, threshold, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanOrEqualsInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	try {
		SetTrueLessThanOrEquals(num_data, data, threshold, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetFalseIfNanOrInfFloat)(
		size_t num_data,
		float const data[], bool result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	try {
		SetFalseIfNanOrInf(num_data, data, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint8ToBool)(
		size_t num_data, uint8_t const data[], bool result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	try {
		ToBool(num_data, data, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint32ToBool)(
		size_t num_data, uint32_t const data[], bool result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	try {
		ToBool(num_data, data, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InvertBool)(
		size_t num_data,
		bool const data[], bool result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	try {
		InvertBool(num_data, data, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}
