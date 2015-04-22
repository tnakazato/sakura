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
#include <cmath>
#include <iostream>

#include "libsakura/localdef.h"
#include "libsakura/sakura.h"

// Vectorization by Compiler
namespace {

template<typename DataType, size_t kNumBounds>
inline void SetTrueIfInRangesInclusiveVector(size_t num_data,
		DataType const *data, DataType const *lower_bounds,
		DataType const *upper_bounds,
		bool *result) {
	uint8_t *result_alias = reinterpret_cast<uint8_t *>(result);
	constexpr DataType kZero(0);

	for (size_t i = 0; i < num_data; ++i) {
		uint8_t is_in_range = 0;
		for (size_t j = 0; j < kNumBounds; ++j) {
			is_in_range |= static_cast<uint8_t>((data[i] - lower_bounds[j])
					* (upper_bounds[j] - data[i]) >= kZero);
		}
		result_alias[i] = is_in_range;
	}
}

template<typename DataType, size_t kNumBounds>
inline void SetTrueIfInRangesInclusiveScalar(size_t num_data,
		DataType const *data, DataType const *lower_bounds,
		DataType const *upper_bounds,
		bool *result) {
	constexpr DataType kZero(0);

	for (size_t i = 0; i < num_data; ++i) {
		bool is_in_range = false;
		for (size_t j = 0; j < kNumBounds; ++j) {
			if (((data[i] - lower_bounds[j]) * (upper_bounds[j] - data[i])
					>= kZero)) {
				is_in_range = true;
				break;
			}
		}
		result[i] = is_in_range;
	}
}

template<typename DataType>
inline void SetTrueIfInRangesInclusiveGeneric(size_t num_data,
		DataType const *data, size_t num_condition,
		DataType const *lower_bounds, DataType const *upper_bounds,
		bool *result) {
	constexpr DataType kZero(0);
	for (size_t i = 0; i < num_data; ++i) {
		bool is_in_range = false;
		for (size_t j = 0; j < num_condition; ++j) {
			DataType lower_value = lower_bounds[j];
			DataType upper_value = upper_bounds[j];
			if (((data[i] - lower_value) * (upper_value - data[i]) >= kZero)) {
				is_in_range = true;
				break;
			}
		}
		result[i] = is_in_range;
	}
}

template<typename DataType, size_t kNumBounds>
inline void SetTrueIfInRangesExclusiveScalar(size_t num_data,
		DataType const *data, DataType const *lower_bounds,
		DataType const *upper_bounds,
		bool *result) {
	constexpr DataType kZero(0);

	for (size_t i = 0; i < num_data; ++i) {
		bool is_in_range = false;
		for (size_t j = 0; j < kNumBounds; ++j) {
			if (((data[i] - lower_bounds[j]) * (upper_bounds[j] - data[i])
					> kZero)) {
				is_in_range = true;
				break;
			}
		}
		result[i] = is_in_range;
	}
}

template<typename DataType>
inline void SetTrueIfInRangesExclusiveGeneric(size_t num_data,
		DataType const *data, size_t num_condition,
		DataType const *lower_bounds, DataType const *upper_bounds,
		bool *result) {
	constexpr DataType kZero(0);
	for (size_t i = 0; i < num_data; ++i) {
		bool is_in_range = false;
		for (size_t j = 0; j < num_condition; ++j) {
			DataType lower_value = lower_bounds[j];
			DataType upper_value = upper_bounds[j];
			if (((data[i] - lower_value) * (upper_value - data[i]) > kZero)) {
				is_in_range = true;
				break;
			}
		}
		result[i] = is_in_range;
	}
}

template<typename DataType>
inline void SetTrueIfGreaterThan(size_t num_data, DataType const *data,
		DataType threshold, bool *result) {
	auto adata = AssumeAligned(data);
	auto aresult = AssumeAligned(result);
	for (size_t i = 0; i < num_data; ++i) {
		aresult[i] = (adata[i] > threshold);
	}
}
//template<typename DataType>
//inline void SetTrueIfGreaterThan(size_t num_data, DataType const *data,
//		DataType threshold, bool *result) {
//	for (size_t i = 0; i < num_data; ++i) {
//		result[i] = (data[i] > threshold);
//	}
//}

template<typename DataType>
inline void SetTrueIfGreaterThanOrEquals(size_t num_data, DataType const *data,
		DataType threshold, bool *result) {
	for (size_t i = 0; i < num_data; ++i) {
		result[i] = (data[i] >= threshold);
	}
}

template<typename DataType>
inline void SetTrueIfLessThan(size_t num_data, DataType const *data,
		DataType threshold, bool *result) {
	for (size_t i = 0; i < num_data; ++i) {
		result[i] = (threshold > data[i]);
	}
}

template<typename DataType>
inline void SetTrueIfLessThanOrEquals(size_t num_data, DataType const *data,
		DataType threshold, bool *result) {
	for (size_t i = 0; i < num_data; ++i) {
		result[i] = (threshold >= data[i]);
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
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	constexpr DataType kZero(0);
	for (size_t i = 0; i < num_data; ++i) {
		result[i] = (data[i] != kZero);
	}
}

inline void InvertBool(size_t num_data, bool const *data, bool *result) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));

	uint8_t const *data8 = reinterpret_cast<uint8_t const *>(data);
	uint8_t *result8 = reinterpret_cast<uint8_t *>(result);
	constexpr uint8_t kTrue8(static_cast<uint8_t>(true));
	STATIC_ASSERT(sizeof(data8[0]) == sizeof(data[0]));
	STATIC_ASSERT(sizeof(result8[0]) == sizeof(result[0]));
	STATIC_ASSERT(sizeof(data[0]) == sizeof(kTrue8));
	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);
	// No operation is done when num_data==0.
	for (size_t i = 0; i < num_data; ++i) {
		result8[i] = (data8[i] ^ kTrue8);
	}
}

template<typename DataType>
void SetTrueIfInRangesInclusive(size_t num_data,
		DataType const data[/*num_data*/], size_t num_condition,
		DataType const lower_bounds[/*num_condition*/],
		DataType const upper_bounds[/*num_condition*/],
		bool result[/*num_data*/]) {
	typedef void (*SetTrueIfInRangesInclusiveFunc)(size_t num_data,
			DataType const *data, DataType const *lower_bounds,
			DataType const *upper_bounds,
			bool *result);
	// Use Scalar version for now
	static SetTrueIfInRangesInclusiveFunc const funcs[] = {
			SetTrueIfInRangesInclusiveScalar<DataType, 0>,
			SetTrueIfInRangesInclusiveScalar<DataType, 1>,
			SetTrueIfInRangesInclusiveScalar<DataType, 2>,
			SetTrueIfInRangesInclusiveScalar<DataType, 3>,
			SetTrueIfInRangesInclusiveScalar<DataType, 4>,
			SetTrueIfInRangesInclusiveScalar<DataType, 5>,
			SetTrueIfInRangesInclusiveScalar<DataType, 6>,
			SetTrueIfInRangesInclusiveScalar<DataType, 7>,
			SetTrueIfInRangesInclusiveScalar<DataType, 8>,
			SetTrueIfInRangesInclusiveScalar<DataType, 9>,
			SetTrueIfInRangesInclusiveScalar<DataType, 10>,
			SetTrueIfInRangesInclusiveScalar<DataType, 11>,
			SetTrueIfInRangesInclusiveScalar<DataType, 12>,
			SetTrueIfInRangesInclusiveScalar<DataType, 13>,
			SetTrueIfInRangesInclusiveScalar<DataType, 14>,
			SetTrueIfInRangesInclusiveScalar<DataType, 15>,
			SetTrueIfInRangesInclusiveScalar<DataType, 16> };

	// So far, only unit8_t version is vectorized
	//std::cout << "Invoking SetTrueIfInRangesInclusiveDefault()" << std::endl;
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
		SetTrueIfInRangesInclusiveGeneric(num_data, data, num_condition,
				lower_bounds, upper_bounds, result);
	}
}

template<typename DataType>
void SetTrueIfInRangesExclusive(size_t num_data,
		DataType const data[/*num_data*/], size_t num_condition,
		DataType const lower_bounds[/*num_condition*/],
		DataType const upper_bounds[/*num_condition*/],
		bool result[/*num_data*/]) {
	typedef void (*SetTrueIfInRangesExclusiveFunc)(size_t num_data,
			DataType const *data, DataType const *lower_bounds,
			DataType const *upper_bounds,
			bool *result);
	// Use Scalar version for now
	static SetTrueIfInRangesExclusiveFunc const funcs[] = {
			SetTrueIfInRangesExclusiveScalar<DataType, 0>,
			SetTrueIfInRangesExclusiveScalar<DataType, 1>,
			SetTrueIfInRangesExclusiveScalar<DataType, 2>,
			SetTrueIfInRangesExclusiveScalar<DataType, 3>,
			SetTrueIfInRangesExclusiveScalar<DataType, 4>,
			SetTrueIfInRangesExclusiveScalar<DataType, 5>,
			SetTrueIfInRangesExclusiveScalar<DataType, 6>,
			SetTrueIfInRangesExclusiveScalar<DataType, 7>,
			SetTrueIfInRangesExclusiveScalar<DataType, 8>,
			SetTrueIfInRangesExclusiveScalar<DataType, 9>,
			SetTrueIfInRangesExclusiveScalar<DataType, 10>,
			SetTrueIfInRangesExclusiveScalar<DataType, 11>,
			SetTrueIfInRangesExclusiveScalar<DataType, 12>,
			SetTrueIfInRangesExclusiveScalar<DataType, 13>,
			SetTrueIfInRangesExclusiveScalar<DataType, 14>,
			SetTrueIfInRangesExclusiveScalar<DataType, 15>,
			SetTrueIfInRangesExclusiveScalar<DataType, 16> };

	// So far, only unit8_t version is vectorized
	//std::cout << "Invoking SetTrueIfInRangesInclusiveDefault()" << std::endl;
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
		SetTrueIfInRangesExclusiveGeneric(num_data, data, num_condition,
				lower_bounds, upper_bounds, result);
	}
}

/* Helper functions to test array parameters*/
/* Returns false if the array is either nullptr or not aligned */
template<typename DataType>
bool IsValidArray(DataType const data_array[]) {
	if (data_array == nullptr) {
		return false;
	} else if (!( LIBSAKURA_SYMBOL(IsAligned)(data_array))) {
		return false;
	}
	return true;
}

/* Test data and result arrays*/
template<typename DataType>
bool CheckArrays(DataType const data[], bool const result[],
LIBSAKURA_SYMBOL(Status) *status) {
	if (!IsValidArray(data) || !IsValidArray(result)) {
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		return false;
	}
	return true;
}

/* Test range parameters */
template<typename DataType>
bool CheckRanges(size_t num_condition,
		DataType const lower_bounds[/*num_condition*/],
		DataType const upper_bounds[/*num_condition*/],
		LIBSAKURA_SYMBOL(Status) *status) {
	if (!IsValidArray(lower_bounds) || !IsValidArray(upper_bounds)) {
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		return false;
	}
#ifndef NDEBUG
	// lower_bounds should be smaller or equals to corresponding upper_bounds.
	for (size_t i = 0; i < num_condition; ++i) {
		if (lower_bounds[i] > upper_bounds[i]) {
			*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
			return false;
		}
	}
#endif
	return true;
}

/* Invoke Bool Filter of each interface type */
template<typename Func, typename DataType>
LIBSAKURA_SYMBOL(Status) DoRangesBoolFilter(Func func, size_t num_data,
		DataType const data[/*num_data*/], size_t num_condition,
		DataType const lower_bounds[/*num_condition*/],
		DataType const upper_bounds[/*num_condition*/],
		bool result[/*num_data*/]) {
	// Check parameter arguments.
	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Status_kOK);
	if (!CheckArrays(data, result, &status)
			|| !CheckRanges(num_condition, lower_bounds, upper_bounds,
					&status)) {
		return status;
	}

	// Now actual operation
	try {
		func(/*num_data, data, num_condition, lower_bounds, upper_bounds, result*/);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

template<typename Func, typename DataType>
LIBSAKURA_SYMBOL(Status) DoBoundaryBoolFilter(Func func, size_t num_data,
		DataType const data[/*num_data*/], DataType threshold,
		bool result[/*num_data*/]) {
	// Check parameter arguments.
	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Status_kOK);
	if (!CheckArrays(data, result, &status)) {
		return status;
	}

	// Now actual operation
	try {
		func(/*num_data, data, threshold, result*/);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);

}

template<typename Func, typename DataType>
LIBSAKURA_SYMBOL(Status) DoSimpleBoolFilter(Func func, size_t num_data,
		DataType const data[/*num_data*/], bool const result[/*num_data*/]) {

	// Check parameter arguments.
	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Status_kOK);
	if (!CheckArrays(data, result, &status)) {
		return status;
	}

	// Now actual operation
	try {
		func(/*num_data, data, result*/);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

} /* anonymous namespace */

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfInRangesInclusiveFloat)(
		size_t num_data, float const data[/*num_data*/], size_t num_condition,
		float const lower_bounds[/*num_condition*/],
		float const upper_bounds[/*num_condition*/], bool result[/*num_data*/])
				noexcept {
	return DoRangesBoolFilter(
			[=] {SetTrueIfInRangesInclusive(num_data, data, num_condition, lower_bounds, upper_bounds, result);},
			num_data, data, num_condition, lower_bounds, upper_bounds, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfInRangesInclusiveInt)(
		size_t num_data, int const data[/*num_data*/], size_t num_condition,
		int const lower_bounds[/*num_condition*/],
		int const upper_bounds[/*num_condition*/], bool result[/*num_data*/])
				noexcept {
	return DoRangesBoolFilter(
			[=] {SetTrueIfInRangesInclusive(num_data, data, num_condition, lower_bounds, upper_bounds, result);},
			num_data, data, num_condition, lower_bounds, upper_bounds, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfInRangesExclusiveFloat)(
		size_t num_data, float const data[/*num_data*/], size_t num_condition,
		float const lower_bounds[/*num_condition*/],
		float const upper_bounds[/*num_condition*/], bool result[/*num_data*/])
				noexcept {
	return DoRangesBoolFilter(
			[=] {SetTrueIfInRangesExclusive(num_data, data, num_condition, lower_bounds, upper_bounds, result);},
			num_data, data, num_condition, lower_bounds, upper_bounds, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfInRangesExclusiveInt)(
		size_t num_data, int const data[/*num_data*/], size_t num_condition,
		int const lower_bounds[/*num_condition*/],
		int const upper_bounds[/*num_condition*/], bool result[/*num_data*/])
				noexcept {
	return DoRangesBoolFilter(
			[=] {SetTrueIfInRangesExclusive(num_data, data, num_condition, lower_bounds, upper_bounds, result);},
			num_data, data, num_condition, lower_bounds, upper_bounds, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanFloat)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]) noexcept {
	return DoBoundaryBoolFilter(
			[=] {SetTrueIfGreaterThan(num_data, data, threshold, result);},
			num_data, data, threshold, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) noexcept {
	return DoBoundaryBoolFilter(
			[=] {SetTrueIfGreaterThan(num_data, data, threshold, result);},
			num_data, data, threshold, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanOrEqualsFloat)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]) noexcept {
	return DoBoundaryBoolFilter(
			[=] {SetTrueIfGreaterThanOrEquals(num_data, data, threshold, result);},
			num_data, data, threshold, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanOrEqualsInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) noexcept {
	return DoBoundaryBoolFilter(
			[=] {SetTrueIfGreaterThanOrEquals(num_data, data, threshold, result);},
			num_data, data, threshold, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanFloat)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]) noexcept {
	return DoBoundaryBoolFilter(
			[=] {SetTrueIfLessThan(num_data, data, threshold, result);},
			num_data, data, threshold, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) noexcept {
	return DoBoundaryBoolFilter(
			[=] {SetTrueIfLessThan(num_data, data, threshold, result);},
			num_data, data, threshold, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanOrEqualsFloat)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]) noexcept {
	return DoBoundaryBoolFilter(
			[=] {SetTrueIfLessThanOrEquals(num_data, data, threshold, result);},
			num_data, data, threshold, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanOrEqualsInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) noexcept {
	return DoBoundaryBoolFilter(
			[=] {SetTrueIfLessThanOrEquals(num_data, data, threshold, result);},
			num_data, data, threshold, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetFalseIfNanOrInfFloat)(
		size_t num_data, float const data[/*num_data*/],
		bool result[/*num_data*/]) noexcept {
	return DoSimpleBoolFilter([=] {SetFalseIfNanOrInf(num_data, data, result);},
			num_data, data, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint8ToBool)(
		size_t num_data, uint8_t const data[/*num_data*/],
		bool result[/*num_data*/]) noexcept {
	return DoSimpleBoolFilter([=] {ToBool(num_data, data, result);}, num_data,
			data, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint32ToBool)(
		size_t num_data, uint32_t const data[/*num_data*/],
		bool result[/*num_data*/]) noexcept {
	return DoSimpleBoolFilter([=] {ToBool(num_data, data, result);}, num_data,
			data, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InvertBool)(
		size_t num_data,
		bool const data[/*num_data*/], bool result[/*num_data*/]) noexcept {
	return DoSimpleBoolFilter([=] {InvertBool(num_data, data, result);},
			num_data, data, result);
}
