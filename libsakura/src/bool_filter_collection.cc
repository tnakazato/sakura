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
	} else if (!(LIBSAKURA_SYMBOL(IsAligned)(data_array))) {
		return false;
	}
	return true;
}

/* Test data and result arrays*/
template<typename DataType>
bool IsValidDataAndResult(DataType const data[], bool const result[]) {
	if (!IsValidArray(data) || !IsValidArray(result)) {
		return false;
	}
	return true;
}

/* Test range parameters */
template<typename DataType>
bool IsValidBounds(size_t num_condition,
		DataType const lower_bounds[/*num_condition*/],
		DataType const upper_bounds[/*num_condition*/]) {
	if (!IsValidArray(lower_bounds) || !IsValidArray(upper_bounds)) {
		return false;
	}
#ifndef NDEBUG
	// lower_bounds should be smaller or equals to corresponding upper_bounds.
	for (size_t i = 0; i < num_condition; ++i) {
		if (lower_bounds[i] > upper_bounds[i]) {
			return false;
		}
	}
#endif
	return true;
}

/* Invoke Bool Filter of different interface types */
/*
 *  @brief Generate a boolean filter based on per array equation function @a func .
 *
 *  This function calls a function, @a func , with arrays, @a data ,
 *  @a lower_bounds , and @a upper_bounds .
 *
 *  @param[in] func A function to return a boolean array, @a result, based on arrays @a data ,
 *  @a lower_bounds , and @a upper_bounds . It should output a boolean @a result array.
 */
template<typename Func, typename DataType>
LIBSAKURA_SYMBOL(Status) DoRangesBoolFilter(Func func, size_t num_data,
		DataType const data[/*num_data*/], size_t num_condition,
		DataType const lower_bounds[/*num_condition*/],
		DataType const upper_bounds[/*num_condition*/],
		bool result[/*num_data*/]) {
	// Check parameter arguments.
	if (!IsValidDataAndResult(data, result)
			|| !IsValidBounds(num_condition, lower_bounds, upper_bounds)) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
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

/*
 *  @brief Generate a boolean filter based per element equation function @a func .
 *
 *  This function loops over elements of an array, @a data, and call a function,
 *  @a func , with each @a data element.
 *
 *  @param[in] func A function to return a boolean value for a @a data element.
 *  It should take a @a data element and return a boolean.
 */
template<typename Func, typename DataType>
LIBSAKURA_SYMBOL(Status) DoElementFuncBoolFilter(Func func, size_t num_data,
		DataType const data[/*num_data*/], bool result[/*num_data*/]) {
	// Check parameter arguments.
	if (!IsValidDataAndResult(data, result)) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}

	// Now actual operation
	try {
		auto adata = AssumeAligned(data);
		auto aresult = AssumeAligned(result);
		// No operation is done when num_data==0.
		for (size_t i = 0; i < num_data; ++i) {
			aresult[i] = func(adata[i]);
		}
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

/*
 *  @brief Generate a boolean filter based on per array equation function @a func .
 *
 *  This function calls a function, @a func , with a @a data array.
 *
 *  @param[in] func A function to return a boolean array, @a result, based on a @a data array.
 *  It should take a @a data array and output a @result array.
 */
template<typename Func, typename DataType>
LIBSAKURA_SYMBOL(Status) DoArrayFuncBoolFilter(Func func, size_t num_data,
		DataType const data[/*num_data*/], bool const result[/*num_data*/]) {

	// Check parameter arguments.
	if (!IsValidDataAndResult(data, result)) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
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

}
/* anonymous namespace */

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
	//Invoke result[i] = (data[i] > threshold)
	return DoElementFuncBoolFilter(
			[threshold](decltype(data[0]) data_value) {
		return (data_value > threshold);
	}, num_data, data, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) noexcept {
	//Invoke result[i] = (data[i] > threshold)
	return DoElementFuncBoolFilter([threshold](decltype(data[0]) data_value) {
		return (data_value > threshold);
	},	num_data, data, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanOrEqualsFloat)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]) noexcept {
	//Invoke result[i] = (data[i] >= threshold)
	return DoElementFuncBoolFilter(
			[threshold](decltype(data[0]) data_value) {
		return (data_value >= threshold);
	}, num_data, data, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanOrEqualsInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) noexcept {
	//Invoke result[i] = (data[i] >= threshold)
	return DoElementFuncBoolFilter(
			[threshold](decltype(data[0]) data_value) {
		return (data_value >= threshold);
	}, num_data, data, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanFloat)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]) noexcept {
	//Invoke result[i] = (data[i] < threshold)
	return DoElementFuncBoolFilter(
			[threshold](decltype(data[0]) data_value) {
		return (data_value < threshold);
	}, num_data, data, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) noexcept {
	//Invoke result[i] = (data[i] < threshold)
	return DoElementFuncBoolFilter([threshold](decltype(data[0]) data_value) {
		return (data_value < threshold);
	}, num_data, data, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanOrEqualsFloat)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]) noexcept {
	//Invoke result[i] = (data[i] <= threshold)
	return DoElementFuncBoolFilter(
			[threshold](decltype(data[0]) data_value) {
		return (data_value <= threshold);
	}, num_data, data, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanOrEqualsInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) noexcept {
	//Invoke result[i] = (data[i] <= threshold)
	return DoElementFuncBoolFilter(
			[threshold](decltype(data[0]) data_value) {
		return (data_value <= threshold);
	}, num_data, data, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetFalseIfNanOrInfFloat)(
		size_t num_data, float const data[/*num_data*/],
		bool result[/*num_data*/]) noexcept {
	//Invoke result[i] = std::isfinite(data[i])
	return DoElementFuncBoolFilter([](decltype(data[0]) data_value) {
		return std::isfinite(data_value);
	}, num_data, data, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint8ToBool)(
		size_t num_data, uint8_t const data[/*num_data*/],
		bool result[/*num_data*/]) noexcept {
	constexpr uint8_t kZero(0);
	//Invoke result[i] = (data[i] != 0)
	return DoElementFuncBoolFilter([kZero](decltype(data[0]) data_value) {
		return (data_value != kZero);
	},	num_data, data, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint32ToBool)(
		size_t num_data, uint32_t const data[/*num_data*/],
		bool result[/*num_data*/]) noexcept {
	constexpr uint8_t kZero(0);
	//Invoke result[i] = (data[i] != 0)
	return DoElementFuncBoolFilter([kZero](decltype(data[0]) data_value) {
		return (data_value != kZero);
	},	num_data, data, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InvertBool)(
		size_t num_data,
		bool const data[/*num_data*/], bool result[/*num_data*/]) noexcept {
	return DoArrayFuncBoolFilter([=] {InvertBool(num_data, data, result);},
			num_data, data, result);
}
