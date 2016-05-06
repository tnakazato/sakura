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
#include <cmath>
#include <iostream>
#include <string>

#include "libsakura/localdef.h"
#include "libsakura/sakura.h"
#include "libsakura/packed_operation.h"

// Vectorization by Compiler
namespace {

#define CHECK_ARGS(x) do { \
	if (!(x)) { \
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument); \
	} \
} while (false)

template<typename DataType, size_t kNumBounds>
struct SetTrueIfInRangesInclusiveVector {
	inline static void process(size_t num_data, DataType const *data,
			DataType const *lower_bounds, DataType const *upper_bounds,
			bool *result) {
		constexpr DataType kZero = 0;

		for (size_t i = 0; i < num_data; ++i) {
			bool is_in_range = false;
			for (size_t j = 0; j < kNumBounds; ++j) {
				is_in_range |= ((data[i] - lower_bounds[j])
						* (upper_bounds[j] - data[i]) >= kZero);
			}
			result[i] = is_in_range;
		}
	}
};

#if defined(__AVX__)
template<size_t kNumBounds>
struct SetTrueIfInRangesInclusiveVector<float, kNumBounds> {
	inline static void process(size_t num_data, float const *data,
			float const *lower_bounds, float const *upper_bounds,
			bool *result) {
		constexpr float kZero = 0.0f;
		constexpr auto kElementsPerLoop = LIBSAKURA_SYMBOL(SimdPacketAVX)::kSize
				/ sizeof(data[0]);
		LIBSAKURA_SYMBOL(SimdPacketAVX) upper, lower, data_packet;
		//pack by data
		const auto data_ptr =
				AssumeAligned(
						reinterpret_cast<LIBSAKURA_SYMBOL(SimdPacketAVX)::RawFloat const *>(data));
		const auto result_ptr = AssumeAligned(
				reinterpret_cast<uint64_t *>(result));
		const auto n = num_data / kElementsPerLoop;
		STATIC_ASSERT(false==0);
		const auto truth = _mm_set1_epi8(true);
		LIBSAKURA_SYMBOL(SimdPacketAVX) is_in_range;
		for (size_t i = 0; i < n; ++i) {
			data_packet.raw_float = data_ptr[i];
			is_in_range.set1(static_cast<float>(0)); // false
			for (size_t j = 0; j < kNumBounds; ++j) {
				lower.set1(lower_bounds[j]);
				upper.set1(upper_bounds[j]);
				//Returns 0xFFFFFFFF if data is in one of ranges, else 0x00000000
				is_in_range =
						LIBSAKURA_SYMBOL(SimdMath)<
						LIBSAKURA_SYMBOL(SimdArchAVX), float>::Or(is_in_range,
								LIBSAKURA_SYMBOL(SimdMath)<
										LIBSAKURA_SYMBOL(SimdArchAVX), float>::And(
										LIBSAKURA_SYMBOL(SimdCompare)<
												LIBSAKURA_SYMBOL(SimdArchAVX),
												float>::LessOrEqual(lower,
												data_packet),
										LIBSAKURA_SYMBOL(SimdCompare)<
												LIBSAKURA_SYMBOL(SimdArchAVX),
												float>::LessOrEqual(data_packet,
												upper)));
			}

			LIBSAKURA_SYMBOL(SimdArchAVX)::PriorArch::PacketType result_prior =
			LIBSAKURA_SYMBOL(SimdConvert)<
			LIBSAKURA_SYMBOL(SimdArchAVX)>::Int32ToLSByte(is_in_range);
			result_ptr[i] = _mm_extract_epi64(
					_mm_and_si128(result_prior.raw_int32, truth), 0);
		}

		// process remaining elements
		const auto end = n * kElementsPerLoop;
		data = &data[end];
		result = &result[end];
		num_data -= end;
		uint8_t *result_alias = reinterpret_cast<uint8_t *>(result);

		for (size_t i = 0; i < num_data; ++i) {
			bool is_in_range = false;
			for (size_t j = 0; j < kNumBounds; ++j) {
				is_in_range |= ((data[i] - lower_bounds[j])
						* (upper_bounds[j] - data[i]) >= kZero);
			}
			result_alias[i] = is_in_range;
		}
	}
};

template<size_t kNumBounds>
struct SetTrueIfInRangesInclusiveVector<int, kNumBounds> {
	inline static void process(size_t num_data, int const *data,
			int const *lower_bounds, int const *upper_bounds,
			bool *result) {
		constexpr int kZero = 0;
#if defined(__AVX2__)
		typedef LIBSAKURA_SYMBOL(SimdPacketAVX) MySimdPacket;
		typedef LIBSAKURA_SYMBOL(SimdArchAVX) MySimdArch;
#else
		typedef LIBSAKURA_SYMBOL(SimdPacketSSE) MySimdPacket;
		typedef LIBSAKURA_SYMBOL(SimdArchSSE) MySimdArch;
#endif
		constexpr auto kElementsPerLoop = MySimdPacket::kSize / sizeof(data[0]);
		MySimdPacket upper, lower, data_packet;
		//pack by data
		const auto data_ptr = AssumeAligned(
				reinterpret_cast<MySimdPacket::RawInt32 const *>(data));
#if defined(__AVX2__)
		const auto result_ptr = AssumeAligned(
				reinterpret_cast<uint64_t *>(result));
#else
		const auto result_ptr = AssumeAligned(
				reinterpret_cast<uint32_t *>(result));
#endif
		const auto n = num_data / kElementsPerLoop;
		STATIC_ASSERT(false==0);
		MySimdPacket is_in_range;
#if defined(__AVX2__)
		const auto truth = _mm_set1_epi8(true);
#else
		const auto truth = _mm_set1_pi8(true);
#endif
		for (size_t i = 0; i < n; ++i) {
			data_packet.raw_int32 = data_ptr[i];
			is_in_range.set1(static_cast<int>(0)); // false
			for (size_t j = 0; j < kNumBounds; ++j) {
				lower.set1(static_cast<int32_t>(lower_bounds[j]));
				upper.set1(static_cast<int32_t>(upper_bounds[j]));
				//Returns 0xFFFFFFFF if data is in one of ranges, else 0x00000000
				is_in_range =
				LIBSAKURA_SYMBOL(SimdMath)<MySimdArch, int32_t>::Or(is_in_range,
						LIBSAKURA_SYMBOL(SimdMath)<MySimdArch, int32_t>::And(
								LIBSAKURA_SYMBOL(SimdCompare)<MySimdArch,
										int32_t>::LessOrEqual(lower,
										data_packet),
								LIBSAKURA_SYMBOL(SimdCompare)<MySimdArch,
										int32_t>::LessOrEqual(data_packet,
										upper)));

			}
			MySimdArch::PriorArch::PacketType result_prior =
					LIBSAKURA_SYMBOL(SimdConvert)<MySimdArch>::Int32ToLSByte(
							is_in_range);
#if defined(__AVX2__)
			result_ptr[i] = _mm_extract_epi64(
					_mm_and_si128(result_prior.raw_int32, truth), 0);
#else
			result_ptr[i] = _mm_cvtsi64_si32(
					_mm_and_si64(result_prior.raw_int64, truth));
#endif
		}

		// process remaining elements
		const auto end = n * kElementsPerLoop;
		data = &data[end];
		result = &result[end];
		num_data -= end;
		uint8_t *result_alias = reinterpret_cast<uint8_t *>(result);

		for (size_t i = 0; i < num_data; ++i) {
			bool is_in_range = false;
			for (size_t j = 0; j < kNumBounds; ++j) {
				is_in_range |= ((data[i] - lower_bounds[j])
						* (upper_bounds[j] - data[i]) >= kZero);
			}
			result_alias[i] = is_in_range;
		}
	}
};
#endif

template<typename DataType>
inline void SetTrueIfInRangesInclusiveGeneric(size_t num_data,
		DataType const *data, size_t num_condition,
		DataType const *lower_bounds, DataType const *upper_bounds,
		bool *result) {
	constexpr DataType kZero = 0;
	for (size_t i = 0; i < num_data; ++i) {
		bool is_in_range = false;
		for (size_t j = 0; j < num_condition; ++j) {
			if (((data[i] - lower_bounds[j]) * (upper_bounds[j] - data[i]) >= kZero)) {
				is_in_range = true;
				break;
			}
		}
		result[i] = is_in_range;
	}
}

template<typename DataType, size_t kNumBounds>
struct SetTrueIfInRangesExclusiveVector {
	inline static void process(size_t num_data,
		DataType const *data, DataType const *lower_bounds,
		DataType const *upper_bounds,
		bool *result) {
	constexpr DataType kZero = 0;

	for (size_t i = 0; i < num_data; ++i) {
		bool is_in_range = false;
		for (size_t j = 0; j < kNumBounds; ++j) {
			is_in_range |= ((data[i] - lower_bounds[j])
					* (upper_bounds[j] - data[i]) > kZero);
		}
		result[i] = is_in_range;
	}
}
};

#if defined(__AVX__)
template<size_t kNumBounds>
struct SetTrueIfInRangesExclusiveVector<float, kNumBounds> {
	inline static void process(size_t num_data, float const *data,
			float const *lower_bounds, float const *upper_bounds,
			bool *result) {
		constexpr float kZero = 0.0f;
		constexpr auto kElementsPerLoop = LIBSAKURA_SYMBOL(SimdPacketAVX)::kSize
				/ sizeof(data[0]);
		LIBSAKURA_SYMBOL(SimdPacketAVX) upper, lower, data_packet;
		//pack by data
		const auto data_ptr =
				AssumeAligned(
						reinterpret_cast<LIBSAKURA_SYMBOL(SimdPacketAVX)::RawFloat const *>(data));
		const auto result_ptr = AssumeAligned(
				reinterpret_cast<uint64_t *>(result));
		const auto n = num_data / kElementsPerLoop;
		STATIC_ASSERT(false==0);
		const auto truth = _mm_set1_epi8(true);
		LIBSAKURA_SYMBOL(SimdPacketAVX) is_in_range;
		for (size_t i = 0; i < n; ++i) {
			data_packet.raw_float = data_ptr[i];
			is_in_range.set1(static_cast<float>(0)); // false
			for (size_t j = 0; j < kNumBounds; ++j) {
				lower.set1(lower_bounds[j]);
				upper.set1(upper_bounds[j]);
				//Returns 0xFFFFFFFF if data is in one of ranges, else 0x00000000
				is_in_range =
						LIBSAKURA_SYMBOL(SimdMath)<
						LIBSAKURA_SYMBOL(SimdArchAVX), float>::Or(is_in_range,
								LIBSAKURA_SYMBOL(SimdMath)<
										LIBSAKURA_SYMBOL(SimdArchAVX), float>::And(
										LIBSAKURA_SYMBOL(SimdCompare)<
												LIBSAKURA_SYMBOL(SimdArchAVX),
												float>::LessThan(lower,
												data_packet),
										LIBSAKURA_SYMBOL(SimdCompare)<
												LIBSAKURA_SYMBOL(SimdArchAVX),
												float>::LessThan(data_packet,
												upper)));
			}

			LIBSAKURA_SYMBOL(SimdArchAVX)::PriorArch::PacketType result_prior =
			LIBSAKURA_SYMBOL(SimdConvert)<
			LIBSAKURA_SYMBOL(SimdArchAVX)>::Int32ToLSByte(is_in_range);
			result_ptr[i] = _mm_extract_epi64(
					_mm_and_si128(result_prior.raw_int32, truth), 0);
		}

		// process remaining elements
		const auto end = n * kElementsPerLoop;
		data = &data[end];
		result = &result[end];
		num_data -= end;
		uint8_t *result_alias = reinterpret_cast<uint8_t *>(result);

		for (size_t i = 0; i < num_data; ++i) {
			bool is_in_range = false;
			for (size_t j = 0; j < kNumBounds; ++j) {
				is_in_range |= ((data[i] - lower_bounds[j])
						* (upper_bounds[j] - data[i]) > kZero);
			}
			result_alias[i] = is_in_range;
		}
	}
};
template<size_t kNumBounds>
struct SetTrueIfInRangesExclusiveVector<int, kNumBounds> {
	inline static void process(size_t num_data, int const *data,
			int const *lower_bounds, int const *upper_bounds,
			bool *result) {
		constexpr int kZero = 0;
#if defined(__AVX2__)
		typedef LIBSAKURA_SYMBOL(SimdPacketAVX) MySimdPacket;
		typedef LIBSAKURA_SYMBOL(SimdArchAVX) MySimdArch;
#else
		typedef LIBSAKURA_SYMBOL(SimdPacketSSE) MySimdPacket;
		typedef LIBSAKURA_SYMBOL(SimdArchSSE) MySimdArch;
#endif
		constexpr auto kElementsPerLoop = MySimdPacket::kSize / sizeof(data[0]);
		MySimdPacket upper, lower, data_packet;
		//pack by data
		const auto data_ptr = AssumeAligned(
				reinterpret_cast<MySimdPacket::RawInt32 const *>(data));
#if defined(__AVX2__)
		const auto result_ptr = AssumeAligned(
				reinterpret_cast<uint64_t *>(result));
#else
		const auto result_ptr = AssumeAligned(
				reinterpret_cast<uint32_t *>(result));
#endif
		const auto n = num_data / kElementsPerLoop;
		STATIC_ASSERT(false==0);
		MySimdPacket is_in_range;
#if defined(__AVX2__)
		const auto truth = _mm_set1_epi8(true);
#else
		const auto truth = _mm_set1_pi8(true);
#endif
		for (size_t i = 0; i < n; ++i) {
			data_packet.raw_int32 = data_ptr[i];
			is_in_range.set1(static_cast<int>(0)); // false
			for (size_t j = 0; j < kNumBounds; ++j) {
				lower.set1(static_cast<int32_t>(lower_bounds[j]));
				upper.set1(static_cast<int32_t>(upper_bounds[j]));
				//Returns 0xFFFFFFFF if data is in one of ranges, else 0x00000000
				is_in_range =
				LIBSAKURA_SYMBOL(SimdMath)<MySimdArch, int32_t>::Or(is_in_range,
						LIBSAKURA_SYMBOL(SimdMath)<MySimdArch, int32_t>::And(
								LIBSAKURA_SYMBOL(SimdCompare)<MySimdArch,
										int32_t>::LessThan(lower,
										data_packet),
								LIBSAKURA_SYMBOL(SimdCompare)<MySimdArch,
										int32_t>::LessThan(data_packet,
										upper)));

			}
			MySimdArch::PriorArch::PacketType result_prior =
					LIBSAKURA_SYMBOL(SimdConvert)<MySimdArch>::Int32ToLSByte(
							is_in_range);
#if defined(__AVX2__)
			result_ptr[i] = _mm_extract_epi64(
					_mm_and_si128(result_prior.raw_int32, truth), 0);
#else
			result_ptr[i] = _mm_cvtsi64_si32(
					_mm_and_si64(result_prior.raw_int64, truth));
#endif
		}

		// process remaining elements
		const auto end = n * kElementsPerLoop;
		data = &data[end];
		result = &result[end];
		num_data -= end;
		uint8_t *result_alias = reinterpret_cast<uint8_t *>(result);

		for (size_t i = 0; i < num_data; ++i) {
			bool is_in_range = false;
			for (size_t j = 0; j < kNumBounds; ++j) {
				is_in_range |= ((data[i] - lower_bounds[j])
						* (upper_bounds[j] - data[i]) > kZero);
			}
			result_alias[i] = is_in_range;
		}
	}
};
#endif

template<typename DataType>
inline void SetTrueIfInRangesExclusiveGeneric(size_t num_data,
		DataType const *data, size_t num_condition,
		DataType const *lower_bounds, DataType const *upper_bounds,
		bool *result) {
	constexpr DataType kZero = 0;
	for (size_t i = 0; i < num_data; ++i) {
		bool is_in_range = false;
		for (size_t j = 0; j < num_condition; ++j) {
			if (((data[i] - lower_bounds[j]) * (upper_bounds[j] - data[i]) > kZero)) {
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
	constexpr uint8_t kTrue8 = static_cast<uint8_t>(true);
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
			SetTrueIfInRangesInclusiveVector<DataType, 0>::process,
			SetTrueIfInRangesInclusiveVector<DataType, 1>::process,
			SetTrueIfInRangesInclusiveVector<DataType, 2>::process,
			SetTrueIfInRangesInclusiveVector<DataType, 3>::process,
			SetTrueIfInRangesInclusiveVector<DataType, 4>::process,
			SetTrueIfInRangesInclusiveVector<DataType, 5>::process,
			SetTrueIfInRangesInclusiveVector<DataType, 6>::process,
			SetTrueIfInRangesInclusiveVector<DataType, 7>::process,
			SetTrueIfInRangesInclusiveVector<DataType, 8>::process,
			SetTrueIfInRangesInclusiveVector<DataType, 9>::process,
			SetTrueIfInRangesInclusiveVector<DataType, 10>::process,
			SetTrueIfInRangesInclusiveVector<DataType, 11>::process,
			SetTrueIfInRangesInclusiveVector<DataType, 12>::process,
			SetTrueIfInRangesInclusiveVector<DataType, 13>::process,
			SetTrueIfInRangesInclusiveVector<DataType, 14>::process,
			SetTrueIfInRangesInclusiveVector<DataType, 15>::process,
			SetTrueIfInRangesInclusiveVector<DataType, 16>::process };

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
			SetTrueIfInRangesExclusiveVector<DataType, 0>::process,
			SetTrueIfInRangesExclusiveVector<DataType, 1>::process,
			SetTrueIfInRangesExclusiveVector<DataType, 2>::process,
			SetTrueIfInRangesExclusiveVector<DataType, 3>::process,
			SetTrueIfInRangesExclusiveVector<DataType, 4>::process,
			SetTrueIfInRangesExclusiveVector<DataType, 5>::process,
			SetTrueIfInRangesExclusiveVector<DataType, 6>::process,
			SetTrueIfInRangesExclusiveVector<DataType, 7>::process,
			SetTrueIfInRangesExclusiveVector<DataType, 8>::process,
			SetTrueIfInRangesExclusiveVector<DataType, 9>::process,
			SetTrueIfInRangesExclusiveVector<DataType, 10>::process,
			SetTrueIfInRangesExclusiveVector<DataType, 11>::process,
			SetTrueIfInRangesExclusiveVector<DataType, 12>::process,
			SetTrueIfInRangesExclusiveVector<DataType, 13>::process,
			SetTrueIfInRangesExclusiveVector<DataType, 14>::process,
			SetTrueIfInRangesExclusiveVector<DataType, 15>::process,
			SetTrueIfInRangesExclusiveVector<DataType, 16>::process };

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

/* Check parameter arguments. Test data and result arrays*/
template<typename DataType>
bool IsValidDataAndResult(DataType const data[], bool const result[]) {
	if (!IsValidArray(data) || !IsValidArray(result)) {
		return false;
	}
	return true;
}

/* Check parameter arguments. Test range parameters */
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
 *  In case an exception is thrown during operation, it aborts if assertion is
 *  enabled. if not, returns kUnknownError status.
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

	CHECK_ARGS(
			IsValidDataAndResult(data, result)
					&& IsValidBounds(num_condition, lower_bounds,
							upper_bounds));

	try {
		func(/*num_data, data, num_condition, lower_bounds, upper_bounds, result*/);
	} catch (...) {
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
 *  In case an exception is thrown during operation, it aborts if assertion is
 *  enabled. if not, returns kUnknownError status.
 *
 *  @param[in] func A function to return a boolean value for a @a data element.
 *  It should take a @a data element and return a boolean.
 */
template<typename Func, typename DataType>
LIBSAKURA_SYMBOL(Status) DoElementFuncBoolFilter(Func func, size_t num_data,
		DataType const data[/*num_data*/], bool result[/*num_data*/]) {

	CHECK_ARGS(IsValidDataAndResult(data, result));

	try {
		auto adata = AssumeAligned(data);
		auto aresult = AssumeAligned(result);
		// No operation is done when num_data==0.
		for (size_t i = 0; i < num_data; ++i) {
			aresult[i] = func(adata[i]);
		}
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

/*
 *  @brief Generate a boolean filter based on per array equation function @a func .
 *
 *  This function calls a function, @a func , with a @a data array.
 *  In case an exception is thrown during operation, it aborts if assertion is
 *  enabled. if not, returns kUnknownError status.
 *
 *  @param[in] func A function to return a boolean array, @a result, based on a @a data array.
 *  It should take a @a data array and output a @result array.
 */
template<typename Func, typename DataType>
LIBSAKURA_SYMBOL(Status) DoArrayFuncBoolFilter(Func func, size_t num_data,
		DataType const data[/*num_data*/], bool const result[/*num_data*/]) {

	CHECK_ARGS(IsValidDataAndResult(data, result));

	try {
		func(/*num_data, data, result*/);
	} catch (...) {
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
	auto operation_for_element = [threshold](decltype(data[0]) data_value) {
		return (data_value > threshold);
	};
	return DoElementFuncBoolFilter(operation_for_element, num_data, data,
			result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) noexcept {
	auto operation_for_element = [threshold](decltype(data[0]) data_value) {
		return (data_value > threshold);
	};
	return DoElementFuncBoolFilter(operation_for_element, num_data, data,
			result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanOrEqualsFloat)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]) noexcept {
	auto operation_for_element = [threshold](decltype(data[0]) data_value) {
		return (data_value >= threshold);
	};
	return DoElementFuncBoolFilter(operation_for_element, num_data, data,
			result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanOrEqualsInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) noexcept {
	auto operation_for_element = [threshold](decltype(data[0]) data_value) {
		return (data_value >= threshold);
	};
	return DoElementFuncBoolFilter(operation_for_element, num_data, data,
			result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanFloat)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]) noexcept {
	auto operation_for_element = [threshold](decltype(data[0]) data_value) {
		return (data_value < threshold);
	};
	return DoElementFuncBoolFilter(operation_for_element, num_data, data,
			result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) noexcept {
	auto operation_for_element = [threshold](decltype(data[0]) data_value) {
		return (data_value < threshold);
	};
	return DoElementFuncBoolFilter(operation_for_element, num_data, data,
			result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanOrEqualsFloat)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]) noexcept {
	auto operation_for_element = [threshold](decltype(data[0]) data_value) {
		return (data_value <= threshold);
	};
	return DoElementFuncBoolFilter(operation_for_element, num_data, data,
			result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanOrEqualsInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) noexcept {
	auto operation_for_element = [threshold](decltype(data[0]) data_value) {
		return (data_value <= threshold);
	};
	return DoElementFuncBoolFilter(operation_for_element, num_data, data,
			result);
}

#if defined(__AVX2__)
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetFalseIfNanOrInfFloat)(
		size_t num_data, float const data[/*num_data*/],
		bool result[/*num_data*/]) noexcept {
	CHECK_ARGS(IsValidDataAndResult(data, result));

	constexpr uint32_t kExponetMask = 0x7F800000;
	STATIC_ASSERT(sizeof(kExponetMask) == sizeof(data[0]));

	constexpr int32_t kZero = 0x80808080;
	constexpr int32_t kLSBytes = 0x0c080400;
	const auto idx0 = _mm256_set_epi32(kZero, kZero, kZero, kLSBytes,
			kZero, kZero, kZero, kLSBytes);
	const auto idx1 = _mm256_set_epi32(kZero, kZero, kLSBytes, kZero,
			kZero, kZero, kLSBytes, kZero);
	const auto idx2 = _mm256_set_epi32(kZero, kLSBytes, kZero, kZero,
			kZero, kLSBytes, kZero, kZero);
	const auto idx3 = _mm256_set_epi32(kLSBytes, kZero, kZero, kZero,
			kLSBytes, kZero, kZero, kZero);
	const auto shuffle = _mm256_set_epi32(7, 3, 6, 2, 5, 1, 4, 0);

	constexpr auto kElementsPerLoop = sizeof(LIBSAKURA_SYMBOL(SimdPacketAVX)) / LIBSAKURA_SYMBOL(SimdPacketAVX)::kNumFloat;
	const auto data_ptr = AssumeAligned(reinterpret_cast<__m256i const (*)[kElementsPerLoop]>(data));
	const auto result_ptr = AssumeAligned<__m256i *>(reinterpret_cast<__m256i *>(result));
	const auto nan_inf_mask = _mm256_set1_epi32(kExponetMask);
	const auto one = _mm256_set1_epi32(true);
	const auto n = num_data / LIBSAKURA_SYMBOL(SimdPacketAVX)::kNumFloat / kElementsPerLoop;
	for (size_t i = 0; i < n; ++i) {
		auto is_valid = _mm256_shuffle_epi8(
				_mm256_andnot_si256(
						_mm256_cmpeq_epi32(nan_inf_mask,
								_mm256_and_si256(data_ptr[i][0], nan_inf_mask)),
						one), idx0);
		is_valid = _mm256_or_si256(
				_mm256_shuffle_epi8(
						_mm256_andnot_si256(
								_mm256_cmpeq_epi32(nan_inf_mask,
										_mm256_and_si256(data_ptr[i][1],
												nan_inf_mask)), one), idx1),
				is_valid);
		is_valid = _mm256_or_si256(
				_mm256_shuffle_epi8(
						_mm256_andnot_si256(
								_mm256_cmpeq_epi32(nan_inf_mask,
										_mm256_and_si256(data_ptr[i][2],
												nan_inf_mask)), one), idx2),
				is_valid);
		is_valid = _mm256_or_si256(
				_mm256_shuffle_epi8(
						_mm256_andnot_si256(
								_mm256_cmpeq_epi32(nan_inf_mask,
										_mm256_and_si256(data_ptr[i][3],
												nan_inf_mask)), one), idx3),
				is_valid);
		result_ptr[i] = _mm256_permutevar8x32_epi32(is_valid, shuffle);
	}

	const auto end = n * LIBSAKURA_SYMBOL(SimdPacketAVX)::kNumFloat * kElementsPerLoop;
	data = &data[end];
	result = &result[end];
	num_data -= end;
	auto operation_for_element = [](decltype(data[0]) data_value) -> bool {
		union {float fvalue; int ivalue;}value;
		value.fvalue = data_value;
		return ((value.ivalue & kExponetMask) != kExponetMask);
	};
	return DoElementFuncBoolFilter(operation_for_element, num_data, data, result);
}
#else
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetFalseIfNanOrInfFloat)(
		size_t num_data, float const data[/*num_data*/],
		bool result[/*num_data*/]) noexcept {
	constexpr uint32_t kExponetMask = 0x7F800000;
	STATIC_ASSERT(sizeof(kExponetMask) == sizeof(data[0]));
	// code 1'
	auto operation_for_element = [](decltype(data[0]) data_value) -> bool {
		union {float fvalue; int ivalue;}value;
		value.fvalue = data_value;
		return ((value.ivalue & kExponetMask) != kExponetMask);
	};
	return DoElementFuncBoolFilter(operation_for_element, num_data, data,
			result);
	// code 2'
//	uint32_t const *data_int = reinterpret_cast<uint32_t const *>(data);
//	auto operation_for_element = [](decltype(data_int[0]) data_value) -> bool {
//		return ((data_value & kExponetMask) != kExponetMask);
//	};
//	return DoElementFuncBoolFilter(operation_for_element, num_data, data_int, result);
}
#endif
#if defined(__AVX2__)
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint8ToBool)(
		size_t num_data, uint8_t const data[/*num_data*/],
		bool result[/*num_data*/]) noexcept {
	CHECK_ARGS(IsValidDataAndResult(data, result));
	STATIC_ASSERT(sizeof(data[0])==sizeof(result[0]));
	STATIC_ASSERT(true==1);
	STATIC_ASSERT(false==0);

	constexpr uint8_t kZero = 0;
	try {
		constexpr auto kElementsPerLoop = LIBSAKURA_SYMBOL(SimdPacketAVX)::kSize / sizeof(data[0]);
		const auto data_ptr = AssumeAligned(reinterpret_cast<__m256i const *>(data));
		const auto result_ptr = AssumeAligned(reinterpret_cast<__m256i *>(result));
		const auto zero = _mm256_set1_epi8(kZero);
		const auto one = _mm256_set1_epi8(1);
		const auto n = num_data / kElementsPerLoop;
		for (size_t i = 0; i < n; ++i) {
			//Returns 0xFF if data==zero, else 0x00
			auto mask = _mm256_cmpeq_epi8(data_ptr[i], zero);
			result_ptr[i] = _mm256_add_epi8(mask, one);
		}
		// process remaining elements
		const auto end = n * kElementsPerLoop;
		data = &data[end];
		result = &result[end];
		num_data -= end;
		auto operation_for_element = [kZero](decltype(data[0]) data_value) {
			return (data_value != kZero);
		};
		return DoElementFuncBoolFilter(operation_for_element, num_data, data, result);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
}
#else
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint8ToBool)(
		size_t num_data, uint8_t const data[/*num_data*/],
		bool result[/*num_data*/]) noexcept {
	constexpr uint8_t kZero = 0;
	auto operation_for_element = [kZero](decltype(data[0]) data_value) {
		return (data_value != kZero);
	};
	return DoElementFuncBoolFilter(operation_for_element, num_data, data,
			result);
}
#endif

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint32ToBool)(
		size_t num_data, uint32_t const data[/*num_data*/],
		bool result[/*num_data*/]) noexcept {
	constexpr uint8_t kZero = 0;
	auto operation_for_element = [kZero](decltype(data[0]) data_value) {
		return (data_value != kZero);
	};
	return DoElementFuncBoolFilter(operation_for_element, num_data, data,
			result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InvertBool)(
		size_t num_data,
		bool const data[/*num_data*/], bool result[/*num_data*/]) noexcept {
	return DoArrayFuncBoolFilter([=] {InvertBool(num_data, data, result);},
			num_data, data, result);
}
