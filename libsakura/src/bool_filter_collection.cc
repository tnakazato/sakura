#include <cassert>
#include <iostream>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory_impl.h"
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

template<typename DataType>
inline void SetTrueGreaterThan(size_t num_data, DataType const *data,
		DataType threshold, bool *result) {
	DataType const zero(static_cast<DataType>(0));
	for (size_t i = 0; i < num_data; ++i) {
		result[i] = (data[i] - threshold) > zero;
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

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {
template<typename DataType>
void ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<DataType>::SetTrueInRangesInclusive(
		size_t num_data, DataType const data[/*num_data*/],
		size_t num_condition, DataType const lower_bounds[/*num_condition*/],
		DataType const upper_bounds[/*num_condition*/],
		bool result[/*num_data*/]) const {
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
void ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<DataType>::SetTrueGreaterThan(
		size_t num_data, DataType const data[/*num_data*/], DataType threshold,
		bool result[/*num_data*/]) const {
	::SetTrueGreaterThan(num_data, data, threshold, result);
}

template<typename DataType>
void ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<DataType>::ToBool(
		size_t num_data, DataType const data[/*num_data*/],
		bool result[/*num_data*/]) const {
	::ToBool(num_data, data, result);
}

template<typename DataType>
void ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<DataType>::InvertBool(
		size_t num_data,
		bool const data[/*num_data*/], bool result[/*num_data*/]) const {
	::InvertBool(num_data, data, result);
}

template class ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<uint8_t> ;
template class ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<uint32_t> ;
template class ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<float> ;
template class ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<int> ;
}
