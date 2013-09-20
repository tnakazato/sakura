#include <cassert>
#include <iostream>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/localdef.h"

#define FORCE_SCALAR 0

// Vectorization by Compiler
namespace {

template<typename DataType>
inline void SetTrueInRangesInclusiveScalar(size_t num_data,
		DataType const *data, size_t num_condition,
		DataType const *lower_bounds, DataType const *upper_bounds,
		bool *result) {
	// So far, only unit8_t version is vectorized
	//std::cout << "Invoking SetTrueInRangesInclusiveDefault()" << std::endl;
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	assert(LIBSAKURA_SYMBOL(IsAligned)(upper_bounds));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lower_bounds));
	static_assert(true == 1, "true==1");
	static_assert(false == 0, "false==0");
	// Initialize result with false
	for (size_t i = 0; i < num_data; ++i) {
		result[i] = false;
	}
	DataType lower_value, upper_value;

	for (size_t j = 0; j < num_condition; ++j) {
		lower_value = lower_bounds[j];
		upper_value = upper_bounds[j];
		for (size_t i = 0; i < num_data; ++i) {
			result[i] = result[i]
					|| ((data[i] - lower_value) * (upper_value - data[i]) >= 0);
		}
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

inline void InvertBoolScalar(size_t num_data, bool const *data, bool *result) {
//	std::cout << "Invoking InvertBoolDefault()" << std::endl;
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	uint8_t true8(static_cast<uint8_t>(true));
	static_assert(sizeof(data[0]) == sizeof(true8), "sizeof(bool)==sizeof(uint8_t)");
	static_assert(true == 1, "true==1");
	static_assert(false == 0, "false==0");
	for (size_t i = 0; i < num_data; ++i) {
		result[i] = (data[i] ^ true8);
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
	::SetTrueInRangesInclusiveScalar(num_data, data, num_condition,
			lower_bounds, upper_bounds, result);
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
	::InvertBoolScalar(num_data, data, result);
}

template class ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<uint8_t> ;
template class ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<uint32_t> ;
template class ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<float> ;
template class ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<int> ;
}
