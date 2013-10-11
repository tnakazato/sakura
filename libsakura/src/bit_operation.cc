#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/localdef.h"

// Vectorization by Compiler
namespace {

template<typename DataType>
inline void OperateBitsAnd(DataType bit_mask, size_t num_data, DataType const *data,
bool const *edit_mask, DataType *result) {

	//std::cout << "Invoking OperateBitsAndDefault()" << std::endl;
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	assert(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
	// cast bool array to uint8_t array
	uint8_t const *mask8 = reinterpret_cast<uint8_t const *>(edit_mask);
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask8));
	STATIC_ASSERT(sizeof(edit_mask[0]) == sizeof(mask8[0]));
	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);

	/* edit_mask = true: (mask8 - 1) = 00...0
	 *                   -> (bit_mask | 00...0) = bit_mask,
	 *           = false: (mask8 - 1) = 11...1
	 *                   -> (bit_mask | 11...1) = 11...1 */
	for (size_t i = 0; i < num_data; ++i) {
		result[i] = data[i] & (bit_mask | (static_cast<DataType>(mask8[i]) - 1));
	}

}

template<typename DataType>
inline void OperateBitsConverseNonImplication(DataType bit_mask, size_t num_data,
		DataType const *data,bool const *edit_mask, DataType *result) {

	//std::cout << "Invoking OperateBitsAndDefault()" << std::endl;
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	assert(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
	// cast bool array to uint8_t array
	uint8_t const *mask8 = reinterpret_cast<uint8_t const *>(edit_mask);
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask8));
	STATIC_ASSERT(sizeof(edit_mask[0]) == sizeof(mask8[0]));
	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);

	/* edit_mask = true: X = (mask8 - 1) = 00...0
	 *                   -> (bit_mask | 00...0) = bit_mask,
	 *                   ~X = 11...1 -> (data[i] ^ ~X) = ~data[i]
	 *           = false: X = (mask8 - 1) = 11...1
	 *                   -> (bit_mask | 11...1) = 11...1
	 *                   ~X = 00...0 -> (data[i] ^ ~X) = data[i] */
	for (size_t i = 0; i < num_data; ++i) {
		DataType X = (static_cast<DataType>(mask8[i]) - 1);
		result[i] = (data[i] ^ ~X) & (bit_mask | X);
	}

}

template<typename DataType>
inline void OperateBitsOr(DataType bit_mask, size_t num_data, DataType const *data,
bool const *edit_mask, DataType *result) {

	//std::cout << "Invoking OperateBitsAndDefault()" << std::endl;
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	assert(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
	// cast bool array to uint8_t array
	uint8_t const *mask8 = reinterpret_cast<uint8_t const *>(edit_mask);
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask8));
	STATIC_ASSERT(sizeof(edit_mask[0]) == sizeof(mask8[0]));
	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);

	/* edit_mask = true: ~(mask8 - 1) = 11...1
	 *                   -> (bit_mask & 11...1) = bit_mask,
	 *           = false: ~(mask8 - 1) = 00...0
	 *                   -> (bit_mask & 00...0) = 00...0 */
	for (size_t i = 0; i < num_data; ++i) {
		result[i] = data[i] | (bit_mask & ~(static_cast<DataType>(mask8[i]) - 1));
	}

}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {
template<typename DataType>
void ADDSUFFIX(BitOperation, ARCH_SUFFIX)<DataType>::OperateBitsAnd(
		DataType bit_mask, size_t num_data, DataType const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], DataType result[/*num_data*/]) const {
	::OperateBitsAnd(bit_mask, num_data, data, edit_mask, result);
}

template<typename DataType>
void ADDSUFFIX(BitOperation, ARCH_SUFFIX)<DataType>::OperateBitsConverseNonImplication(
		DataType bit_mask, size_t num_data, DataType const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], DataType result[/*num_data*/]) const {
	::OperateBitsConverseNonImplication(bit_mask, num_data, data, edit_mask, result);
}

template<typename DataType>
void ADDSUFFIX(BitOperation, ARCH_SUFFIX)<DataType>::OperateBitsOr(
		DataType bit_mask, size_t num_data, DataType const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], DataType result[/*num_data*/]) const {
	::OperateBitsOr(bit_mask, num_data, data, edit_mask, result);
}

template class ADDSUFFIX(BitOperation, ARCH_SUFFIX)<uint8_t> ;
template class ADDSUFFIX(BitOperation, ARCH_SUFFIX)<uint32_t> ;
}
