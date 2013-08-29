#include <cassert>
#include <iostream>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/localdef.h"

#define FORCE_SCALAR 1

#if FORCE_SCALAR
// Scalar implementation
namespace {

template<typename DataType>
inline void SetTrueInRangesInclusiveScalar(size_t num_data, DataType const *data,
		size_t num_condition, DataType const *lower_bounds,
		DataType const *upper_bounds, bool *result) {

	std::cout << "Invoking SetTrueInRangesInclusiveScalar()" << std::endl;
//	for (size_t i=0; i < num_in ; ++i){
//		out[i] = edit_mask[i] ? (in[i] & bit_mask) : in[i];
//	}

}

} /* anonymous namespace */

#else /* FORCE_SCALAR */

// Vectorization by Compiler
namespace {

template<typename DataType>
inline void SetTrueInRangesInclusiveDefault(size_t num_data, DataType const *data,
		size_t num_condition, DataType const *lower_bounds,
		DataType const *upper_bounds, bool *result) {

	std::cout << "Invoking SetTrueInRangesInclusiveDefault()" << std::endl;
//	assert(LIBSAKURA_SYMBOL(IsAligned)(in));
//	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
//	assert(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
//	// cast bool array to uint8_t array
//	uint8_t const *mask8 = reinterpret_cast<uint8_t const *>(edit_mask);
//	assert(LIBSAKURA_SYMBOL(IsAligned)(mask8));
//	assert(sizeof(edit_mask[0]) == sizeof(mask8[0]));
//	assert(true == 1);
//	assert(false == 0);
//
//	/* edit_mask = true: (mask8 - 1) = 00...0
//	 *                   -> (bit_mask | 00...0) = bit_mask,
//	 *           = false: (mask8 - 1) = 11...1
//	 *                   -> (bit_mask | 11...1) = 11...1 */
//	for (size_t i=0; i < num_in ; ++i){
//		out[i] = in[i] & (bit_mask | (static_cast<DataType>(mask8[i]) - 1));
//	}

}

} /* anonymous namespace */

#endif /* FORCE_SCALAR */

namespace LIBSAKURA_PREFIX {
template<typename DataType>
void ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<DataType>::SetTrueInRangesInclusive(size_t num_data,
		DataType const data[/*num_data*/], size_t num_condition,
		DataType const lower_bounds[/*num_condition*/],
		DataType const upper_bounds[/*num_condition*/],
		bool result[/*num_data*/]) const {
#if FORCE_SCALAR
	SetTrueInRangesInclusiveScalar(num_data, data, num_condition,
			lower_bounds, upper_bounds, result);
#else
	SetTrueInRangesInclusiveDefault(num_data, data, num_condition,
			lower_bounds, upper_bounds, result);
#endif
}

template class ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<float>;
template class ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<int>;
}
