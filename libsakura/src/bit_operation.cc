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
inline void OperateBitsAnd(DataType bit_mask, size_t num_in, DataType const *in,
bool const *edit_mask, DataType *out) {

	//std::cout << "Invoking OperateBitsAndDefault()" << std::endl;
	assert(LIBSAKURA_SYMBOL(IsAligned)(in));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	assert(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
	// cast bool array to uint8_t array
	uint8_t const *mask8 = reinterpret_cast<uint8_t const *>(edit_mask);
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask8));
	static_assert(sizeof(edit_mask[0]) == sizeof(mask8[0]), "sizeof(bool)==sizeof(uint8_t)");
	static_assert(true == 1, "true==1");
	static_assert(false == 0, "false==0");

	/* edit_mask = true: (mask8 - 1) = 00...0
	 *                   -> (bit_mask | 00...0) = bit_mask,
	 *           = false: (mask8 - 1) = 11...1
	 *                   -> (bit_mask | 11...1) = 11...1 */
	for (size_t i = 0; i < num_in; ++i) {
		out[i] = in[i] & (bit_mask | (static_cast<DataType>(mask8[i]) - 1));
	}

}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {
template<typename DataType>
void ADDSUFFIX(BitOperation, ARCH_SUFFIX)<DataType>::OperateBitsAnd(
		DataType bit_mask, size_t num_in, DataType const in[/*num_in*/],
		bool const edit_mask[/*num_in*/], DataType out[/*num_in*/]) const {
	::OperateBitsAnd(bit_mask, num_in, in, edit_mask, out);
}

template class ADDSUFFIX(BitOperation, ARCH_SUFFIX)<uint8_t> ;
template class ADDSUFFIX(BitOperation, ARCH_SUFFIX)<uint32_t> ;
}
