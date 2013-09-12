#include <cassert>
#include <iostream>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/localdef.h"

#define FORCE_SCALAR 0

#if FORCE_SCALAR
// Scalar implementation
namespace {

template<typename DataType>
inline void SetTrueInRangesInclusiveScalar(size_t num_data, DataType const *data,
		size_t num_condition, DataType const *lower_bounds,
		DataType const *upper_bounds, bool *result) {
	std::cout << "Invoking SetTrueInRangesInclusiveScalar()" << std::endl;
	// Initialize result with false
	for (size_t i=0; i < num_data ; ++i){
			result[i] = false;
		}
	DataType lower_value, upper_value;
	for (size_t i=0; i < num_data ; ++i){
		//std::cout << "i = " << i << std::endl;
		for (size_t j=0; j < num_condition; ++j){
			lower_value = lower_bounds[j];
			upper_value = upper_bounds[j];
			//std::cout << "Searching range [ " << lower_value << ", " << upper_value << " ]" << std::endl;
			if ( result[i] ) continue;
			result[i] = ( (data[i] - lower_value) * (upper_value - data[i]) >= 0 );
		}
	}
}

template<typename DataType>
inline void ToBoolScalar(size_t num_data, DataType const *in, bool *out) {
//	std::cout << "Invoking ToBoolScalar()" << std::endl;
	DataType const zero(static_cast<DataType>(0));
	for (size_t i=0; i < num_data ; ++i){
			out[i] = (in[i] != zero);
		}
}

inline void InvertBoolScalar(size_t num_data, bool const *in, bool *out) {
//	std::cout << "Invoking InvertBoolScalar()" << std::endl;
	for (size_t i=0; i < num_data ; ++i){
			out[i] = !in[i];
		}
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
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	assert(LIBSAKURA_SYMBOL(IsAligned)(upper_bounds));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lower_bounds));
	assert(true == 1);
	assert(false == 0);
	// Initialize result with false
	for (size_t i=0; i < num_data ; ++i){
		result[i] = false;
	}
	DataType lower_value, upper_value;


	for (size_t j=0; j < num_condition; ++j){
		lower_value = lower_bounds[j];
		upper_value = upper_bounds[j];
		for (size_t i=0; i < num_data ; ++i){
			result[i] = result[i] || ( (data[i] - lower_value) * (upper_value - data[i]) >= 0 );
		}
	}
//	for (size_t i=0; i < num_data ; ++i){
//		for (size_t j=0; j < num_condition; ++j){
//			lower_value = lower_bounds[j];
//			upper_value = upper_bounds[j];
//			result[i] = result[i] || ( (data[i] - lower_value) * (upper_value - data[i]) >= 0 );
//		}
//	}
//	// cast bool array to uint8_t array
//	uint8_t const *mask8 = reinterpret_cast<uint8_t const *>(edit_mask);
//	assert(sizeof(edit_mask[0]) == sizeof(mask8[0]));
//
//	/* edit_mask = true: (mask8 - 1) = 00...0
//	 *                   -> (bit_mask | 00...0) = bit_mask,
//	 *           = false: (mask8 - 1) = 11...1
//	 *                   -> (bit_mask | 11...1) = 11...1 */
//	for (size_t i=0; i < num_in ; ++i){
//		out[i] = in[i] & (bit_mask | (static_cast<DataType>(mask8[i]) - 1));
//	}

}

template<typename DataType>
inline void ToBoolDefault(size_t num_data, DataType const *in, bool *out) {
//	std::cout << "Invoking ToBoolDefault()" << std::endl;
	assert(LIBSAKURA_SYMBOL(IsAligned)(in));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	DataType const zero(static_cast<DataType>(0));
	for (size_t i=0; i < num_data ; ++i){
			out[i] = (in[i] != zero);
		}
}

inline void InvertBoolDefault(size_t num_data, bool const *in, bool *out) {
//	std::cout << "Invoking InvertBoolDefault()" << std::endl;
	assert(LIBSAKURA_SYMBOL(IsAligned)(in));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	assert(true == 1);
	assert(false == 0);
	for (size_t i=0; i < num_data ; ++i){
			out[i] = (in[i] ^ true);
		}
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

template<typename DataType>
void ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<DataType>::ToBool(size_t num_in,
		DataType const in[/*num_in*/], bool out[/*num_in*/]) const {
#if FORCE_SCALAR
	ToBoolScalar(num_in, in, out);
#else
	ToBoolDefault(num_in, in, out);
#endif
}

template<typename DataType>
void ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<DataType>::InvertBool(size_t num_in,
		bool const in[/*num_in*/], bool out[/*num_in*/]) const {
#if FORCE_SCALAR
	InvertBoolScalar(num_in, in, out);
#else
	InvertBoolDefault(num_in, in, out);
#endif
}

template class ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<uint8_t>;
template class ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<uint32_t>;
template class ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<float>;
template class ADDSUFFIX(BoolFilterCollection, ARCH_SUFFIX)<int>;
}
