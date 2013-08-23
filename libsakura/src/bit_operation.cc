#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/localdef.h"

#define FORCE_SCALAR 0

#if FORCE_SCALAR
// Scalar implementation
namespace {

template<typename DataType>
inline void OperateBitsAndScalar(DataType bit_mask, size_t num_in, DataType const *in,
		bool const *edit_mask, DataType *out) {

	//std::cout << "Invoking OperateBitsAndScalar()" << std::endl;
	assert(LIBSAKURA_SYMBOL(IsAligned)(in));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	assert(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
	for (size_t i=0; i < num_in ; ++i){
		out[i] = edit_mask[i] ? (in[i] & bit_mask) : in[i];
	}

}

} /* anonymous namespace */

#else /* FORCE_SCALAR */

// Vectorization by Compiler
namespace {

template<typename DataType>
inline void OperateBitsAndDefault(DataType bit_mask, size_t num_in, DataType const *in,
		bool const *edit_mask, DataType *out) {

	//std::cout << "Invoking OperateBitsAndDefault()" << std::endl;
	assert(LIBSAKURA_SYMBOL(IsAligned)(in));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	assert(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
	// cast bool array to uint8_t array
	uint8_t const *mask8 = reinterpret_cast<uint8_t const *>(edit_mask);
	assert(sizeof(edit_mask[0]) == sizeof(mask8[0]));
	assert(true == 1);
	assert(false == 0);

	/* edit_mask = true: (mask8 - 1) = 00...0
	 *                   -> (bit_mask | 00...0) = bit_mask,
	 *           = false: (mask8 - 1) = 11...1
	 *                   -> (bit_mask | 11...1) = 11...1 */
	for (size_t i=0; i < num_in ; ++i){
		out[i] = in[i] & (bit_mask | (static_cast<DataType>(mask8[i]) - 1));
	}

}

} /* anonymous namespace */

#endif /* FORCE_SCALAR */

namespace LIBSAKURA_PREFIX {
template<typename DataType>
void ADDSUFFIX(BitOperation, ARCH_SUFFIX)<DataType>::OperateBitsAnd(DataType bit_mask, size_t num_in,
		DataType const in[/*num_in*/], bool const edit_mask[/*num_in*/],
		DataType out[/*num_in*/]) const {
#if FORCE_SCALAR
	OperateBitsAndScalar(bit_mask, num_in, in, edit_mask, out);
#else
	OperateBitsAndDefault(bit_mask, num_in, in, edit_mask, out);
#endif
}

template class ADDSUFFIX(BitOperation, ARCH_SUFFIX)<uint8_t>;
template class ADDSUFFIX(BitOperation, ARCH_SUFFIX)<uint32_t>;
}

//#define FORCE_EIGEN 0
//
//#if defined(__AVX__) && (! FORCE_EIGEN)
//#include <immintrin.h>
//#include <cstdint>
//
//namespace {
//
//void OperateBitsAndSimd(uint8_t bit_mask, size_t num_in,
//		uint8_t const in[], bool const edit_mask[], uint8_t out[]) {
//	std::cout << "OperateBitsAndSimd function for uint8_t is called. This function is not implemented yet." << std::endl;
//}
//
//void OperateBitsAndSimd(uint32_t bit_mask, size_t num_in,
//		uint32_t const in[], bool const edit_mask[], uint32_t out[]) {
//	std::cout << "OperateBitsAndSimd function for uint32_t is called. This function is not implemented yet." << std::endl;
//}
//
//} /* anonymous namespace */
//
//#else /* defined(__AVX__) */

//#define EIGEN_DENSEBASE_PLUGIN "eigen_binary_visitor_plugin.h"
//#include <Eigen/Core>
//
//using ::Eigen::Map;
//using ::Eigen::Array;
//using ::Eigen::Dynamic;
//using ::Eigen::Aligned;
//
//namespace {
//
//template<typename DataType>
//inline void OperateBitsAndEigen(DataType bit_mask, size_t num_in, DataType const *in,
//		bool const *edit_mask, DataType *out) {
//	//std::cout << "OperateBitsAndEigen function is called" << std::endl;
//
//	assert(LIBSAKURA_SYMBOL(IsAligned)(in));
//	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
//	assert(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
//	Map<Array<DataType, Dynamic, 1>, Aligned> in_(const_cast<DataType *>(in),
//			num_in);
//	Map<Array<DataType, Dynamic, 1>, Aligned> out_(const_cast<DataType *>(out),
//			num_in);
//	Map<Array<bool, Dynamic, 1>, Aligned> edit_mask_(
//			const_cast<bool *>(edit_mask), num_in);
//
//	for (size_t i=0; i < num_in ; ++i){
//		out[i] = edit_mask_[i] ? (in_[i] & bit_mask) : in_[i];
//	}
//	//out_ = edit_mask_ ? (in_ & bit_mask) : in_;
//}
//
//} /* anonymous namespace */
//
//#endif /* defined(__AVX__) */
//
//namespace LIBSAKURA_PREFIX {
//template<typename DataType>
//void ADDSUFFIX(BitOperation, ARCH_SUFFIX)<DataType>::OperateBitsAnd(DataType bit_mask, size_t num_in,
//		DataType const in[/*num_in*/], bool const edit_mask[/*num_in*/],
//		DataType out[/*num_in*/]) const {
//#if defined( __AVX__) && (! FORCE_EIGEN)
//	OperateBitsAndSimd(bit_mask, num_in, in, edit_mask, out);
//#else
//	OperateBitsAndEigen(bit_mask, num_in, in, edit_mask, out);
//#endif
//}
//
//}

