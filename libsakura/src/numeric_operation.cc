#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/localdef.h"

#define FORCE_EIGEN 0

#if defined(__AVX__) && (! FORCE_EIGEN)
#include <immintrin.h>
#include <cstdint>

namespace {

void OperateFloatSubtractionSimd(size_t num_in, float const in1[],
		float const in2[], float out[]) {
	std::cout << "OperateFloatSubtractionSimd function is called. This function is not implemented yet." << std::endl;
}

} /* anonymous namespace */

#else /* defined(__AVX__) */

#define EIGEN_DENSEBASE_PLUGIN "eigen_binary_visitor_plugin.h"
#include <Eigen/Core>

using ::Eigen::Map;
using ::Eigen::Array;
using ::Eigen::Dynamic;
using ::Eigen::Aligned;

namespace {

inline void OperateFloatSubtractionEigen(size_t num_in, float const *in1,
		float const *in2, float *out) {
	//std::cout << "OperateFloatSubtractionEigen function is called" << std::endl;

	assert(LIBSAKURA_SYMBOL(IsAligned)(in1));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in2));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	Map<Array<float, Dynamic, 1>, Aligned> in1_(const_cast<float *>(in1),
			num_in);
	Map<Array<float, Dynamic, 1>, Aligned> in2_(const_cast<float *>(in2),
			num_in);
	Map<Array<float, Dynamic, 1>, Aligned> out_(const_cast<float *>(out),
			num_in);

	for (size_t i=0; i < num_in ; i++){
		out[i] = in1_[i] - in2_[i];
	}
}

} /* anonymous namespace */

#endif /* defined(__AVX__) */

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::OperateFloatSubtraction(size_t num_in,
		float const in1[/*num_in*/], float const in2[/*num_in*/],
		float out[/*num_in*/]) const {
#if defined( __AVX__) && (! FORCE_EIGEN)
	OperateFloatSubtractionSimd(num_in, in1, in2, out);
#else
	OperateFloatSubtractionEigen(num_in, in1, in2, out);
#endif
}

}
