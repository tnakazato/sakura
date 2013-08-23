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

void OperateLogicalAndSimd(size_t num_in, bool const in1[],
		bool const in2[], bool out[]) {
	std::cout << "OperateLogicalAndSimd function is called. This function is not implemented yet." << std::endl;
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

inline void OperateLogicalAndEigen(size_t num_in, bool const *in1,
		bool const *in2, bool *out) {
	//std::cout << "OperateLogicalAndEigen function is called" << std::endl;

	assert(LIBSAKURA_SYMBOL(IsAligned)(in1));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in2));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	Map<Array<bool, Dynamic, 1>, Aligned> in1_(const_cast<bool *>(in1),
			num_in);
	Map<Array<bool, Dynamic, 1>, Aligned> in2_(const_cast<bool *>(in2),
			num_in);
	Map<Array<bool, Dynamic, 1>, Aligned> out_(const_cast<bool *>(out),
			num_in);

	out_ = in1_ && in2_;
}

} /* anonymous namespace */

#endif /* defined(__AVX__) */

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(LogicalOperation, ARCH_SUFFIX)::OperateLogicalAnd(size_t num_in,
		bool const in1[/*num_in*/], bool const in2[/*num_in*/],
		bool out[/*num_in*/]) const {
#if defined( __AVX__) && (! FORCE_EIGEN)
	OperateLogicalAndSimd(num_in, in1, in2, out);
#else
	OperateLogicalAndEigen(num_in, in1, in2, out);
#endif
}

}
