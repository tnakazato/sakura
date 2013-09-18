#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/localdef.h"

namespace {

void OperateLogicalAnd(size_t num_in, bool const *in1,
		bool const *in2, bool *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(in1));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in2));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	for (size_t i = 0; i < num_in; ++i) {
		out[i] = in1[i] && in2[i];
	}
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(LogicalOperation, ARCH_SUFFIX)::OperateLogicalAnd(size_t num_in,
		bool const in1[/*num_in*/], bool const in2[/*num_in*/],
		bool out[/*num_in*/]) const {
	::OperateLogicalAnd(num_in, in1, in2, out);
}

}
