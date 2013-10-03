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

template<class DataType>
inline void ApplyPositionSwitchCalibration(size_t num_scaling_factor,
		DataType const scaling_factor[/*num_scaling_factor*/], size_t num_data,
		DataType const target[/*num_data*/], DataType const reference[/*num_data*/],
		DataType result[/*num_data*/]) {
	assert(num_scaling_factor > 0);
	assert(num_scaling_factor == 1 || num_scaling_factor >= num_data);
	assert(scaling_factor != nullptr);
	assert(target != nullptr);
	assert(reference != nullptr);
	assert(result != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(scaling_factor));
	assert(LIBSAKURA_SYMBOL(IsAligned)(target));
	assert(LIBSAKURA_SYMBOL(IsAligned)(reference));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));

	if (num_scaling_factor == 1) {
		DataType const constant_scaling_factor = scaling_factor[0];
		for (size_t i = 0; i < num_data; ++i) {
			result[i] = constant_scaling_factor * (target[i] - reference[i]) / reference[i];
		}
	}
	else {
		for (size_t i = 0; i < num_data; ++i) {
			result[i] = scaling_factor[i] * (target[i] - reference[i]) / reference[i];
		}
	}
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {
template<class DataType>
void ADDSUFFIX(ApplyCalibration, ARCH_SUFFIX)<DataType>::ApplyPositionSwitchCalibration(
		size_t num_scaling_factor,
		DataType const scaling_factor[/*num_scaling_factor*/], size_t num_data,
		DataType const target[/*num_data*/], DataType const reference[/*num_data*/],
		DataType result[/*num_data*/]) const {
	::ApplyPositionSwitchCalibration(num_scaling_factor, scaling_factor,
			num_data, target, reference, result);
}

template class ADDSUFFIX(ApplyCalibration, ARCH_SUFFIX)<float> ;
} /* namespace LIBSAKURA_PREFIX */
