#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/localdef.h"
#include "libsakura/packed_operation.h"

#include <Eigen/Core>

using ::Eigen::Map;
using ::Eigen::Array;
using ::Eigen::Dynamic;
using ::Eigen::Aligned;

namespace {

#if defined(__AVX__)
#include <immintrin.h>

template<typename Arch, typename DataType, typename Context>
struct PacketAction {
	static inline void prologue(Context *context) {
	}
	static inline void action(size_t idx,
			typename Arch::PacketType const *scaling_factor,
			typename Arch::PacketType const *target,
			typename Arch::PacketType const *reference,
			typename Arch::PacketType *result, Context *context) {
		*result = LIBSAKURA_SYMBOL(SimdMath)<Arch, DataType>::Mul(
		LIBSAKURA_SYMBOL(SimdMath)<Arch, DataType>::Mul(*scaling_factor,
		LIBSAKURA_SYMBOL(SimdMath)<Arch, DataType>::Sub(*target, *reference)),
		LIBSAKURA_SYMBOL(SimdMath)<Arch, DataType>::Reciprocal(*reference));
	}
	static inline void epilogue(Context *context) {
	}
};

template<typename ScalarType, typename Context>
struct ScalarAction {
	static inline void prologue(Context *context) {
	}
	static inline void action(size_t idx, ScalarType const*scaling_factor,
			ScalarType const*target, ScalarType const*reference,
			ScalarType *result, Context *context) {
		*result = *scaling_factor * (*target - *reference) / *reference;
	}
	static inline void epilogue(Context *context) {
	}
};

template<class DataType>
inline void ApplyPositionSwitchCalibrationSimd(size_t num_scaling_factor,
		DataType const scaling_factor[/*num_scaling_factor*/], size_t num_data,
		DataType const target[/*num_data*/],
		DataType const reference[/*num_data*/], DataType result[/*num_data*/]) {
	if (num_scaling_factor == 1) {
		DataType const constant_scaling_factor = scaling_factor[0];
		for (size_t i = 0; i < num_data; ++i) {
			result[i] = constant_scaling_factor * (target[i] - reference[i])
					/ reference[i];
		}
	} else {
		typedef LIBSAKURA_SYMBOL(SimdArchNative) Arch;
		LIBSAKURA_SYMBOL(SimdIterate)<Arch, DataType const, Arch,
				DataType const, Arch, DataType const, Arch, DataType,
				PacketAction<Arch, DataType, void>,
				ScalarAction<DataType, void>, void>(num_data, scaling_factor,
				target, reference, result, (void*) nullptr);
	}
}

template<>
inline void ApplyPositionSwitchCalibrationSimd<float>(
		size_t num_scaling_factor,
		float const scaling_factor[/*num_scaling_factor*/], size_t num_data,
		float const target[/*num_data*/], float const reference[/*num_data*/],
		float result[/*num_data*/]) {
	if (num_scaling_factor == 1) {
		float const constant_scaling_factor = scaling_factor[0];
		for (size_t i = 0; i < num_data; ++i) {
			result[i] = constant_scaling_factor * (target[i] - reference[i])
					/ reference[i];
		}
	} else {
		for (size_t i = 0; i < num_data; i += 8) {
			__m256 s = _mm256_load_ps(&scaling_factor[i]);
			__m256 t = _mm256_load_ps(&target[i]);
			__m256 f = _mm256_load_ps(&reference[i]);
			__m256 r = _mm256_div_ps(_mm256_mul_ps(s, _mm256_sub_ps(t, f)), f);
			_mm256_store_ps(&result[i], r);
		}
		for (size_t i = std::max((size_t) 8, num_data - num_data % 8);
				i < num_data; ++i) {
			result[i] = scaling_factor[i] * (target[i] - reference[i])
					/ reference[i];
		}
		// TODO: loop on num_data % 8
		// TODO: use rcp and mul instead of div
	}
}
#endif
template<class DataType>
inline void ApplyPositionSwitchCalibration(size_t num_scaling_factor,
		DataType const scaling_factor[/*num_scaling_factor*/], size_t num_data,
		DataType const target[/*num_data*/],
		DataType const reference[/*num_data*/], DataType result[/*num_data*/]) {
	if (num_scaling_factor == 1) {
		DataType const constant_scaling_factor = scaling_factor[0];
		for (size_t i = 0; i < num_data; ++i) {
			result[i] = constant_scaling_factor * (target[i] - reference[i])
					/ reference[i];
		}
	} else {
#if defined(__AVX__)
		ApplyPositionSwitchCalibrationSimd(num_scaling_factor, scaling_factor,
				num_data, target, reference, result);
#else
		for (size_t i = 0; i < num_data; ++i) {
			result[i] = scaling_factor[i] * (target[i] - reference[i])
			/ reference[i];
		}
#endif
	}
}

template<class DataType>
inline void ApplyPositionSwitchCalibrationEigen(size_t num_scaling_factor,
		DataType const scaling_factor[/*num_scaling_factor*/], size_t num_data,
		DataType const target[/*num_data*/],
		DataType const reference[/*num_data*/], DataType result[/*num_data*/]) {
	if (num_scaling_factor == 1) {
		Map<Array<DataType, Dynamic, 1>, Aligned> eigen_result(result,
				num_data);
		Map<Array<DataType, Dynamic, 1>, Aligned> eigen_target(
				const_cast<DataType *>(target), num_data);
		Map<Array<DataType, Dynamic, 1>, Aligned> eigen_reference(
				const_cast<DataType *>(reference), num_data);
		DataType const constant_scaling_factor = scaling_factor[0];
		eigen_result = constant_scaling_factor
				* (eigen_target - eigen_reference) / eigen_reference;
	} else {
		Map<Array<DataType, Dynamic, 1>, Aligned> eigen_result(result,
				num_data);
		Map<Array<DataType, Dynamic, 1>, Aligned> eigen_target(
				const_cast<DataType *>(target), num_data);
		Map<Array<DataType, Dynamic, 1>, Aligned> eigen_reference(
				const_cast<DataType *>(reference), num_data);
		Map<Array<DataType, Dynamic, 1>, Aligned> eigen_scaling_factor(
				const_cast<DataType *>(scaling_factor), num_data);
		eigen_result = eigen_scaling_factor * (eigen_target - eigen_reference)
				/ eigen_reference;
	}
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {
template<class DataType>
void ADDSUFFIX(ApplyCalibration, ARCH_SUFFIX)<DataType>::ApplyPositionSwitchCalibration(
		size_t num_scaling_factor,
		DataType const scaling_factor[/*num_scaling_factor*/], size_t num_data,
		DataType const target[/*num_data*/],
		DataType const reference[/*num_data*/],
		DataType result[/*num_data*/]) const {
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
	::ApplyPositionSwitchCalibration(num_scaling_factor, scaling_factor,
			num_data, target, reference, result);
//	::ApplyPositionSwitchCalibrationEigen(num_scaling_factor, scaling_factor,
//			num_data, target, reference, result);
}

template class ADDSUFFIX(ApplyCalibration, ARCH_SUFFIX)<float> ;
} /* namespace LIBSAKURA_PREFIX */
