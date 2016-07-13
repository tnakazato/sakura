/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2016
 * National Astronomical Observatory of Japan
 * 2-21-1, Osawa, Mitaka, Tokyo, 181-8588, Japan.
 * 
 * This file is part of Sakura.
 * 
 * Sakura is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the 
 * Free Software Foundation, either version 3 of the License, or (at your 
 * option) any later version.
 * 
 * Sakura is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License 
 * along with Sakura.  If not, see <http://www.gnu.org/licenses/>.
 * @SAKURA_LICENSE_HEADER_END@
 */

#include <cassert>

#include <Eigen/Core>

#include "libsakura/localdef.h"
#include "libsakura/logger.h"
#include "libsakura/packed_type.h"
#include "libsakura/sakura.h"

using ::Eigen::Map;
using ::Eigen::Array;
using ::Eigen::Dynamic;
using ::Eigen::Aligned;

namespace {

// a logger for this module
auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("normalization");

template<class DataType>
struct InPlaceImpl {
	static void CalibrateData(size_t num_data,
			DataType const scaling_factor[/*num_data*/],
			DataType const reference[/*num_data*/],
			DataType result[/*num_data*/]) {
		Map<Array<DataType, Dynamic, 1>, Aligned> eigen_result(result,
				num_data);
		Map<Array<DataType, Dynamic, 1>, Aligned> eigen_reference(
				const_cast<DataType *>(reference), num_data);
		Map<Array<DataType, Dynamic, 1>, Aligned> eigen_scaling_factor(
				const_cast<DataType *>(scaling_factor), num_data);
		eigen_result = eigen_scaling_factor * (eigen_result - eigen_reference)
				/ eigen_reference;
	}
	static void CalibrateData(DataType scaling_factor, size_t num_data,
			DataType const reference[/*num_data*/],
			DataType result[/*num_data*/]) {
		Map<Array<DataType, Dynamic, 1>, Aligned> eigen_result(result,
				num_data);
		Map<Array<DataType, Dynamic, 1>, Aligned> eigen_reference(
				const_cast<DataType *>(reference), num_data);
		DataType const constant_scaling_factor = scaling_factor;
		eigen_result = constant_scaling_factor
				* (eigen_result - eigen_reference) / eigen_reference;
	}
};

#if defined(__AVX__) && !defined(ARCH_SCALAR)
#include <immintrin.h>

template<>
struct InPlaceImpl<float> {
	typedef __m256 SimdType;
	static void CalibrateData(size_t num_data,
			float const scaling_factor[/*num_data*/],
			float const reference[/*num_data*/], float result[/*num_data*/]) {
		constexpr size_t kNumFloat = LIBSAKURA_SYMBOL(SimdPacketAVX)::kNumFloat;
		size_t num_packed_operation = num_data / kNumFloat;
		size_t num_data_packed = num_packed_operation * kNumFloat;
		assert(num_data_packed <= num_data);
		size_t num_extra = num_data - num_data_packed;
		auto const packed_reference =
				reinterpret_cast<SimdType const *>(reference);
		auto packed_result = reinterpret_cast<SimdType *>(result);
		auto const reference_extra = &reference[num_data_packed];
		auto result_extra = &result[num_data_packed];
		auto const packed_factor =
				reinterpret_cast<SimdType const *>(scaling_factor);
		IterateSimd([packed_factor] (size_t i) {return packed_factor[i];},
				num_packed_operation, packed_reference, packed_result);
		auto const factor_extra = &scaling_factor[num_data_packed];
		IterateExtra([factor_extra] (size_t i) {return factor_extra[i];},
				num_extra, reference_extra, result_extra);
	}
	static void CalibrateData(float scaling_factor, size_t num_data,
			float const reference[/*num_data*/], float result[/*num_data*/]) {
		constexpr size_t kNumFloat = LIBSAKURA_SYMBOL(SimdPacketAVX)::kNumFloat;
		size_t num_packed_operation = num_data / kNumFloat;
		size_t num_data_packed = num_packed_operation * kNumFloat;
		assert(num_data_packed <= num_data);
		size_t num_extra = num_data - num_data_packed;
		auto const packed_reference =
				reinterpret_cast<SimdType const *>(reference);
		auto packed_result = reinterpret_cast<SimdType *>(result);
		auto const reference_extra = &reference[num_data_packed];
		auto result_extra = &result[num_data_packed];
		SimdType const packed_scalar_factor = _mm256_broadcast_ss(
				&scaling_factor);
		IterateSimd(
				[packed_scalar_factor] (size_t i) {return packed_scalar_factor;},
				num_packed_operation, packed_reference, packed_result);
		IterateExtra([scaling_factor] (size_t i) {return scaling_factor;},
				num_extra, reference_extra, result_extra);
	}
private:
	template<class Feeder>
	static void IterateExtra(Feeder factor_feeder, size_t num_data,
			float const *reference, float *result) {
		for (size_t i = 0; i < num_data; ++i) {
			result[i] = factor_feeder(i) * (result[i] - reference[i])
					/ reference[i];
		}
	}

	template<class Feeder>
	static void IterateSimd(Feeder factor_feeder, size_t num_data,
			SimdType const *packed_reference, SimdType *packed_result) {
		for (size_t i = 0; i < num_data; ++i) {
			// Here, we don't use _mm256_rcp_ps with _mm256_mul_ps instead of
			// _mm256_div_ps since the former loses accuracy (worse than
			// documented).
			packed_result[i] = _mm256_div_ps(
					_mm256_mul_ps(factor_feeder(i),
							_mm256_sub_ps(packed_result[i],
									packed_reference[i])), packed_reference[i]);
		}
	}
};
#endif

template<class DataType>
struct DefaultImpl {
	static void CalibrateData(size_t num_data, DataType const *scaling_factor,
			DataType const *target, DataType const *reference,
			DataType *result) {
		Iterate([scaling_factor] (size_t i) {return scaling_factor[i];},
				num_data, target, reference, result);
	}
	static void CalibrateData(DataType scaling_factor, size_t num_data,
			DataType const *target, DataType const *reference,
			DataType *result) {
		Iterate([scaling_factor] (size_t i) {return scaling_factor;}, num_data,
				target, reference, result);
	}
private:
	template<class Feeder>
	static void Iterate(Feeder factor_feeder, size_t num_data,
			DataType const *target_arg, DataType const *reference_arg,
			DataType *result_arg) {
		DataType const *__restrict target = AssumeAligned(target_arg);
		DataType const *__restrict reference = AssumeAligned(reference_arg);
		DataType *__restrict result = AssumeAligned(result_arg);
		for (size_t i = 0; i < num_data; ++i) {
			result[i] = factor_feeder(i) * (target[i] - reference[i])
					/ reference[i];
		}
	}
};

template<class DataType>
void CalibrateData(size_t num_data, DataType const scaling_factor[/*num_data*/],
		DataType const target[/*num_data*/],
		DataType const reference[/*num_data*/], DataType result[/*num_data*/]) {
	assert(scaling_factor != nullptr);
	assert(target != nullptr);
	assert(reference != nullptr);
	assert(result != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(scaling_factor));
	assert(LIBSAKURA_SYMBOL(IsAligned)(target));
	assert(LIBSAKURA_SYMBOL(IsAligned)(reference));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	if (target == result) {
		InPlaceImpl<DataType>::CalibrateData(num_data, scaling_factor,
				reference, result);
	} else {
		DefaultImpl<DataType>::CalibrateData(num_data, scaling_factor, target,
				reference, result);
	}
}

template<class DataType>
void CalibrateData(DataType scaling_factor, size_t num_data,
		DataType const target[/*num_data*/],
		DataType const reference[/*num_data*/], DataType result[/*num_data*/]) {
	assert(target != nullptr);
	assert(reference != nullptr);
	assert(result != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(target));
	assert(LIBSAKURA_SYMBOL(IsAligned)(reference));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	if (target == result) {
		InPlaceImpl<DataType>::CalibrateData(scaling_factor, num_data,
				reference, result);
	} else {
		DefaultImpl<DataType>::CalibrateData(scaling_factor, num_data, target,
				reference, result);
	}
}

} /* anonymous namespace */

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CalibrateDataWithArrayScalingFloat)(
		size_t num_data, float const scaling_factor[/*num_data*/],
		float const target[/*num_data*/], float const reference[/*num_data*/],
		float result[/*num_data*/]) noexcept {
	if (num_data == 0) {
		// Nothing to do
		return LIBSAKURA_SYMBOL(Status_kOK);
	}

	if (scaling_factor == nullptr|| target == nullptr || reference == nullptr || result == nullptr) {
		// null pointer
		LOG4CXX_ERROR(logger, "Input pointers are null");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}

	if (!LIBSAKURA_SYMBOL(IsAligned)(scaling_factor)
			|| !LIBSAKURA_SYMBOL(IsAligned)(target)
			|| !LIBSAKURA_SYMBOL(IsAligned)(reference)
			|| !LIBSAKURA_SYMBOL(IsAligned)(result)) {
		// array is not aligned
		LOG4CXX_ERROR(logger, "Arrays are not aligned");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}

	try {
		CalibrateData(num_data, scaling_factor, target, reference, result);
		return LIBSAKURA_SYMBOL(Status_kOK);
	} catch (...) {
		// any exception is thrown during interpolation
		assert(false);
		LOG4CXX_ERROR(logger, "Aborted due to unknown error");
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CalibrateDataWithConstScalingFloat)(
		float scaling_factor, size_t num_data, float const target[/*num_data*/],
		float const reference[/*num_data*/], float result[/*num_data*/])
				noexcept {
	if (num_data == 0) {
		// Nothing to do
		return LIBSAKURA_SYMBOL(Status_kOK);
	}

	if (target == nullptr|| reference == nullptr || result == nullptr) {
		// null pointer
		LOG4CXX_ERROR(logger, "Input pointers are null");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}

	if (!LIBSAKURA_SYMBOL(IsAligned)(target)
			|| !LIBSAKURA_SYMBOL(IsAligned)(reference)
			|| !LIBSAKURA_SYMBOL(IsAligned)(result)) {
		// array is not aligned
		LOG4CXX_ERROR(logger, "Arrays are not aligned");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}

	try {
		CalibrateData(scaling_factor, num_data, target, reference, result);
		return LIBSAKURA_SYMBOL(Status_kOK);
	} catch (...) {
		// any exception is thrown during interpolation
		assert(false);
		LOG4CXX_ERROR(logger, "Aborted due to unknown error");
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
}
