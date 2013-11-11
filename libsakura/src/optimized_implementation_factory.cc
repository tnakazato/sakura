// Author: Kohji Nakamura <k.nakamura@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution

#include <cstdio>
#include <cassert>
#include <cstring>
#include <cstdint>

#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/logger.h"
#include "libsakura/localdef.h"

namespace {
using ::LIBSAKURA_PREFIX::ApplyCalibration;
using ::LIBSAKURA_PREFIX::Baseline;
using ::LIBSAKURA_PREFIX::BitOperation;
using ::LIBSAKURA_PREFIX::BoolFilterCollection;
using ::LIBSAKURA_PREFIX::Convolution;
using ::LIBSAKURA_PREFIX::Interpolation;
using ::LIBSAKURA_PREFIX::NumericOperation;
using ::LIBSAKURA_PREFIX::Gridding;
using ::LIBSAKURA_PREFIX::Statistics;

// a logger for this module
auto logger = LIBSAKURA_PREFIX::Logger::GetLogger(
		"optimized_implementation_factory");

struct CPURegister {
	uint32_t eax, ebx, ecx, edx;
};

void GetCpuId(CPURegister &reg) {
	__asm__ volatile("cpuid"
			: "=a" (reg.eax),
			"=b" (reg.ebx),
			"=c" (reg.ecx),
			"=d" (reg.edx)
			: "a" (reg.eax)
	);
}

struct SimdFeature {
	bool mmx, sse, sse2, sse3, ssse3, sse4_1, sse4_2, avx;
};

struct SIMD_FLAGS {
	uint32_t ecx, edx;bool SimdFeature::*flag;
} simd_flags[] = { { 0, 0x0800000, &SimdFeature::mmx }, { 0, 0x2000000,
		&SimdFeature::sse }, { 0, 0x4000000, &SimdFeature::sse2 }, { 0x00000001,
		0, &SimdFeature::sse3 }, { 0x00000200, 0, &SimdFeature::ssse3 }, {
		0x00080000, 0, &SimdFeature::sse4_1 }, { 0x00100000, 0,
		&SimdFeature::sse4_2 }, { 0x10000000, 0, &SimdFeature::avx }, };

void GetCpuFeature(SimdFeature &simd_feature) {
	CPURegister reg;
	reg.eax = 1;
	GetCpuId(reg);
#if 0
	printf("eax: %08x\n", reg.eax);
	printf("ebx: %08x\n", reg.ebx);
	printf("ecx: %08x\n", reg.ecx);
	printf("edx: %08x\n", reg.edx);
#endif

	static SimdFeature const simd_all_false = { };
	simd_feature = simd_all_false;
	for (size_t i = 0; i < ELEMENTSOF(simd_flags); ++i) {
		if ((reg.ecx & simd_flags[i].ecx) == simd_flags[i].ecx
				&& (reg.edx & simd_flags[i].edx) == simd_flags[i].edx) {
			simd_feature.*simd_flags[i].flag = true;
		}
	}
}

::LIBSAKURA_PREFIX::ApplyCalibrationDefault<float> const apply_calibration_default;
::LIBSAKURA_PREFIX::BaselineDefault const baseline_default;
::LIBSAKURA_PREFIX::BitOperationDefault<uint8_t> const bit_operation_default_uint8;
::LIBSAKURA_PREFIX::BitOperationDefault<uint32_t> const bit_operation_default_uint32;
::LIBSAKURA_PREFIX::BoolFilterCollectionDefault<float> const bool_filter_collection_default_float;
::LIBSAKURA_PREFIX::BoolFilterCollectionDefault<int> const bool_filter_collection_default_int;
::LIBSAKURA_PREFIX::BoolFilterCollectionDefault<uint8_t> const bool_filter_collection_default_uint8;
::LIBSAKURA_PREFIX::BoolFilterCollectionDefault<uint32_t> const bool_filter_collection_default_uint32;
::LIBSAKURA_PREFIX::ConvolutionDefault const convolution_default;
::LIBSAKURA_PREFIX::InterpolationDefault<double, float> const interpolation_default;
::LIBSAKURA_PREFIX::NumericOperationDefault const numeric_operation_default;
::LIBSAKURA_PREFIX::GriddingDefault const gridding_default;
::LIBSAKURA_PREFIX::StatisticsDefault const statistics_default;

class OptimizedImplementationFactoryDefault: public ::LIBSAKURA_PREFIX::OptimizedImplementationFactory {
public:
	virtual char const *GetName() const {
		return "Default";
	}
	virtual ApplyCalibration<float> const *GetApplyCalibrationImpl() const {
		return &apply_calibration_default;
	}
	virtual BitOperation<uint8_t> const *GetBitOperationImplUint8() const {
		return &bit_operation_default_uint8;
	}
	virtual BitOperation<uint32_t> const *GetBitOperationImplUint32() const {
		return &bit_operation_default_uint32;
	}
	virtual BoolFilterCollection<float> const *GetBoolFilterCollectionImplFloat() const {
		return &bool_filter_collection_default_float;
	}
	virtual BoolFilterCollection<int> const *GetBoolFilterCollectionImplInt() const {
		return &bool_filter_collection_default_int;
	}
	virtual BoolFilterCollection<uint8_t> const *GetBoolFilterCollectionImplUint8() const {
		return &bool_filter_collection_default_uint8;
	}
	virtual BoolFilterCollection<uint32_t> const *GetBoolFilterCollectionImplUint32() const {
		return &bool_filter_collection_default_uint32;
	}
	virtual Baseline const *GetBaselineImpl() const {
		return &baseline_default;
	}
	virtual Convolution const *GetConvolutionImpl() const {
		return &convolution_default;
	}
	virtual Gridding const *GetGriddingImpl() const {
		return &gridding_default;
	}
	virtual Interpolation<double, float> const *GetInterpolationImpl() const {
		return &interpolation_default;
	}
	virtual NumericOperation const *GetNumericOperationImpl() const {
		return &numeric_operation_default;
	}
	virtual Statistics const *GetStatisticsImpl() const {
		return &statistics_default;
	}

} default_factory;

::LIBSAKURA_PREFIX::ApplyCalibrationAfterSandyBridge<float> const apply_calibration_after_sandy_bridge;
::LIBSAKURA_PREFIX::BaselineAfterSandyBridge const baseline_after_sandy_bridge;
::LIBSAKURA_PREFIX::BitOperationAfterSandyBridge<uint8_t> const bit_operation_after_sandy_bridge_uint8;
::LIBSAKURA_PREFIX::BitOperationAfterSandyBridge<uint32_t> const bit_operation_after_sandy_bridge_uint32;
::LIBSAKURA_PREFIX::BoolFilterCollectionAfterSandyBridge<float> const bool_filter_collection_after_sandy_bridge_float;
::LIBSAKURA_PREFIX::BoolFilterCollectionAfterSandyBridge<int> const bool_filter_collection_after_sandy_bridge_int;
::LIBSAKURA_PREFIX::BoolFilterCollectionAfterSandyBridge<uint8_t> const bool_filter_collection_after_sandy_bridge_uint8;
::LIBSAKURA_PREFIX::BoolFilterCollectionAfterSandyBridge<uint32_t> const bool_filter_collection_after_sandy_bridge_uint32;
::LIBSAKURA_PREFIX::ConvolutionAfterSandyBridge const convolution_after_sandy_bridge;
::LIBSAKURA_PREFIX::InterpolationAfterSandyBridge<double, float> const interpolation_after_sandy_bridge;
::LIBSAKURA_PREFIX::NumericOperationAfterSandyBridge const numeric_operation_after_sandy_bridge;
::LIBSAKURA_PREFIX::GriddingAfterSandyBridge const gridding_after_sandy_bridge;
::LIBSAKURA_PREFIX::StatisticsAfterSandyBridge const statistics_after_sandy_bridge;

class OptimizedImplementationFactoryAfterSandyBridge: public ::LIBSAKURA_PREFIX::OptimizedImplementationFactory {
public:
	virtual char const *GetName() const {
		return "AfterSandyBridge";
	}
	virtual ApplyCalibration<float> const *GetApplyCalibrationImpl() const {
		return &apply_calibration_after_sandy_bridge;
	}
	virtual Baseline const *GetBaselineImpl() const {
		return &baseline_after_sandy_bridge;
	}
	virtual BitOperation<uint8_t> const *GetBitOperationImplUint8() const {
		return &bit_operation_after_sandy_bridge_uint8;
	}
	virtual BitOperation<uint32_t> const *GetBitOperationImplUint32() const {
		return &bit_operation_after_sandy_bridge_uint32;
	}
	virtual BoolFilterCollection<float> const *GetBoolFilterCollectionImplFloat() const {
		return &bool_filter_collection_after_sandy_bridge_float;
	}
	virtual BoolFilterCollection<int> const *GetBoolFilterCollectionImplInt() const {
		return &bool_filter_collection_after_sandy_bridge_int;
	}
	virtual BoolFilterCollection<uint8_t> const *GetBoolFilterCollectionImplUint8() const {
		return &bool_filter_collection_after_sandy_bridge_uint8;
	}
	virtual BoolFilterCollection<uint32_t> const *GetBoolFilterCollectionImplUint32() const {
		return &bool_filter_collection_after_sandy_bridge_uint32;
	}
	virtual Convolution const *GetConvolutionImpl() const {
		// return &convolution_after_sandy_bridge;
		return &convolution_after_sandy_bridge;
	}
	virtual Gridding const *GetGriddingImpl() const {
		return &gridding_after_sandy_bridge;
	}
	virtual Interpolation<double, float> const *GetInterpolationImpl() const {
		return &interpolation_after_sandy_bridge;
//		return &interpolation_default;
	}
	virtual NumericOperation const *GetNumericOperationImpl() const {
		return &numeric_operation_after_sandy_bridge;
	}
	virtual Statistics const *GetStatisticsImpl() const {
		return &statistics_after_sandy_bridge;
	}
} after_sandy_bridge;

}

namespace LIBSAKURA_PREFIX {
OptimizedImplementationFactory const *OptimizedImplementationFactory::factory_ =
		nullptr;

void OptimizedImplementationFactory::InitializeFactory(char const *simd_spec) {
	SimdFeature simd_feature;
	GetCpuFeature(simd_feature);

	if (strcmp(simd_spec, "adaptive") == 0) {
		if (simd_feature.avx) {
			factory_ = &after_sandy_bridge;
		} else {
			factory_ = &default_factory;
		}
	} else if (simd_feature.avx && strcmp(simd_spec, "avx") == 0) {
		factory_ = &after_sandy_bridge;
	} else if (strcmp(simd_spec, "disabled") == 0) {
		factory_ = &default_factory;
	} else {
		factory_ = &default_factory;
	}
	if (::LIBSAKURA_PREFIX::Logger::IsDebugEnabled(logger)) {
		std::ostringstream os;
		os << "SIMD implementation: " << factory_->GetName() << std::endl;
		Logger::Debug(logger, os.str().c_str());
	}
}

void OptimizedImplementationFactory::CleanUpFactory() {
	factory_ = nullptr;
}

OptimizedImplementationFactory const *
OptimizedImplementationFactory::GetFactory() {
	assert("Call LIBSAKURA_SYMBOL(Initialize)() first." && factory_ != nullptr);
	return factory_;
}
}
