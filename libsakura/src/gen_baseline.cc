#include <cassert>
#include <cstdlib>
#include <climits>
#include <iostream>

#include "libsakura/localdef.h"
#include "libsakura/logger.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/sakura.h"

namespace {
auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("baseline");
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateBaselineContext)(
LIBSAKURA_SYMBOL(BaselineType) const baseline_type, uint16_t const order,
		size_t const num_data,
		LIBSAKURA_SYMBOL(BaselineContext) **context) {
	if ((baseline_type != LIBSAKURA_SYMBOL(BaselineType_kPolynomial))
			&& (baseline_type != LIBSAKURA_SYMBOL(BaselineType_kChebyshev))
			&& (baseline_type != LIBSAKURA_SYMBOL(BaselineType_kCubicSpline))
			&& (baseline_type != LIBSAKURA_SYMBOL(BaselineType_kSinusoid))) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (context == nullptr) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}

	auto baselineop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBaselineImpl();
	try {
		baselineop->CreateBaselineContext(baseline_type, order, num_data,
				context);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const std::invalid_argument &e) {
		LOG4CXX_ERROR(logger, "Order must be smaller than num_data.");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DestroyBaselineContext)(
LIBSAKURA_SYMBOL(BaselineContext) *context) {
	if (context == nullptr) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}

	auto baselineop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBaselineImpl();
	try {
		baselineop->DestroyBaselineContext(context);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetBestFitBaseline)(
		size_t num_data, float const data[], bool const mask[],
		LIBSAKURA_SYMBOL(BaselineContext) const *context, float out[]) {
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (context == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	auto baselineop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBaselineImpl();
	try {
		baselineop->GetBestFitBaseline(num_data, data, mask, context, out);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaseline)(
		size_t num_data, float const data[], bool const mask[],
		LIBSAKURA_SYMBOL(BaselineContext) const *context,
		float clipping_threshold_sigma, uint16_t num_fitting_max,
		bool get_residual, bool final_mask[], float out[]) {
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (context == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (clipping_threshold_sigma <= 0.0)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (final_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(final_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	auto baselineop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBaselineImpl();
	try {
		baselineop->SubtractBaseline(num_data, data, mask, context,
				clipping_threshold_sigma, num_fitting_max, get_residual,
				final_mask, out);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselinePolynomial)(
		size_t num_data, float const data[], bool const mask[], uint16_t order,
		float clipping_threshold_sigma, uint16_t num_fitting_max,
		bool get_residual, bool final_mask[], float out[]) {
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (order >= num_data)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (clipping_threshold_sigma <= 0.0)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (final_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(final_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	auto baselineop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBaselineImpl();
	try {
		baselineop->SubtractBaselinePolynomial(num_data, data, mask, order,
				clipping_threshold_sigma, num_fitting_max, get_residual,
				final_mask, out);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const std::invalid_argument &e) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}

	return LIBSAKURA_SYMBOL(Status_kOK);
}
