/*
 * @SAKURA_LICENSE_HEADER_START@
 * @SAKURA_LICENSE_HEADER_END@
 */
#include <cassert>
#include <cstdlib>
#include <climits>
#include <iostream>

#include "libsakura/localdef.h"
#include "libsakura/logger.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/sakura.h"
#include "baseline.h"

namespace {
auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("baseline");
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateBaselineContext)(
LIBSAKURA_SYMBOL(BaselineType) const baseline_type, uint16_t const order,
		size_t const num_data,
		LIBSAKURA_SYMBOL(BaselineContext) **context) {
	if (baseline_type >= LIBSAKURA_SYMBOL(BaselineType_kNumType)) {
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
	} catch (const std::runtime_error &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kNG);
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
		LIBSAKURA_SYMBOL(BaselineContext) const *context, float out[],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
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
	if (num_data != context->num_basis_data)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_data < context->num_bases)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (baseline_status == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	auto baselineop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBaselineImpl();
	try {
		baselineop->GetBestFitBaseline(num_data, data, mask, context, out,
				baseline_status);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const std::runtime_error &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kNG);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaseline)(
		size_t num_data, float const data[], bool const mask[],
		LIBSAKURA_SYMBOL(BaselineContext) const *context,
		float clip_threshold_sigma, uint16_t num_fitting_max,
		bool get_residual, bool final_mask[], float out[],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
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
	if (num_data != context->num_basis_data)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_data < context->num_bases)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (clip_threshold_sigma <= 0.0f)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (final_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(final_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (baseline_status == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	auto baselineop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBaselineImpl();
	try {
		baselineop->SubtractBaseline(num_data, data, mask, context,
				clip_threshold_sigma, num_fitting_max, get_residual, final_mask,
				out, baseline_status);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const std::runtime_error &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kNG);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselinePolynomial)(
		size_t num_data, float const data[], bool const mask[], uint16_t order,
		float clip_threshold_sigma, uint16_t num_fitting_max,
		bool get_residual, bool final_mask[], float out[],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	if (num_data <= order)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (clip_threshold_sigma <= 0.0f)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (final_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(final_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (baseline_status == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	auto baselineop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBaselineImpl();
	try {
		baselineop->SubtractBaselinePolynomial(num_data, data, mask, order,
				clip_threshold_sigma, num_fitting_max, get_residual, final_mask,
				out, baseline_status);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const std::runtime_error &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kNG);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}
