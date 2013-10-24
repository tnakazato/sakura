#include <iostream>
#include <cassert>
#include <cstdlib>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselinePolynomial)(
		size_t num_data, float const in_data[], bool const in_mask[],
		size_t order, float clipping_threshold_sigma,
		size_t num_fitting_max, bool get_residual, float out_data[]) {
	if (in_data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (in_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out_data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(in_data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(in_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out_data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	auto baselineop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBaselineImpl();
	try {
		baselineop->SubtractBaselinePolynomial(
				num_data, in_data, in_mask, order,
				clipping_threshold_sigma, num_fitting_max,
				get_residual, out_data);
	} catch (...) {
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetBaselineModel)(
		size_t num_data, size_t order, double out[]) {
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	auto baselineop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBaselineImpl();
	try {
		baselineop->GetBaselineModel(num_data, order, out);
	} catch (...) {
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DoSubtractBaseline)(
		size_t num_data, float const in_data[], bool const in_mask[],
		size_t num_model, double const model_data[], float clipping_threshold_sigma,
		size_t num_fitting_max, bool get_residual, float out[]) {
	if (in_data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (in_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (model_data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(in_data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(in_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(model_data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	auto baselineop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBaselineImpl();
	try {
		baselineop->DoSubtractBaseline(
			num_data, in_data, in_mask, num_model, model_data,
			clipping_threshold_sigma, num_fitting_max, get_residual, out);
	} catch (...) {
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}
