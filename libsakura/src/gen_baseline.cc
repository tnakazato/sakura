#include <iostream>
#include <cassert>
#include <cstdlib>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselinePolynomial)(
		size_t num_chan, float const in_data[], bool const in_mask[], int order,
		float clipping_threshold_sigma, unsigned int num_fitting_max,
		bool get_residual, float out_data[]) {

	auto baselineop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBaselineImpl();

	baselineop->SubtractBaselinePolynomial(num_chan, in_data, in_mask, order,
			clipping_threshold_sigma, num_fitting_max, get_residual,
			out_data);
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetBaselineModel)(
		size_t num_chan, int order, double out[]) {

	auto baselineop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBaselineImpl();

	baselineop->GetBaselineModel(num_chan, order, out);
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DoSubtractBaseline)(
		size_t num_chan, float const in_data[], bool const in_mask[],
		size_t num_model, double model_data[], float clipping_threshold_sigma,
		unsigned int num_fitting_max, bool get_residual, float out[]) {

	auto baselineop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBaselineImpl();

	baselineop->DoSubtractBaseline(
			num_chan, in_data, in_mask, num_model, model_data,
			clipping_threshold_sigma, num_fitting_max, get_residual, out);
	return LIBSAKURA_SYMBOL(Status_kOK);
}
