#include <cassert>
#include <cstdlib>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateConvolve1DContext)(
		size_t num_data, LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
		size_t kernel_width, bool use_fft,
		LIBSAKURA_SYMBOL(Convolve1DContext) **context) {
	if (num_data > 0) {
		try {
			auto convolutionop =
					::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetConvolutionImpl();
			convolutionop->CreateConvolve1DContext(num_data, kernel_type,
					kernel_width, use_fft, context);
		} catch (...) {
			assert(false); // no exception should not be raised for the current implementation.
			return LIBSAKURA_SYMBOL(Status_kUnknownError);
		}
		return LIBSAKURA_SYMBOL(Status_kOK);
	} else {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Convolve1D)(
LIBSAKURA_SYMBOL(Convolve1DContext) *context, size_t num_data,
		float input_data[/*num_data*/],
		bool const mask[/*num_data*/], float output_data[/*num_data*/]) {
	if (context != nullptr) {
		try {
			auto convolutionop =
					::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetConvolutionImpl();
			convolutionop->Convolve1D(context, num_data, input_data, mask,
					output_data);
		} catch (...) {
			assert(false); // no exception should not be raised for the current implementation.
			return LIBSAKURA_SYMBOL(Status_kUnknownError);
		}
		return LIBSAKURA_SYMBOL(Status_kOK);
	} else {
		return LIBSAKURA_SYMBOL(Status_kNG);
	}
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(
LIBSAKURA_SYMBOL(Convolve1DContext) *context) {
	try {
		auto convolutionop =
				::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetConvolutionImpl();
		convolutionop->DestroyConvolve1DContext(context);
	} catch (...) {
		assert(false); // no exception should not be raised for the current implementation.
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}
