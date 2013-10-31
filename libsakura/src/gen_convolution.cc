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
	try {
		auto createconv =
				::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetConvolutionImpl();
		createconv->CreateConvolve1DContext(num_data, kernel_type,
				kernel_width, use_fft, context);
	} catch (...) {
		assert(false); // no exception should not be raised for the current implementation.
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Convolve1D)(
LIBSAKURA_SYMBOL(Convolve1DContext) **context,
        size_t num_data,
		float input_data[/*num_data*/],
		bool const input_flag[/*num_data*/],
		float output_data[/*num_data*/]) {

	try {
		auto createconv =
				::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetConvolutionImpl();
		createconv->Convolve1D(context, num_data, input_data, input_flag,
				output_data);
	} catch (...) {
		assert(false); // no exception should not be raised for the current implementation.
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(
LIBSAKURA_SYMBOL(Convolve1DContext) *context) {
	try {
		auto createconv =
				::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetConvolutionImpl();
		createconv->DestroyConvolve1DContext(context);
	} catch (...) {
		assert(false); // no exception should not be raised for the current implementation.
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}
