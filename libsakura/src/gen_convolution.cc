#include <cassert>
#include <cstdlib>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"
#include "libsakura/logger.h"

namespace {
auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("Convolution");
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateConvolve1DContext)(
		size_t num_data, LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
		size_t kernel_width, bool use_fft,
		LIBSAKURA_SYMBOL(Convolve1DContext) **context) {
	if (num_data < 1) {
		LOG4CXX_ERROR(logger, "ERROR: num_data must be > 0\n");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (kernel_type != LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian)
			&& kernel_type != LIBSAKURA_SYMBOL(Convolve1DKernelType_kBoxcar)
			&& kernel_type != LIBSAKURA_SYMBOL(Convolve1DKernelType_kHanning)
			&& kernel_type != LIBSAKURA_SYMBOL(Convolve1DKernelType_kHamming)) {
		LOG4CXX_ERROR(logger, "ERROR: Invalid Kernel Tyep\n");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (kernel_width < 1) {
		LOG4CXX_ERROR(logger, "ERROR: kernel_width must be > 0\n");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
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
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Convolve1D)(
LIBSAKURA_SYMBOL(Convolve1DContext) *context, size_t num_data,
	float input_data[/*num_data*/],
	bool const mask[/*num_data*/], float output_data[/*num_data*/]) {
	if (num_data < 1) {
		LOG4CXX_ERROR(logger, "ERROR: num_data must be > 0\n");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
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
