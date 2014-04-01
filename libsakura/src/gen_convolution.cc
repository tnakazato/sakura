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
	if (context == nullptr) {
		LOG4CXX_ERROR(logger, "context should not be NULL");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (!(0 < num_data && num_data <= INT_MAX)) {
		LOG4CXX_ERROR(logger, "num_data must be '0 < num_data <= INT_MAX'");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (kernel_type != LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian)
			&& kernel_type != LIBSAKURA_SYMBOL(Convolve1DKernelType_kBoxcar)
			&& kernel_type != LIBSAKURA_SYMBOL(Convolve1DKernelType_kHanning)
			&& kernel_type != LIBSAKURA_SYMBOL(Convolve1DKernelType_kHamming)) {
		LOG4CXX_ERROR(logger, "Invalid Kernel Type");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (!(0 < kernel_width)) {
		LOG4CXX_ERROR(logger, "kernel_width must be >0");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	try {
		auto convolutionop =
				::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetConvolutionImpl();
		convolutionop->CreateConvolve1DContext(num_data, kernel_type,
				kernel_width, use_fft, context);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const std::invalid_argument &e) {
		LOG4CXX_ERROR(logger, e.what());
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

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Convolve1D)(
LIBSAKURA_SYMBOL(Convolve1DContext) const *context, size_t num_data,
		float const input_data[/*num_data*/], float output_data[/*num_data*/]) {
	if (context == nullptr) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (!(0 < num_data && num_data <= INT_MAX)) {
		LOG4CXX_ERROR(logger, "num_data must be '0 < num_data <= INT_MAX'");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (!( LIBSAKURA_SYMBOL(IsAligned)(input_data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(output_data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	try {
		auto convolutionop =
				::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetConvolutionImpl();
		convolutionop->Convolve1D(context, num_data, input_data, output_data);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const std::invalid_argument &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(
LIBSAKURA_SYMBOL(Convolve1DContext) *context) {
	if (context == nullptr) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	try {
		auto convolutionop =
				::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetConvolutionImpl();
		convolutionop->DestroyConvolve1DContext(context);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}
