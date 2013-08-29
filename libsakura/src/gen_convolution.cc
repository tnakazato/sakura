#include <cassert>
#include <cstdlib>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateConvolve1DContext)(
		             size_t num_channel,LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
		             size_t kernel_width,bool use_fft,LIBSAKURA_SYMBOL(Convole1DContext) **context){

	auto createconv =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetConvolutionImpl();
	createconv->CreateConvolve1DContext(num_channel,kernel_type,kernel_width,use_fft,context);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(
		             LIBSAKURA_SYMBOL(Convole1DContext) **context){

	auto createconv =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetConvolutionImpl();
	createconv->DestroyConvolve1DContext(context);

	return LIBSAKURA_SYMBOL(Status_kOK);
}
