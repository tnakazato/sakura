#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/localdef.h"

#define FORCE_EIGEN 0

#if defined(__AVX__) && (! FORCE_EIGEN)
#include <immintrin.h>
#include <cstdint>

namespace {

void CreateConvole1DContextSimd(size_t num_channel,LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
        size_t kernel_width,bool use_fft,LIBSAKURA_SYMBOL(Convole1DContext) **context) {
	std::cout << "This function is not implemented yet." << std::endl;
}

} /* anonymous namespace */

#else /* defined(__AVX__) */

#define EIGEN_DENSEBASE_PLUGIN "eigen_binary_visitor_plugin.h"
#include <Eigen/Core>

using ::Eigen::Map;
using ::Eigen::Array;
using ::Eigen::Dynamic;
using ::Eigen::Aligned;

namespace {

inline void CreateConvole1DContextEigen(size_t num_channel,LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
        size_t kernel_width,bool use_fft,LIBSAKURA_SYMBOL(Convole1DContext) **context) {
	std::cout << " CreateConvole1DContextEigen function is called" << std::endl;

}

} /* anonymous namespace */

#endif /* defined(__AVX__) */

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(Convolution, ARCH_SUFFIX)::CreateConvole1DContext(size_t num_channel,LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
        size_t kernel_width,bool use_fft,LIBSAKURA_SYMBOL(Convole1DContext) **context) const {
#if defined( __AVX__) && (! FORCE_EIGEN)
	CreateConvole1DContextSimd(num_channel,kernel_type,kernel_width,use_fft,context);
#else
	CreateConvole1DContextEigen(num_channel,kernel_type,kernel_width,use_fft,context);
#endif
}

}
