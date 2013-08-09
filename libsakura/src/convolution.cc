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

void CreateConvolve1DContextSimd(size_t num_channel,LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
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

inline void CreateConvolve1DContextEigen(size_t num_channel,LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
        size_t kernel_width,bool use_fft,LIBSAKURA_SYMBOL(Convole1DContext) **context) {
	std::cout << " CreateConvole1DContextEigen function is called" << std::endl;

    //const double ln2 = 0.6931471805599453094172321;
    //const double pi = 3.141592653589793238462643;
    //const double sigma = width / sqrt(double(8.0) * ln2);
    //const double refPix = double(nPixels)/2;
    //double norm;
	// norm = 1.0 / (sigma * sqrt(2.0 * pi));
    //const Gaussian1D<double> gauss(norm, refPix, double(width));
    //for (uint j=0; j<nPixels; j++) kernel[j] = gauss(double(j));
}

} /* anonymous namespace */

#endif /* defined(__AVX__) */

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(Convolution, ARCH_SUFFIX)::CreateConvolve1DContext(size_t num_channel,LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
        size_t kernel_width,bool use_fft,LIBSAKURA_SYMBOL(Convole1DContext) **context) const {
#if defined( __AVX__) && (! FORCE_EIGEN)
	CreateConvolve1DContextSimd(num_channel,kernel_type,kernel_width,use_fft,context);
#else
	CreateConvolve1DContextEigen(num_channel,kernel_type,kernel_width,use_fft,context);
#endif
}

}
