#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fftw3.h>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/localdef.h"

#define FORCE_EIGEN 0

extern "C" {
struct LIBSAKURA_SYMBOL(Convole1DContext) {
	fftwf_plan plan_real_to_complex_float;
	fftwf_plan plan_complex_to_real_float;
	fftwf_complex *complex_array;
	float fft_applied_kernel[0];
};
}

#if defined(__AVX__) && (! FORCE_EIGEN)
#include <immintrin.h>
#include <cstdint>

namespace {

void CreateConvolve1DContextSimd(size_t num_channel,
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type, size_t kernel_width,
		bool use_fft, LIBSAKURA_SYMBOL(Convole1DContext) **context) {
	std::cout << "This function is not implemented yet." << std::endl;
}
void DestroyConvolve1DContextSimd(
		LIBSAKURA_SYMBOL(Convole1DContext) **context) {
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
			size_t kernel_width,bool use_fft,LIBSAKURA_SYMBOL(Convole1DContext)** context) {
		std::cout << " CreateConvole1DContextEigen function is called" << std::endl;

		const float ln2 = 0.6931471805599453094172321;
		const float pi = 3.141592653589793238462643;
		const float sigma = kernel_width / sqrt(float(8.0) * ln2);
		const float center = float(num_channel)/2;
		float height = 1.0 / (sigma * sqrt(2.0 * pi));
		float fwhm2int(float(1.0)/sqrt(log(float(16.0))));

		(*context)=(LIBSAKURA_SYMBOL(Convole1DContext) *)malloc(sizeof(struct LIBSAKURA_SYMBOL(Convole1DContext)) + sizeof(float)*num_channel);

		for (uint j=0; j<num_channel; j++) {
			float value = (j - center)/kernel_width/fwhm2int;
			(*context)->fft_applied_kernel[j] = height * exp(-(value*value));
		}

		(*context)->complex_array = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * num_channel );
		(*context)->plan_real_to_complex_float = fftwf_plan_dft_r2c_1d(num_channel/2 + 1, (*context)->fft_applied_kernel, (*context)->complex_array,FFTW_MEASURE);
		//(*context)->plan_complex_to_real_float = fftwf_plan_dft_c2r_1d(num_channel, (*f)->complex_array,(*f)->fft_applied_kernel,FFTW_MEASURE);

		fftwf_execute( (*context)-> plan_real_to_complex_float);
	}

	inline void DestroyConvolve1DContextEigen(LIBSAKURA_SYMBOL(Convole1DContext)** context) {
		fftwf_destroy_plan( (*context)->plan_real_to_complex_float);
		fftwf_destroy_plan( (*context)->plan_complex_to_real_float);
		fftwf_free( (*context)->complex_array);
		free(*context);
	}

} /* anonymous namespace */

#endif /* defined(__AVX__) */

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(Convolution, ARCH_SUFFIX)::CreateConvolve1DContext(
		size_t num_channel, LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
		size_t kernel_width, bool use_fft,
		LIBSAKURA_SYMBOL(Convole1DContext) **context) const {
#if defined( __AVX__) && (! FORCE_EIGEN)
	CreateConvolve1DContextSimd(num_channel, kernel_type, kernel_width, use_fft,
			context);
#else
	CreateConvolve1DContextEigen(num_channel,kernel_type,kernel_width,use_fft,context);
#endif
}

void ADDSUFFIX(Convolution, ARCH_SUFFIX)::DestroyConvolve1DContext(
		LIBSAKURA_SYMBOL(Convole1DContext) **context) const {
#if defined( __AVX__) && (! FORCE_EIGEN)
	DestroyConvolve1DContextSimd(context);
#else
	DestroyConvolve1DContextEigen(context);
#endif
}

}
