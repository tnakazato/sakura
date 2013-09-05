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
struct LIBSAKURA_SYMBOL(Convolve1DContext) {
	fftwf_plan plan_real_to_complex_float;
	fftwf_plan plan_complex_to_real_float;
	fftwf_complex *fft_applied_complex_kernel;
	float input_real_kernel[0];
};
}

#if defined(__AVX__) && (! FORCE_EIGEN)
#include <immintrin.h>
#include <cstdint>

namespace {

void CreateConvolve1DContextSimd(size_t num_channel,
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type, size_t kernel_width,
		bool use_fft, LIBSAKURA_SYMBOL(Convolve1DContext) **context) {
	std::cout << "This function implementation is in progress." << std::endl;
}
void Convolve1DSimd(LIBSAKURA_SYMBOL(Convolve1DContext) *context,
		size_t num_channel,
		float const input_spectrum[/*num_channels*/],bool const input_flag[/*num_channels*/],
		float output_spectrum[/*num_channels*/]) {
	std::cout << "This function implementation is in progress." << std::endl;
}
void DestroyConvolve1DContextSimd(
		LIBSAKURA_SYMBOL(Convolve1DContext) *context) {
	std::cout << "This function implementation is in progress." << std::endl;
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
			size_t kernel_width,bool use_fft,LIBSAKURA_SYMBOL(Convolve1DContext)** context) {
		std::cout << " CreateConvolve1DContextEigen function is called" << std::endl;

		const float ln2 = 0.6931471805599453094172321;
		const float pi = 3.141592653589793238462643;
		const float sigma = kernel_width / sqrt(float(8.0) * ln2);
		const float center = float(num_channel)/2;
		float height = 1.0 / (sigma * sqrt(2.0 * pi));
		float fwhm2int(float(1.0)/sqrt(log(float(16.0))));

		(*context)=(LIBSAKURA_SYMBOL(Convolve1DContext) *)malloc(sizeof(struct LIBSAKURA_SYMBOL(Convolve1DContext)) + sizeof(float)*num_channel);
        if(*context == NULL){
        	std::cout << " context wasn't created properly " << std::endl;
        	exit(1);
        }

		for (uint j=0; j<num_channel; ++j) {
			float value = (j - center)/kernel_width/fwhm2int;
			(*context)->input_real_kernel[j] = height * exp(-(value*value));
		}

		(*context)->fft_applied_complex_kernel = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * num_channel );
		(*context)->plan_real_to_complex_float = fftwf_plan_dft_r2c_1d(num_channel, (*context)->input_real_kernel, (*context)->fft_applied_complex_kernel,FFTW_ESTIMATE);

		fftwf_execute( (*context)-> plan_real_to_complex_float);
	}

	inline void Convolve1DEigen(LIBSAKURA_SYMBOL(Convolve1DContext) *context,
				     size_t num_channel,float input_spectrum[/*num_channels*/],bool const input_flag[/*num_channels*/],
				     float output_spectrum[/*num_channels*/]) {

		fftwf_complex *input_complex_spectrum;
	    fftwf_complex *output_complex_spectrum;
		fftwf_plan plan_r2cf_input_spectrum;
		fftwf_plan plan_c2rf_output_spectrum;

		input_complex_spectrum  = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (num_channel/2 + 1) );
		output_complex_spectrum = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (num_channel/2 + 1) );

		if(input_complex_spectrum == NULL ||  output_complex_spectrum == NULL){
        	std::cout << " input/output _complex_spectrum weren't created properly " << std::endl;
        	exit(1);
		}

		plan_r2cf_input_spectrum = fftwf_plan_dft_r2c_1d(num_channel, input_spectrum, input_complex_spectrum,FFTW_ESTIMATE);
		fftwf_execute(plan_r2cf_input_spectrum);

		for(uint i=0; i < num_channel/2 + 1; ++i){
		      output_complex_spectrum[i][0] = context->fft_applied_complex_kernel[i][0] * input_complex_spectrum[i][0]
		                                    - context->fft_applied_complex_kernel[i][1] * input_complex_spectrum[i][1];
		      output_complex_spectrum[i][1] = context->fft_applied_complex_kernel[i][0] * input_complex_spectrum[i][1]
		                                    + context->fft_applied_complex_kernel[i][1] * input_complex_spectrum[i][0];
		}

		plan_c2rf_output_spectrum = fftwf_plan_dft_c2r_1d(num_channel, output_complex_spectrum , output_spectrum, FFTW_ESTIMATE);
		fftwf_execute(plan_c2rf_output_spectrum);

		fftwf_destroy_plan(plan_r2cf_input_spectrum );
		fftwf_destroy_plan(plan_c2rf_output_spectrum );
		fftwf_free(input_complex_spectrum);
	    fftwf_free(output_complex_spectrum);
	}

	inline void DestroyConvolve1DContextEigen(LIBSAKURA_SYMBOL(Convolve1DContext)* context) {
		fftwf_destroy_plan( (context)->plan_real_to_complex_float);
		fftwf_destroy_plan( (context)->plan_complex_to_real_float);
		fftwf_free( (context)->fft_applied_complex_kernel);
		free(context);
	}

} /* anonymous namespace */

#endif /* defined(__AVX__) */

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(Convolution, ARCH_SUFFIX)::CreateConvolve1DContext(
		size_t num_channel, LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
		size_t kernel_width, bool use_fft,
		LIBSAKURA_SYMBOL(Convolve1DContext) **context) const {
#if defined( __AVX__) && (! FORCE_EIGEN)
	CreateConvolve1DContextSimd(num_channel, kernel_type, kernel_width, use_fft,
			context);
#else
	CreateConvolve1DContextEigen(num_channel,kernel_type,kernel_width,use_fft,context);
#endif
}

void ADDSUFFIX(Convolution, ARCH_SUFFIX)::Convolve1D(
		LIBSAKURA_SYMBOL(Convolve1DContext) *context,
		size_t num_channel, float input_spectrum[/*num_channels*/],bool const input_flag[/*num_channels*/],
		float output_spectrum[/*num_channels*/]) const {
#if defined( __AVX__) && (! FORCE_EIGEN)
	Convolve1DSimd(context,num_channel,input_spectrum,input_flag,output_spectrum);
#else
	Convolve1D(context,num_channel,input_spectrum,input_flag,output_spectrum);
#endif
}

void ADDSUFFIX(Convolution, ARCH_SUFFIX)::DestroyConvolve1DContext(
		LIBSAKURA_SYMBOL(Convolve1DContext) *context) const {
#if defined( __AVX__) && (! FORCE_EIGEN)
	DestroyConvolve1DContextSimd(context);
#else
	DestroyConvolve1DContextEigen(context);
#endif
}

}
