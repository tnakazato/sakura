#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fftw3.h>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/localdef.h"

extern "C" {
struct LIBSAKURA_SYMBOL(Convolve1DContext) {
	fftwf_plan plan_real_to_complex_float;
	fftwf_plan plan_complex_to_real_float;
	fftwf_complex *input_complex_spectrum;
	fftwf_complex *output_complex_spectrum;
	fftwf_complex *fft_applied_complex_kernel;
	size_t num_channel;
	float input_real_array[0];
};
}

namespace {

inline void CreateConvolve1DContextDefault(size_t num_channel,
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type, size_t kernel_width,
		bool use_fft, LIBSAKURA_SYMBOL(Convolve1DContext)** context) {
	std::cout << " CreateConvolve1DContextEigen function is called"
			<< std::endl;

	float const sqrt_ln2_mul8_over_2pi = .939437278699651333772340328410;
	float const sqrt_ln16 = 1.66510922231539551270632928979040;
	float const height = sqrt_ln2_mul8_over_2pi / kernel_width;
	float center = num_channel / 2.f;

	switch (kernel_type) {
	case LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian):
		std::cout << "Kernel : Gaussian" << std::endl;

		(*context) = (LIBSAKURA_SYMBOL(Convolve1DContext) *) malloc(
				sizeof(struct LIBSAKURA_SYMBOL(Convolve1DContext))
						+ sizeof(float) * num_channel);
		assert(LIBSAKURA_SYMBOL(IsAligned)(*context));
		if (*context == NULL) {
			std::cout << " context memory wasn't allocated properly "
					<< std::endl;
			exit(1);
		}

		(*context)->num_channel = num_channel;

		for (size_t j = 0; j < num_channel; ++j) {
			float value = (j - center) * sqrt_ln16 / kernel_width;
			(*context)->input_real_array[j] = height * exp(-(value * value));
		}
		assert(LIBSAKURA_SYMBOL(IsAligned)((*context)->input_real_array));

		if (use_fft) {
			(*context)->fft_applied_complex_kernel =
					(fftwf_complex*) fftwf_malloc(
							sizeof(fftwf_complex) * (num_channel / 2 + 1));
			(*context)->input_complex_spectrum = (fftwf_complex*) fftwf_malloc(
					sizeof(fftwf_complex) * (num_channel / 2 + 1));
			(*context)->output_complex_spectrum = (fftwf_complex*) fftwf_malloc(
					sizeof(fftwf_complex) * (num_channel / 2 + 1));
			if ((*context)->fft_applied_complex_kernel == NULL
					|| (*context)->input_complex_spectrum == NULL
					|| (*context)->output_complex_spectrum == NULL) {
				std::cout
						<< " context->fft_applied_complex_kernel or (*context)->input_complex_spectrum memory wasn't allocated properly "
						<< std::endl;
				exit(1);
			}
			assert(
					LIBSAKURA_SYMBOL(IsAligned)((*context)->fft_applied_complex_kernel));
			assert(
					LIBSAKURA_SYMBOL(IsAligned)((*context)->input_complex_spectrum));
			assert(
					LIBSAKURA_SYMBOL(IsAligned)((*context)->output_complex_spectrum));

			(*context)->plan_real_to_complex_float = fftwf_plan_dft_r2c_1d(
					num_channel, (*context)->input_real_array,
					(*context)->fft_applied_complex_kernel, FFTW_ESTIMATE);
			if ((*context)->plan_real_to_complex_float != NULL)
				fftwf_execute((*context)->plan_real_to_complex_float);

			(*context)->plan_real_to_complex_float = fftwf_plan_dft_r2c_1d(
					num_channel, (*context)->input_real_array,
					(*context)->input_complex_spectrum, FFTW_ESTIMATE);
			(*context)->plan_complex_to_real_float = fftwf_plan_dft_c2r_1d(
					num_channel, (*context)->output_complex_spectrum,
					(*context)->input_real_array, FFTW_ESTIMATE);
			if ((*context)->plan_real_to_complex_float == NULL
					|| (*context)->plan_complex_to_real_float == NULL) {
				std::cout
						<< "(*context)->plan_real_to_complex_float or (*context)->plan_complex_to_real_float wasn't created properly "
						<< std::endl;
				exit(1);
			}
		} else
			std::cout << "with out fft; not implemented yet " << std::endl;
		break;
	case LIBSAKURA_SYMBOL(Convolve1DKernelType_kBoxcar):
		std::cout << "Kernel : Boxcar" << std::endl;
		break;
	case LIBSAKURA_SYMBOL(Convolve1DKernelType_kHanning):
		std::cout << "Kernel : Hanning" << std::endl;
		break;
	case LIBSAKURA_SYMBOL(Convolve1DKernelType_kHamming):
		std::cout << "Kernel : Hamming" << std::endl;
		break;
	default:
		std::cout << "Kernel : no" << std::endl;
		break;
	}
}

inline void Convolve1DDefault(LIBSAKURA_SYMBOL(Convolve1DContext) **context,
		float input_spectrum[/*num_channels*/],
		bool const input_flag[/*num_channels*/],
		float output_spectrum[/*num_channels*/]) {
	std::cout << " Convolve1DEigen function is called" << std::endl;

	assert(LIBSAKURA_SYMBOL(IsAligned)(input_spectrum));
	assert(LIBSAKURA_SYMBOL(IsAligned)(output_spectrum));

	uint i;
	for (i = 0; i < (*context)->num_channel; ++i) {
		(*context)->input_real_array[i] = input_spectrum[i];
	}

	if ((*context)->plan_real_to_complex_float != NULL)
		fftwf_execute((*context)->plan_real_to_complex_float);

	float scale = 1.0 / (*context)->num_channel;
	for (uint i = 0; i < (*context)->num_channel / 2 + 1; ++i) {
		(*context)->output_complex_spectrum[i][0] =
				((*context)->fft_applied_complex_kernel[i][0]
						* (*context)->input_complex_spectrum[i][0]
						- (*context)->fft_applied_complex_kernel[i][1]
								* (*context)->input_complex_spectrum[i][1])
						* scale;
		(*context)->output_complex_spectrum[i][1] =
				((*context)->fft_applied_complex_kernel[i][0]
						* (*context)->input_complex_spectrum[i][1]
						+ (*context)->fft_applied_complex_kernel[i][1]
								* (*context)->input_complex_spectrum[i][0])
						* scale;
	}

	if ((*context)->plan_complex_to_real_float != NULL)
		fftwf_execute((*context)->plan_complex_to_real_float);

	for (i = 0; i < (*context)->num_channel; ++i) {
		output_spectrum[i] = (*context)->input_real_array[i];
	}
}

inline void DestroyConvolve1DContextDefault(
		LIBSAKURA_SYMBOL(Convolve1DContext)* context) {
	std::cout << " DestroyConvolve1DContextEigen function is called"
			<< std::endl;
	if (context->plan_real_to_complex_float != NULL) {
		fftwf_destroy_plan(context->plan_real_to_complex_float);
		context->plan_real_to_complex_float = NULL;
	}
	if (context->plan_complex_to_real_float != NULL) {
		fftwf_destroy_plan(context->plan_complex_to_real_float);
		context->plan_complex_to_real_float = NULL;
	}
	if (context->fft_applied_complex_kernel != NULL) {
		fftwf_free(context->fft_applied_complex_kernel);
		context->fft_applied_complex_kernel = NULL;
	}
	if (context->input_complex_spectrum != NULL) {
		fftwf_free(context->input_complex_spectrum);
		context->input_complex_spectrum = NULL;
	}
	if (context->output_complex_spectrum != NULL) {
		fftwf_free(context->output_complex_spectrum);
		context->output_complex_spectrum = NULL;
	}
	if (context != NULL) {
		free(context);
		context = NULL;
	}
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(Convolution, ARCH_SUFFIX)::CreateConvolve1DContext(
		size_t num_channel, LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
		size_t kernel_width, bool use_fft,
		LIBSAKURA_SYMBOL(Convolve1DContext) **context) const {
	CreateConvolve1DContextDefault(num_channel, kernel_type, kernel_width,
			use_fft, context);
}

void ADDSUFFIX(Convolution, ARCH_SUFFIX)::Convolve1D(
LIBSAKURA_SYMBOL(Convolve1DContext) **context,
		float input_spectrum[/*num_channels*/],
		bool const input_flag[/*num_channels*/],
		float output_spectrum[/*num_channels*/]) const {

	Convolve1DDefault(context, input_spectrum, input_flag, output_spectrum);

}

void ADDSUFFIX(Convolution, ARCH_SUFFIX)::DestroyConvolve1DContext(
LIBSAKURA_SYMBOL(Convolve1DContext) *context) const {
	DestroyConvolve1DContextDefault(context);

}

}
