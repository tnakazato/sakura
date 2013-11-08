#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fftw3.h>
#include <memory>

#include <libsakura/sakura.h>
#include <libsakura/optimized_implementation_factory_impl.h>
#include <libsakura/localdef.h>
#include <libsakura/logger.h>
#include <libsakura/memory_manager.h>

extern "C" {
struct LIBSAKURA_SYMBOL(Convolve1DContext) {
	size_t num_data;
	fftwf_plan plan_real_to_complex_float;
	fftwf_plan plan_complex_to_real_float;
	fftwf_complex *fft_applied_complex_input_data;
	fftwf_complex *multiplied_complex_data;
	fftwf_complex *fft_applied_complex_kernel;
	float real_array[0];
};
}

namespace {

inline void CreateConvolve1DContextDefault(size_t num_data,
LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type, size_t kernel_width,
bool use_fft, LIBSAKURA_SYMBOL(Convolve1DContext)** context) {

	float const sqrt_ln16 = 1.66510922231539551270632928979040;
	float const height = .939437278699651333772340328410 / kernel_width; // sqrt((8*ln2)/(2*pi)) / kernel_width
	std::unique_ptr<float[]> work_real_array(new float[num_data]);

	size_t loop_max = num_data / 2;
	work_real_array[0] = height;

	float center = (num_data - 1) / 2.f;
	size_t middle = (num_data - 1) / 2;
	switch (kernel_type) {
	case LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian):
		for (size_t j = 0; j < loop_max; ++j) {
			float value = (j - center) * sqrt_ln16 / kernel_width;
			work_real_array[middle + 1 + j] = work_real_array[middle - j] =
					height * exp(-(value * value));
		}
		break;
	case LIBSAKURA_SYMBOL(Convolve1DKernelType_kBoxcar):
		// BoxCar
		break;
	case LIBSAKURA_SYMBOL(Convolve1DKernelType_kHanning):
		// Hanning
		break;
	case LIBSAKURA_SYMBOL(Convolve1DKernelType_kHamming):
		// Hamming
		break;
	default:
		break;
	}

	if (use_fft) {
		//kernel array
		fftwf_complex *work_fft_applied_complex_kernel =
				(fftwf_complex*) fftwf_malloc(
						sizeof(fftwf_complex) * (num_data / 2 + 1));
		if (work_fft_applied_complex_kernel == nullptr) {
			throw "";
		}

		// create plan for kernel
		fftwf_plan work_plan_real_to_complex_float = fftwf_plan_dft_r2c_1d(
				num_data, work_real_array.get(),
				work_fft_applied_complex_kernel,
				FFTW_ESTIMATE);

		if (work_plan_real_to_complex_float == nullptr) {
			throw "";
		} else {
			// alloc real_array[j]
			*context = (LIBSAKURA_SYMBOL(Convolve1DContext) *) malloc(
					sizeof(LIBSAKURA_SYMBOL(Convolve1DContext))
							+ sizeof(float) * num_data);
			if (*context == nullptr) {
				throw "";
			}
			(*context)->plan_real_to_complex_float =
					work_plan_real_to_complex_float;
			fftwf_execute((*context)->plan_real_to_complex_float); // (4) fftwf_execute kernel
		}

		(*context)->fft_applied_complex_kernel =
				work_fft_applied_complex_kernel;

		//fft_applied_complex_input_data
		fftwf_complex *work_fft_applied_complex_input_data =
				(fftwf_complex*) fftwf_malloc(
						sizeof(fftwf_complex) * (num_data / 2 + 1));
		if (work_fft_applied_complex_input_data == nullptr) {
			throw "";
		} else {
			(*context)->fft_applied_complex_input_data =
					work_fft_applied_complex_input_data;
		}

		// create plan
		work_plan_real_to_complex_float = fftwf_plan_dft_r2c_1d(num_data,
				(*context)->real_array,
				(*context)->fft_applied_complex_input_data, FFTW_ESTIMATE);
		if (work_plan_real_to_complex_float == nullptr) {
			throw "";
		} else {
			(*context)->plan_real_to_complex_float =
					work_plan_real_to_complex_float;
		}

		//multiplied_complex_data
		fftwf_complex *work_multiplied_complex_data =
				(fftwf_complex*) fftwf_malloc(
						sizeof(fftwf_complex) * (num_data / 2 + 1));

		if (work_multiplied_complex_data == nullptr) {
			throw "";
		} else {
			(*context)->multiplied_complex_data = work_multiplied_complex_data;
		}

		// create plan
		fftwf_plan work_plan_complex_to_real_float = fftwf_plan_dft_c2r_1d(
				num_data, (*context)->multiplied_complex_data,
				(*context)->real_array, FFTW_ESTIMATE);

		if (work_plan_complex_to_real_float == nullptr) {
			throw "";
		} else {
			(*context)->plan_complex_to_real_float =
					work_plan_complex_to_real_float;
		}

		// copy from work to (*context)->real_array[j]
		for (size_t j = 0; j < num_data; ++j) {
			(*context)->real_array[j] = work_real_array[j];
		}

		(*context)->num_data = num_data;

	} else {
		// not implemented yet
	}

}

inline void Convolve1DDefault(LIBSAKURA_SYMBOL(Convolve1DContext) *context,
		size_t num_data, float input_data[/*num_data*/],
		bool const input_flag[/*num_data*/], float output_data[/*num_data*/]) {

	for (size_t i = 0; i < num_data; ++i) {
		(context)->real_array[i] = input_data[i];
	}

	if ((context)->plan_real_to_complex_float != nullptr)
		fftwf_execute((context)->plan_real_to_complex_float);

	float scale = 1.0 / num_data;
	for (size_t i = 0; i < num_data / 2 + 1; ++i) {
		(context)->multiplied_complex_data[i][0] =
				((context)->fft_applied_complex_kernel[i][0]
						* (context)->fft_applied_complex_input_data[i][0]
						- (context)->fft_applied_complex_kernel[i][1]
								* (context)->fft_applied_complex_input_data[i][1])
						* scale;
		(context)->multiplied_complex_data[i][1] =
				((context)->fft_applied_complex_kernel[i][0]
						* (context)->fft_applied_complex_input_data[i][1]
						+ (context)->fft_applied_complex_kernel[i][1]
								* (context)->fft_applied_complex_input_data[i][0])
						* scale;
	}

	if ((context)->plan_complex_to_real_float != nullptr)
		fftwf_execute((context)->plan_complex_to_real_float);

	for (size_t i = 0; i < num_data; ++i) {
		output_data[i] = (context)->real_array[i];
	}
}

inline void DestroyConvolve1DContextDefault(
LIBSAKURA_SYMBOL(Convolve1DContext)* context) {
	std::cout << " DestroyConvolve1DContextEigen function is called"
			<< std::endl;
	if (context->plan_real_to_complex_float != nullptr) {
		fftwf_destroy_plan(context->plan_real_to_complex_float);
		context->plan_real_to_complex_float = nullptr;
	}
	if (context->plan_complex_to_real_float != nullptr) {
		fftwf_destroy_plan(context->plan_complex_to_real_float);
		context->plan_complex_to_real_float = nullptr;
	}
	if (context->fft_applied_complex_kernel != nullptr) {
		fftwf_free(context->fft_applied_complex_kernel);
		context->fft_applied_complex_kernel = nullptr;
	}
	if (context->fft_applied_complex_input_data != nullptr) {
		fftwf_free(context->fft_applied_complex_input_data);
		context->fft_applied_complex_input_data = nullptr;
	}
	if (context->multiplied_complex_data != nullptr) {
		fftwf_free(context->multiplied_complex_data);
		context->multiplied_complex_data = nullptr;
	}
	if (context != nullptr) {
		free(context);
		context = nullptr;
	}
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(Convolution, ARCH_SUFFIX)::CreateConvolve1DContext(
		size_t num_data, LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
		size_t kernel_width, bool use_fft,
		LIBSAKURA_SYMBOL(Convolve1DContext) **context) const {
	CreateConvolve1DContextDefault(num_data, kernel_type, kernel_width, use_fft,
			context);
}

void ADDSUFFIX(Convolution, ARCH_SUFFIX)::Convolve1D(
LIBSAKURA_SYMBOL(Convolve1DContext) *context, size_t num_data,
		float input_data[/*num_data*/],
		bool const input_flag[/*num_data*/],
		float output_data[/*num_data*/]) const {
	Convolve1DDefault(context, num_data, input_data, input_flag, output_data);

}

void ADDSUFFIX(Convolution, ARCH_SUFFIX)::DestroyConvolve1DContext(
LIBSAKURA_SYMBOL(Convolve1DContext) *context) const {
	DestroyConvolve1DContextDefault(context);
}

}
