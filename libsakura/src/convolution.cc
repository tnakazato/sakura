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
	float *real_array;
};
}

namespace {

inline void FreeFFTArray(fftwf_complex *ptr) {
	if (ptr != nullptr) {
		fftwf_free(ptr);
		ptr = nullptr;
	}
}

inline fftwf_complex* AllocateFFTArray(size_t num_data) {
	return (fftwf_complex*) fftwf_malloc(
			sizeof(fftwf_complex) * (num_data / 2 + 1));
}

inline LIBSAKURA_SYMBOL(Status) CreateConvolve1DContext(size_t num_data,
LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type, size_t kernel_width,
bool use_fft, LIBSAKURA_SYMBOL(Convolve1DContext)** context) {
	float const reciprocal_of_denominator = 1.66510922231539551270632928979040
			/ kernel_width; // sqrt(log(16))/ kernel_width
	float const height = .939437278699651333772340328410 / kernel_width; // sqrt(8*log(2)/(2*M_PI)) / kernel_width

	// real array for kernel
	std::unique_ptr<float[], LIBSAKURA_PREFIX::Memory> real_array_kernel(
			static_cast<float*>(LIBSAKURA_PREFIX::Memory::Allocate(
					sizeof(float) * num_data)), LIBSAKURA_PREFIX::Memory());
	if (real_array_kernel == nullptr) {
		*context = nullptr;
		return LIBSAKURA_SYMBOL(Status_kNG);
	}
	size_t loop_max = num_data / 2;
	real_array_kernel[0] = height;
	float center = (num_data - 1) / 2.f;
	size_t middle = (num_data - 1) / 2;
	switch (kernel_type) {
	case LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian):
		for (size_t j = 0; j < loop_max; ++j) {
			float value = (j - center) * reciprocal_of_denominator;
			real_array_kernel[middle + 1 + j] = real_array_kernel[middle - j] =
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
		// fft applied array for kernel
		std::unique_ptr<fftwf_complex[], decltype(&FreeFFTArray)> fft_applied_complex_kernel(
				AllocateFFTArray(num_data), FreeFFTArray);
		if (fft_applied_complex_kernel == nullptr) {
			*context = nullptr;
			return LIBSAKURA_SYMBOL(Status_kNG);
		}
		// real array for input and output data
		std::unique_ptr<float[], LIBSAKURA_PREFIX::Memory> real_array(
				static_cast<float*>(LIBSAKURA_PREFIX::Memory::Allocate(
						sizeof(float) * num_data)), LIBSAKURA_PREFIX::Memory());
		if (real_array == nullptr) {
			*context = nullptr;
			return LIBSAKURA_SYMBOL(Status_kNG);
		}
		// fft applied array foe input_data
		std::unique_ptr<fftwf_complex[], decltype(&FreeFFTArray)> fft_applied_complex_input_data(
				AllocateFFTArray(num_data), FreeFFTArray);
		if (fft_applied_complex_input_data == nullptr) {
			*context = nullptr;
			return LIBSAKURA_SYMBOL(Status_kNG);
		}
		// fft applied array for multiplied data
		std::unique_ptr<fftwf_complex[], decltype(&FreeFFTArray)> multiplied_complex_data(
				AllocateFFTArray(num_data), FreeFFTArray);
		if (multiplied_complex_data == nullptr) {
			*context = nullptr;
			return LIBSAKURA_SYMBOL(Status_kNG);
		}
		// plan of fft for kernel
		fftwf_plan plan_real_to_complex_float_kernel = fftwf_plan_dft_r2c_1d(
				num_data, real_array_kernel.get(),
				fft_applied_complex_kernel.get(),
				FFTW_ESTIMATE);
		if (plan_real_to_complex_float_kernel == nullptr) {
			*context = nullptr;
			return LIBSAKURA_SYMBOL(Status_kNG);
		}
		// create context
		*context =
				static_cast<LIBSAKURA_SYMBOL(Convolve1DContext) *>(LIBSAKURA_PREFIX::Memory::Allocate(
						sizeof(LIBSAKURA_SYMBOL(Convolve1DContext))));
		if (*context == nullptr) {
			fftwf_destroy_plan(plan_real_to_complex_float_kernel);
			plan_real_to_complex_float_kernel = nullptr;
			*context = nullptr;
			return LIBSAKURA_SYMBOL(Status_kNG);
		}

		// plan of fft for input data
		(*context)->plan_real_to_complex_float = fftwf_plan_dft_r2c_1d(num_data,
				real_array.get(), fft_applied_complex_input_data.get(),
				FFTW_ESTIMATE);
		if ((*context)->plan_real_to_complex_float == nullptr) {
			fftwf_destroy_plan(plan_real_to_complex_float_kernel);
			plan_real_to_complex_float_kernel = nullptr;
			LIBSAKURA_PREFIX::Memory::Free(context);
			*context = nullptr;
			return LIBSAKURA_SYMBOL(Status_kNG);
		}

		// plan of ifft for output data
		(*context)->plan_complex_to_real_float = fftwf_plan_dft_c2r_1d(num_data,
				multiplied_complex_data.get(), real_array.get(), FFTW_ESTIMATE);
		if ((*context)->plan_complex_to_real_float == nullptr) {
			fftwf_destroy_plan(plan_real_to_complex_float_kernel);
			plan_real_to_complex_float_kernel = nullptr;
			fftwf_destroy_plan((*context)->plan_real_to_complex_float);
			(*context)->plan_real_to_complex_float = nullptr;
			LIBSAKURA_PREFIX::Memory::Free(context);
			*context = nullptr;
			return LIBSAKURA_SYMBOL(Status_kNG);
		}
		(*context)->real_array = real_array.release(); // real_array
		(*context)->fft_applied_complex_input_data = // fft_applied_complex_input_data
				fft_applied_complex_input_data.release();
		(*context)->multiplied_complex_data = multiplied_complex_data.release(); // multiplied_complex_data
		fftwf_execute(plan_real_to_complex_float_kernel); // execute fft for kernel
		(*context)->fft_applied_complex_kernel = // fft applied kernel
				fft_applied_complex_kernel.release();
		(*context)->num_data = num_data; // number of data
	} else {
		// not implemented yet
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

inline LIBSAKURA_SYMBOL(Status) Convolve1D(
		LIBSAKURA_SYMBOL(Convolve1DContext) *context, size_t num_data,
		float input_data[/*num_data*/],
		bool const input_flag[/*num_data*/], float output_data[/*num_data*/]) {
	if ((context)->num_data != num_data) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	} else if (context != nullptr) {
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
	return LIBSAKURA_SYMBOL(Status_kOK);
}

inline void DestroyConvolve1DContext(
LIBSAKURA_SYMBOL(Convolve1DContext)* context) {
	if (context != nullptr) {
		if (context->plan_real_to_complex_float != nullptr) {
			fftwf_destroy_plan(context->plan_real_to_complex_float);
			context->plan_real_to_complex_float = nullptr;
		}
		if (context->plan_complex_to_real_float != nullptr) {
			fftwf_destroy_plan(context->plan_complex_to_real_float);
			context->plan_complex_to_real_float = nullptr;
		}
		if (context->fft_applied_complex_kernel != nullptr) {
			FreeFFTArray(context->fft_applied_complex_kernel);
		}
		if (context->fft_applied_complex_input_data != nullptr) {
			FreeFFTArray(context->fft_applied_complex_input_data);
		}
		if (context->multiplied_complex_data != nullptr) {
			FreeFFTArray(context->multiplied_complex_data);
		}
		LIBSAKURA_PREFIX::Memory::Free(context);
		context = nullptr;
	}
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(Convolution, ARCH_SUFFIX)::CreateConvolve1DContext(
		size_t num_data, LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
		size_t kernel_width, bool use_fft,
		LIBSAKURA_SYMBOL(Convolve1DContext) **context) const {
	::CreateConvolve1DContext(num_data, kernel_type, kernel_width, use_fft,
			context);
}

void ADDSUFFIX(Convolution, ARCH_SUFFIX)::Convolve1D(
LIBSAKURA_SYMBOL(Convolve1DContext) *context, size_t num_data,
		float input_data[/*num_data*/],
		bool const input_flag[/*num_data*/],
		float output_data[/*num_data*/]) const {
	::Convolve1D(context, num_data, input_data, input_flag, output_data);

}

void ADDSUFFIX(Convolution, ARCH_SUFFIX)::DestroyConvolve1DContext(
LIBSAKURA_SYMBOL(Convolve1DContext) *context) const {
	::DestroyConvolve1DContext(context);
}

}
