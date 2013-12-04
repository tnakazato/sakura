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
	bool use_fft;
	size_t num_data;
	size_t kernel_width;
	fftwf_plan plan_real_to_complex_float;
	fftwf_plan plan_complex_to_real_float;
	fftwf_complex *fft_applied_complex_kernel;
	float *real_array;
};
}

namespace {

inline fftwf_complex* AllocateFFTArray(size_t num_data) {
	return (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * num_data);
}

inline void FreeFFTArray(fftwf_complex *ptr) {
	if (ptr != nullptr) {
		fftwf_free(ptr);
	}
}

inline void DestroyFFTPlan(fftwf_plan ptr) {
	if (ptr != nullptr) {
		fftwf_destroy_plan(ptr);
	}
}

inline void ApplyMaskByZero(size_t num_data, float const *input_data,
bool const *mask, float *output_data) {
	for (size_t i = 0; i < num_data; ++i) {
		if (mask[i]) {
			output_data[i] = input_data[i];
		} else {
			output_data[i] = 0.0;
		}
	}
}

inline bool FindFalseToGetStart(size_t num_data, float const input_data,
bool const mask, size_t* X1, float* Y1, size_t count) {
	if (!mask) {
		X1[count] = num_data; //- 1;
		Y1[count] = input_data;
		return true;
	}
	return false;
}
inline bool FindTrueToGetEnd(size_t num_data, float const input_data,
bool const mask, size_t* X2, float* Y2, size_t count) {
	if (mask) {
		X2[count] = num_data;
		Y2[count] = input_data;
		return true;
	}
	return false;
}

inline void ApplyMaskByLinearInterpolation(size_t num_data,
		float const *input_data,
		bool const *mask, float *output_data) {
	size_t X1[num_data], X2[num_data];
	float Y1[num_data], Y2[num_data];
	bool find_zero = true;
	if (!mask[0]) {
		find_zero = false;
	}
	size_t count = 0;
	for (size_t i = 1; i < num_data; ++i) { // finding (x1,y1),(x2,y2)
		if (find_zero) {
			if (FindFalseToGetStart(i - 1, input_data[i - 1], mask[i], X1, Y1,
					count)) {
				find_zero = false;
			}
		} else {
			if (FindTrueToGetEnd(i, input_data[i], mask[i], X2, Y2, count)) {
				find_zero = true;
				count++;
			} else if (i == (num_data - 1)) // achived final ch
				FindTrueToGetEnd(i, Y1[count], true, X2, Y2, count);
		}
	}
	size_t count_max = count;
	size_t n = 0;
	bool hit = false;
	for (size_t i = 0; i < num_data; ++i) { // calculate interpolation
		if (!hit) {
			output_data[i] = input_data[i];
			if (i == X1[n]) {
				hit = true;
			}
		} else if (hit && i < X2[n]) {
			output_data[i] = i * (Y2[n] - Y1[n]) / (X2[n] - X1[n])
					+ (X2[n] * Y1[n] - X1[n] * Y2[n]) / (X2[n] - X1[n]);
		} else if (hit && i == X2[n]) {
			hit = false;
			if (n < count_max)
				n++;
			output_data[i] = input_data[i];
		}
	}
}

inline void Create1DGaussianKernel(size_t num_data, size_t kernel_width,
		float* output_data) {
	float const reciprocal_of_denominator = 1.66510922231539551270632928979040
			/ kernel_width; // sqrt(log(16))/ kernel_width
	float const height = .939437278699651333772340328410 / kernel_width; // sqrt(8*log(2)/(2*M_PI)) / kernel_width
	size_t loop_max = num_data / 2;
	output_data[0] = height;
	float center = (num_data - 1) / 2.f;
	size_t middle = (num_data - 1) / 2;
	for (size_t j = 0; j < loop_max; ++j) {
		float value = (j - center) * reciprocal_of_denominator;
		output_data[middle + 1 + j] = output_data[middle - j] = height
				* exp(-(value * value));
	}
}

inline void Create1DKernel(size_t num_data,
LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type, size_t kernel_width,
		float* output_data) {
	switch (kernel_type) {
	case LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian): // Gaussian
		Create1DGaussianKernel(num_data, kernel_width, output_data);
		break;
	case LIBSAKURA_SYMBOL(Convolve1DKernelType_kBoxcar): // Boxcar
		break;
	case LIBSAKURA_SYMBOL(Convolve1DKernelType_kHanning): // Hanning
		break;
	case LIBSAKURA_SYMBOL(Convolve1DKernelType_kHamming): // Hamming
		break;
	default:
		assert(false);
		break;
	}
}

inline void CalculateConvolutionWithoutFFT(size_t num_data,
		float const *input_data, float const *input_kernel,
		float *output_data) {
	for (size_t j = 0; j < num_data; ++j) {
		float center = input_data[j] * input_kernel[0];
		float right = 0.0, left = 0.0;
		for (size_t i = 0; ((j + 1 + i < num_data) && (i < num_data / 2 - 1));
				++i) {
			left += input_data[j + 1 + i] * input_kernel[num_data - 1 - i];
		}
		for (size_t k = 0; k < j && (k < num_data / 2); ++k) {
			right += input_data[j - 1 - k] * input_kernel[k + 1];
		}
		output_data[j] = (left + center + right);
	}
}

inline void CreateConvolve1DContext(size_t num_data,
LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type, size_t kernel_width,
bool use_fft, LIBSAKURA_SYMBOL(Convolve1DContext)** context) {
	*context = nullptr;
	if (use_fft) {
		// real array for kernel
		std::unique_ptr<float[], LIBSAKURA_PREFIX::Memory> real_array_kernel(
				static_cast<float*>(LIBSAKURA_PREFIX::Memory::Allocate(
						sizeof(float) * (num_data + (2 * kernel_width)))),
				LIBSAKURA_PREFIX::Memory());
		if (real_array_kernel == nullptr) {
			throw std::bad_alloc();
		}
		// real array for input and output data
		std::unique_ptr<float[], LIBSAKURA_PREFIX::Memory> real_array(
				static_cast<float*>(LIBSAKURA_PREFIX::Memory::Allocate(
						sizeof(float) * (num_data + (2 * kernel_width)))),
				LIBSAKURA_PREFIX::Memory());
		if (real_array == nullptr) {
			throw std::bad_alloc();
		}
		// fft applied array for kernel
		std::unique_ptr<fftwf_complex[], decltype(&FreeFFTArray)> fft_applied_complex_kernel(
				AllocateFFTArray((num_data + (2 * kernel_width)) / 2 + 1),
				FreeFFTArray);
		if (fft_applied_complex_kernel == nullptr) {
			throw std::bad_alloc();
		}
		// fft applied array for input_data
		std::unique_ptr<fftwf_complex[], decltype(&FreeFFTArray)> fft_applied_complex_input_data(
				AllocateFFTArray((num_data + (2 * kernel_width)) / 2 + 1),
				FreeFFTArray);
		if (fft_applied_complex_input_data == nullptr) {
			throw std::bad_alloc();
		}
		// fft applied array for multiplied data
		std::unique_ptr<fftwf_complex[], decltype(&FreeFFTArray)> multiplied_complex_data(
				AllocateFFTArray((num_data + (2 * kernel_width)) / 2 + 1),
				FreeFFTArray);
		if (multiplied_complex_data == nullptr) {
			throw std::bad_alloc();
		}
		// create fft plan for kernel
		fftwf_plan plan_real_to_complex_float_kernel = fftwf_plan_dft_r2c_1d(
				num_data + (2 * kernel_width), real_array_kernel.get(),
				fft_applied_complex_kernel.get(),
				FFTW_ESTIMATE | FFTW_PRESERVE_INPUT); // added FFTW_PRESERVE_INPUT
		ScopeGuard guard_for_fft_plan_kernel([&]() {
			DestroyFFTPlan(plan_real_to_complex_float_kernel);
		});
		if (plan_real_to_complex_float_kernel == nullptr) {
			guard_for_fft_plan_kernel.Disable();
			throw std::bad_alloc();
		}
		// create fft plan for input data
		fftwf_plan plan_real_to_complex_float = fftwf_plan_dft_r2c_1d(
				num_data + (2 * kernel_width), real_array.get(),
				fft_applied_complex_input_data.get(),
				FFTW_ESTIMATE | FFTW_PRESERVE_INPUT); // added FFTW_PRESERVE_INPUT
		ScopeGuard guard_for_fft_plan([&]() {
			DestroyFFTPlan(plan_real_to_complex_float);
		});
		if (plan_real_to_complex_float == nullptr) {
			guard_for_fft_plan.Disable();
			throw std::bad_alloc();
		}
		// create ifft plan for output data
		fftwf_plan plan_complex_to_real_float = fftwf_plan_dft_c2r_1d(
				num_data + (2 * kernel_width), multiplied_complex_data.get(),
				real_array.get(),
				FFTW_ESTIMATE | FFTW_DESTROY_INPUT); // added FFTW_DESTROY_INPUT
		ScopeGuard guard_for_ifft_plan([&]() {
			DestroyFFTPlan(plan_complex_to_real_float);
		});
		if (plan_complex_to_real_float == nullptr) {
			guard_for_ifft_plan.Disable();
			throw std::bad_alloc();
		}
		Create1DKernel(num_data + (2 * kernel_width), kernel_type, kernel_width,
				real_array_kernel.get());
		fftwf_execute(plan_real_to_complex_float_kernel); // execute fft for kernel
		real_array_kernel.reset(nullptr);
		guard_for_fft_plan_kernel.CleanUpNow();
		// work_context
		std::unique_ptr<LIBSAKURA_SYMBOL(Convolve1DContext),
		LIBSAKURA_PREFIX::Memory> work_context(
				static_cast<LIBSAKURA_SYMBOL(Convolve1DContext)*>(LIBSAKURA_PREFIX::Memory::Allocate(
						sizeof(LIBSAKURA_SYMBOL(Convolve1DContext)))),
				LIBSAKURA_PREFIX::Memory());
		if (work_context == nullptr) {
			throw std::bad_alloc();
		}
		work_context->fft_applied_complex_kernel = // fft applied kernel
				fft_applied_complex_kernel.release();
		work_context->real_array = real_array.release(); // real_array
		work_context->plan_real_to_complex_float = plan_real_to_complex_float; // plan for real to complex
		guard_for_fft_plan.Disable();
		work_context->plan_complex_to_real_float = plan_complex_to_real_float; // plan for complex to real
		guard_for_ifft_plan.Disable();
		work_context->num_data = num_data; // number of data
		work_context->use_fft = use_fft; // use_fft flag
		work_context->kernel_width = kernel_width;
		*context = work_context.release(); // context
	} else { // without fft
		// real array for kernel
		std::unique_ptr<float[], LIBSAKURA_PREFIX::Memory> real_array_kernel(
				static_cast<float*>(LIBSAKURA_PREFIX::Memory::Allocate(
						sizeof(float) * num_data)), LIBSAKURA_PREFIX::Memory());
		if (real_array_kernel == nullptr) {
			throw std::bad_alloc();
		}
		Create1DKernel(num_data, kernel_type, kernel_width,
				real_array_kernel.get());
		// work_context
		std::unique_ptr<LIBSAKURA_SYMBOL(Convolve1DContext),
		LIBSAKURA_PREFIX::Memory> work_context(
				static_cast<LIBSAKURA_SYMBOL(Convolve1DContext)*>(LIBSAKURA_PREFIX::Memory::Allocate(
						sizeof(LIBSAKURA_SYMBOL(Convolve1DContext)))),
				LIBSAKURA_PREFIX::Memory());
		if (work_context == nullptr) {
			throw std::bad_alloc();
		}
		work_context->fft_applied_complex_kernel = nullptr;
		work_context->plan_complex_to_real_float = nullptr;
		work_context->plan_real_to_complex_float = nullptr;
		work_context->real_array = real_array_kernel.release(); // real_array for kernel
		work_context->num_data = num_data; // number of data
		work_context->use_fft = use_fft; // use_fft flag
		work_context->kernel_width = kernel_width;
		*context = work_context.release(); // context
	}
}

inline void Convolve1D(
LIBSAKURA_SYMBOL(Convolve1DContext) const *context, size_t num_data,
		float const input_data[/*num_data*/],
		bool const mask[/*num_data*/], float output_data[/*num_data*/]) {
	if (context->num_data != num_data) {
		throw std::invalid_argument(
				"num_data doesn't equal to context->num_data");
	}
	// real array for kernel
	std::unique_ptr<float[], LIBSAKURA_PREFIX::Memory> real_array(
			static_cast<float*>(LIBSAKURA_PREFIX::Memory::Allocate(
					sizeof(float) * (num_data + (2 * context->kernel_width)))),
			LIBSAKURA_PREFIX::Memory());
	if (real_array == nullptr) {
		throw std::bad_alloc();
	}
	if (context->use_fft) { // with fft
		// fft applied array for input_data (complex)
		std::unique_ptr<fftwf_complex[], decltype(&FreeFFTArray)> fft_applied_complex_input_data(
				AllocateFFTArray(
						(num_data + (2 * context->kernel_width)) / 2 + 1),
				FreeFFTArray);
		if (fft_applied_complex_input_data == nullptr) {
			throw std::bad_alloc();
		}
		// fft applied array for multiplied data (kernel * input_data) (complex)
		std::unique_ptr<fftwf_complex[], decltype(&FreeFFTArray)> multiplied_complex_data(
				AllocateFFTArray(
						(num_data + (2 * context->kernel_width)) / 2 + 1),
				FreeFFTArray);
		if (multiplied_complex_data == nullptr) {
			throw std::bad_alloc();
		}
		for (size_t i = 0; i < (num_data + (2 * context->kernel_width)); ++i) {
			real_array[i] = input_data[i];
		}
		for (size_t i = 0; i < 2 * context->kernel_width; ++i) {
			real_array[num_data + i] = 0.0; // zero padding to avoid polution
		}
		fftwf_execute_dft_r2c(context->plan_real_to_complex_float,
				real_array.get(), //const_cast<float*>(input_data),
				fft_applied_complex_input_data.get());
		float scale = 1.0 / (num_data + (2 * context->kernel_width));
		for (size_t i = 0; i < (num_data + (2 * context->kernel_width)) / 2 + 1;
				++i) {
			multiplied_complex_data[i][0] =
					(context->fft_applied_complex_kernel[i][0]
							* fft_applied_complex_input_data[i][0]
							- context->fft_applied_complex_kernel[i][1]
									* fft_applied_complex_input_data[i][1])
							* scale;
			multiplied_complex_data[i][1] =
					(context->fft_applied_complex_kernel[i][0]
							* fft_applied_complex_input_data[i][1]
							+ context->fft_applied_complex_kernel[i][1]
									* fft_applied_complex_input_data[i][0])
							* scale;
		}
		fftwf_execute_dft_c2r(context->plan_complex_to_real_float,
				multiplied_complex_data.get(), output_data);
	} else { // without fft
		CalculateConvolutionWithoutFFT(num_data, const_cast<float*>(input_data),
				context->real_array, output_data);
	}
}

inline void DestroyConvolve1DContext(
LIBSAKURA_SYMBOL(Convolve1DContext)* context) {
	if (context != nullptr) {
		if (context->plan_real_to_complex_float != nullptr) {
			DestroyFFTPlan(context->plan_real_to_complex_float);
		}
		if (context->plan_complex_to_real_float != nullptr) {
			DestroyFFTPlan(context->plan_complex_to_real_float);
		}
		if (context->fft_applied_complex_kernel != nullptr) {
			FreeFFTArray(context->fft_applied_complex_kernel);
		}
		if (context->real_array != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(context->real_array);
		}
		LIBSAKURA_PREFIX::Memory::Free(context);
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
LIBSAKURA_SYMBOL(Convolve1DContext) const *context, size_t num_data,
		float const input_data[/*num_data*/],
		bool const mask[/*num_data*/], float output_data[/*num_data*/]) const {
	::Convolve1D(context, num_data, input_data, mask, output_data);

}

void ADDSUFFIX(Convolution, ARCH_SUFFIX)::DestroyConvolve1DContext(
LIBSAKURA_SYMBOL(Convolve1DContext) *context) const {
	::DestroyConvolve1DContext(context);
}

}
