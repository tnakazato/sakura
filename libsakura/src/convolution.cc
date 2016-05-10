/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2016
 * National Astronomical Observatory of Japan
 * 2-21-1, Osawa, Mitaka, Tokyo, 181-8588, Japan.
 * 
 * This file is part of Sakura.
 * 
 * Sakura is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the 
 * Free Software Foundation, either version 3 of the License, or (at your 
 * option) any later version.
 * 
 * Sakura is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License 
 * along with Sakura.  If not, see <http://www.gnu.org/licenses/>.
 * @SAKURA_LICENSE_HEADER_END@
 */
/**
 * convolution.cc
 *
 *  Created on: 2013/11/25
 *      Author: shinnosuke kawakami
 */

#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fftw3.h>
#include <memory>
#include <climits>
#include <fstream>
#include <sstream>

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
#include <libsakura/logger.h>
#include <libsakura/memory_manager.h>

extern "C" {
struct LIBSAKURA_SYMBOL(Convolve1DContextFloat) {
	size_t kernel_width;
	fftwf_plan plan_real_to_complex_float;
	fftwf_plan plan_complex_to_real_float;
	float *real_array;
	float *imag_array;
	void *real_array_work;
	void *imag_array_work;
	float *real_kernel_array;
	float *imag_kernel_array;
	void *real_kernel_array_work;
	void *imag_kernel_array_work;
};
}

namespace {

auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("Convolution");

inline fftwf_complex* AllocateFFTArray(size_t num_data) {
	return reinterpret_cast<fftwf_complex*>(fftwf_malloc(
			sizeof(fftwf_complex) * num_data));
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

/**
 * @brief Create 1 dimensional Gaussian kernel
 *
 * @param num_kernel number of elements in kernel
 * @param use_fft
 * @param kernel_width FWHM of the kernel
 * @param kernel output kernel
 */
inline void Create1DGaussianKernelFloat(size_t num_kernel, float kernel_width,
		float* kernel) {
	assert(num_kernel > 0);
	assert(kernel_width > 0.0f);

	// special treatment of the case that num_kernel is 1
	if (num_kernel == 1) {
		kernel[0] = 1.0;
		return;
	}

	assert((2 * num_kernel - 1) >= 0);
	double const reciprocal_of_denominator = 1.66510922231539551270632928979040
			/ static_cast<double>(kernel_width); // sqrt(log(16)) / kernel_width
	double const peak_value = .939437278699651333772340328410
			/ static_cast<double>(kernel_width); // sqrt(8*log(2)/(2*M_PI)) / kernel_width
	size_t peak_location = (num_kernel) / 2;

	auto kernel_aligned = AssumeAligned(kernel);

	double kernel_sum = 0.0;
	for (size_t i = 0; i < num_kernel; ++i) {
		double value = (static_cast<double>(i)
				- static_cast<double>(peak_location))
				* reciprocal_of_denominator;
		kernel_aligned[i] = peak_value * exp(-(value * value));
		kernel_sum += kernel_aligned[i];
	}

	// Normalization
	for (size_t i = 0; i < num_kernel; ++i) {
		kernel_aligned[i] /= kernel_sum;
	}
}

/**
 *
 * @param num_data
 * @param src
 * @param dst
 * @tparam T
 */
template<typename T>
inline void RevertData1D(size_t num_data, T const src[/*num_data*/],
		T dst[/*num_data*/]) {
	for (size_t i = 0; i < num_data; ++i) {
		dst[i] = src[num_data - 1 - i];
	}
}

/**
 *
 * @param num_data
 * @param input_data_arg
 * @param num_kernel
 * @param kernel_arg
 * @param output_data_arg
 */
inline void ConvolutionWithoutFFT(size_t num_kernel, float const *kernel_arg,
		size_t num_data, float const *input_data_arg,
		bool const *input_mask_arg, float *output_data_arg, float *output_weight_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(kernel_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(input_data_arg));
	uint8_t const *mask8 = reinterpret_cast<uint8_t const *>(input_mask_arg);
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask8));
	assert(LIBSAKURA_SYMBOL(IsAligned)(output_data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(output_weight_arg));
	// Need to revert kernel array in direct convolution.
	SIMD_ALIGN
	float revert_kernel[num_kernel];
	assert(LIBSAKURA_SYMBOL(IsAligned)(revert_kernel));
	RevertData1D(num_kernel, kernel_arg, revert_kernel);

	auto input_data = AssumeAligned(input_data_arg);
	auto input_mask = AssumeAligned(mask8);
	auto kernel = AssumeAligned(revert_kernel);
	auto output_data = AssumeAligned(output_data_arg);
	auto weight_data = AssumeAligned(output_weight_arg);
	STATIC_ASSERT(sizeof(input_mask_arg[0]) == sizeof(input_mask[0]));
	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);
	// the center index in *reverted* kernel
	size_t center_index = num_kernel - 1 - (num_kernel / 2);
	// the number of elements in the left half of kernel (excluding center)
	size_t left_half = center_index;
	// the number of elements in the right half of kernel (excluding center)
	size_t right_half = num_kernel - left_half - 1;
	for (size_t i = 0; i < num_data; ++i) {
		size_t lower_index = (i >= center_index) ? i - left_half : 0;
		size_t kernel_start = (i >= center_index) ? 0 : center_index - i;
		size_t kernel_end =
				(num_data - 1 - i >= right_half) ?
						num_kernel :
						num_kernel - (right_half - (num_data - 1 - i)); // end index +1
		size_t num_convolve = kernel_end - kernel_start;
		assert(num_convolve <= num_kernel);
		double value = 0.0;
		double weight = 0.0;
		for (size_t j = 0; j < num_convolve; ++j) {
			double local_weight = static_cast<double>(kernel[kernel_start + j])
					* static_cast<double>(input_mask[lower_index + j]);
			value += (local_weight
					* static_cast<double>(input_data[(lower_index + j)]));
			weight += local_weight;
		}
		output_data[i] = (weight == 0.0 ? 0.0 : value / weight);
	}
}

/**
 *
 * @param num_data
 * @param src
 * @param dst
 * @tparam T
 */
template<typename T>
inline void FlipData1D(size_t num_data, T const src[], T dst[]) {
	size_t offset2 = (num_data + 1) / 2;
	size_t offset1 = num_data / 2;
	for (size_t i = 0; i < offset2; ++i) {
		dst[i] = src[i + offset1];
	}
	for (size_t i = offset2; i < num_data; ++i) {
		dst[i] = src[i - offset2];
	}
}

inline void CreateConvolve1DContextFFTFloat(size_t num_kernel,
		float const kernel[],
		LIBSAKURA_SYMBOL(Convolve1DContextFloat)** context) {

	assert(context != nullptr);
	assert(kernel != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(kernel));
	std::unique_ptr<LIBSAKURA_SYMBOL(Convolve1DContextFloat),
	LIBSAKURA_PREFIX::Memory> work_context(
			static_cast<LIBSAKURA_SYMBOL(Convolve1DContextFloat)*>(LIBSAKURA_PREFIX::Memory::Allocate(
					sizeof(LIBSAKURA_SYMBOL(Convolve1DContextFloat)))),
			LIBSAKURA_PREFIX::Memory());
	if (work_context == nullptr) {
		throw std::bad_alloc();
	}
	float *real_kernel_array = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> real_kernel_array_work(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(float) * num_kernel, &real_kernel_array));
	assert(LIBSAKURA_SYMBOL(IsAligned)(real_kernel_array));
	float *imag_kernel_array = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> imag_kernel_array_work(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(float) * num_kernel, &imag_kernel_array));
	assert(LIBSAKURA_SYMBOL(IsAligned)(imag_kernel_array));
	FlipData1D(num_kernel, kernel, real_kernel_array);
	float *real_array = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> real_array_work(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(float) * num_kernel, &real_array));
	assert(LIBSAKURA_SYMBOL(IsAligned)(real_array));
	float *imag_array = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> imag_array_work(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(float) * num_kernel, &imag_array));
	assert(LIBSAKURA_SYMBOL(IsAligned)(imag_array));

	// Create plan for forward/backward FFT
	// make use of guru interface
	int rank = 1;
	fftw_iodim dims[1];
	dims[0].n = num_kernel;
	dims[0].is = 1;
	dims[0].os = 1;
	fftw_iodim howmany_dims[1];
	howmany_dims[0].n = 1;
	howmany_dims[0].is = 1;
	howmany_dims[0].os = 1;
	fftwf_plan plan_real_to_complex_float = fftwf_plan_guru_split_dft_r2c(rank,
			dims, rank, howmany_dims, real_array, real_kernel_array, imag_array,
			FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
	ScopeGuard guard_for_fft_plan([&]() {
		DestroyFFTPlan(plan_real_to_complex_float);
	});
	if (plan_real_to_complex_float == nullptr) {
		guard_for_fft_plan.Disable();
		throw std::bad_alloc();
	}
	fftwf_plan plan_complex_to_real_float = fftwf_plan_guru_split_dft_c2r(rank,
			dims, rank, howmany_dims, real_array, imag_array, imag_kernel_array,
			FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
	ScopeGuard guard_for_ifft_plan([&]() {
		DestroyFFTPlan(plan_complex_to_real_float);
	});
	if (plan_complex_to_real_float == nullptr) {
		guard_for_ifft_plan.Disable();
		throw std::bad_alloc();
	}

	// Obtain Fourier transform of kernel
	fftwf_plan plan_real_to_complex_float_kernel =
			fftwf_plan_guru_split_dft_r2c(rank, dims, rank, howmany_dims,
					real_kernel_array, real_kernel_array, imag_kernel_array,
					FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
	ScopeGuard guard_for_fft_plan_kernel([&]() {
		DestroyFFTPlan(plan_real_to_complex_float_kernel);
	});
	if (plan_real_to_complex_float_kernel == nullptr) {
		guard_for_fft_plan_kernel.Disable();
		throw std::bad_alloc();
	}
	fftwf_execute(plan_real_to_complex_float_kernel);
	guard_for_fft_plan_kernel.CleanUpNow();

	work_context->kernel_width = num_kernel; //kernel_width;
	work_context->plan_real_to_complex_float = plan_real_to_complex_float;
	work_context->plan_complex_to_real_float = plan_complex_to_real_float;
	work_context->real_array = real_array;
	work_context->imag_array = imag_array;
	work_context->real_array_work = real_array_work.release();
	work_context->imag_array_work = imag_array_work.release();
	work_context->real_kernel_array = real_kernel_array;
	work_context->imag_kernel_array = imag_kernel_array;
	work_context->real_kernel_array_work = real_kernel_array_work.release();
	work_context->imag_kernel_array_work = imag_kernel_array_work.release();
	guard_for_fft_plan.Disable();
	guard_for_ifft_plan.Disable();
	*context = work_context.release();
}

inline void ConvolutionWithFFT(
LIBSAKURA_SYMBOL(Convolve1DContextFloat) const *context, size_t num_data,
		float const input_data_arg[/*num_data*/],
		float output_data_arg[/*num_data*/]) {
	assert(context->num_data == num_data);
	assert(LIBSAKURA_SYMBOL(IsAligned)(input_data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(input_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(output_data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(output_mask_arg));
	auto input_data = AssumeAligned(input_data_arg);
	auto output_data = AssumeAligned(output_data_arg);
	size_t num_fft_data = num_data / 2 + 1;
	auto real_input_array = context->real_array;
	auto imag_input_array = context->imag_array;
	auto real_kernel_array = context->real_kernel_array;
	auto imag_kernel_array = context->imag_kernel_array;

	// FFT input array
	fftwf_execute_split_dft_r2c(context->plan_real_to_complex_float,
			const_cast<float *>(input_data), real_input_array,
			imag_input_array);

	// Complex multiplication in Fourier domain
	float scale = 1.0f / static_cast<float>(num_data);
	for (size_t i = 0; i < num_fft_data; ++i) {
		auto real_data = real_input_array[i];
		auto imag_data = imag_input_array[i];
		auto real_kernel = real_kernel_array[i];
		auto imag_kernel = imag_kernel_array[i];
		real_input_array[i] =
				(real_kernel * real_data - imag_kernel * imag_data) * scale;
		imag_input_array[i] =
				(real_kernel * imag_data + imag_kernel * real_data) * scale;
	}

	// Inverse FFT
	fftwf_execute_split_dft_c2r(context->plan_complex_to_real_float,
			real_input_array, imag_input_array, output_data);
}

inline void DestroyConvolve1DContextFloat(
LIBSAKURA_SYMBOL(Convolve1DContextFloat)* context) {
	if (context != nullptr) {
		DestroyFFTPlan(context->plan_real_to_complex_float);
		DestroyFFTPlan(context->plan_complex_to_real_float);
		if (context->real_kernel_array_work != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(context->real_kernel_array_work);
		}
		if (context->real_array_work != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(context->real_array_work);
		}
		if (context->imag_kernel_array_work != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(context->imag_kernel_array_work);
		}
		if (context->imag_array_work != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(context->imag_array_work);
		}
		LIBSAKURA_PREFIX::Memory::Free(context);
	}
}

} /* anonymous namespace */

#define CHECK_ARGS_WITH_MESSAGE(x,msg) do { \
	if (!(x)) { \
		LOG4CXX_ERROR(logger, msg); \
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument); \
	} \
} while (false)

#define CHECK_ARGS(x) do { \
	if (!(x)) { \
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument); \
	} \
} while (false)

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateGaussianKernelFloat)(
		float kernel_width, size_t num_kernel, float kernel[]) noexcept {
	CHECK_ARGS(0 < num_kernel && num_kernel <= INT_MAX);
	CHECK_ARGS(kernel_width > 0.0f);
	CHECK_ARGS(std::isfinite(kernel_width));
	CHECK_ARGS(kernel != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned(kernel)));
	try {
		Create1DGaussianKernelFloat(num_kernel, kernel_width, kernel);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const std::runtime_error &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kNG);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateConvolve1DContextFFTFloat)(
		size_t num_kernel, float const kernel[/*num_kernel*/],
		LIBSAKURA_SYMBOL(Convolve1DContextFloat) **context) noexcept {
	CHECK_ARGS_WITH_MESSAGE(0 < num_kernel && num_kernel <= INT_MAX,
			"num_kernel must satisfy '0 < num_kernel <= INT_MAX'");
	CHECK_ARGS_WITH_MESSAGE(kernel != nullptr, "kernel should not be NULL");
	CHECK_ARGS_WITH_MESSAGE(LIBSAKURA_SYMBOL(IsAligned)(kernel),
			"kernel should be aligned");
	CHECK_ARGS_WITH_MESSAGE(context != nullptr, "context should not be NULL");
	try {
		CreateConvolve1DContextFFTFloat(num_kernel, kernel, context);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const std::runtime_error &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kNG);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Convolve1DFloat)(
		size_t num_kernel, float const kernel[/*num_kernel*/], size_t num_data,
		float const input_data[/*num_data*/],
		bool const input_mask[/*num_data*/], float output_data[/*num_data*/],
		float output_weight[/*num_data*/]) noexcept {
	CHECK_ARGS_WITH_MESSAGE(0 < num_kernel && num_kernel <= INT_MAX,
			"num_data must be 0 < num_kernel <= INT_MAX");
	CHECK_ARGS_WITH_MESSAGE(0 < num_data && num_data <= INT_MAX,
			"num_data must be 0 < num_data <= INT_MAX");
	CHECK_ARGS_WITH_MESSAGE(
			kernel != nullptr && LIBSAKURA_SYMBOL(IsAligned)(kernel),
			"invalid kernel");
	CHECK_ARGS_WITH_MESSAGE(
			input_data != nullptr && LIBSAKURA_SYMBOL(IsAligned)(input_data),
			"invalid input_data");
	CHECK_ARGS_WITH_MESSAGE(
			input_mask != nullptr && LIBSAKURA_SYMBOL(IsAligned)(input_mask),
			"invalid input_mask");
	CHECK_ARGS_WITH_MESSAGE(
			output_data != nullptr && LIBSAKURA_SYMBOL(IsAligned)(output_data),
			"invalid output_data");
	CHECK_ARGS_WITH_MESSAGE(
			output_weight != nullptr && LIBSAKURA_SYMBOL(IsAligned)(output_weight),
			"invalid output_weight");
	CHECK_ARGS_WITH_MESSAGE(input_data != output_data,
			"in-place operation is not supported in convolution without FFT. ");
	//assert(fftw_alignment_of((double *)input_data) == 0);
	//assert(fftw_alignment_of((double *)output_data) == 0);
	try {
		ConvolutionWithoutFFT(num_kernel, kernel, num_data, input_data,
				input_mask, output_data, output_weight);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Convolve1DFFTFloat)(
LIBSAKURA_SYMBOL(Convolve1DContextFloat) const *context, size_t num_data,
		float const input_data[/*num_data*/], float output_data[/*num_data*/])
				noexcept {
	CHECK_ARGS_WITH_MESSAGE(context != nullptr, "context should not be NULL");
	CHECK_ARGS_WITH_MESSAGE(0 < num_data && num_data <= INT_MAX,
			"num_data must be 0 < num_data <= INT_MAX");
	CHECK_ARGS_WITH_MESSAGE(
			input_data != nullptr && LIBSAKURA_SYMBOL(IsAligned)(input_data),
			"invalid input_data");
	CHECK_ARGS_WITH_MESSAGE(
			output_data != nullptr && LIBSAKURA_SYMBOL(IsAligned)(output_data),
			"invalid output_data");
	CHECK_ARGS_WITH_MESSAGE(context->kernel_width == num_data,
			"using FFT for convolution. num_data must be equal to the one in the context");
	//assert(fftw_alignment_of((double *)input_data) == 0);
	//assert(fftw_alignment_of((double *)output_data) == 0);
	try {
		ConvolutionWithFFT(context, num_data, input_data, output_data);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DestroyConvolve1DContextFloat)(
LIBSAKURA_SYMBOL(Convolve1DContextFloat) *context) noexcept {
	if (context == nullptr) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	try {
		DestroyConvolve1DContextFloat(context);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}
