/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2015
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
	bool use_fft;
	size_t num_data;
	size_t kernel_width;
	fftwf_plan plan_real_to_complex_float;
	fftwf_plan plan_complex_to_real_float;
	fftwf_complex *fft_applied_complex_kernel;
	float *real_array;
	void *real_array_work;
	float *real_kernel_array;
	void *real_kernel_array_work;
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
	assert((2 * num_kernel - 1) >= 0);
	assert(kernel_width != 0);
	double const reciprocal_of_denominator = 1.66510922231539551270632928979040
			/ static_cast<double>(kernel_width); // sqrt(log(16)) / kernel_width
	double const peak_value = .939437278699651333772340328410
			/ static_cast<double>(kernel_width); // sqrt(8*log(2)/(2*M_PI)) / kernel_width
	size_t peak_location = (num_kernel) / 2;

	for (size_t i = 0; i < num_kernel; ++i) {
		double value = (static_cast<double>(i)
				- static_cast<double>(peak_location))
				* reciprocal_of_denominator;
		kernel[i] = peak_value * exp(-(value * value));
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
inline void ConvolutionWithoutFFT(size_t num_data, float const *input_data_arg,
		size_t num_kernel, float const *kernel_arg, float *output_data_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(input_data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(kernel_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(output_data_arg));
	auto input_data = AssumeAligned(input_data_arg);
	auto kernel = AssumeAligned(kernel_arg);
	auto output_data = AssumeAligned(output_data_arg);

	size_t center_index = num_kernel / 2;
	size_t left_half = num_kernel / 2;
	size_t right_half = num_kernel - left_half - 1;
	for (size_t i = 0; i < num_data; ++i) {
		size_t lower_index = (i >= center_index) ? i - center_index : 0;
		size_t kernel_start = (i >= left_half) ? 0 : left_half - i;
		size_t kernel_end =
				(num_data - i - 1 >= right_half) ?
						num_kernel :
						num_kernel - (right_half - (num_data - i - 1));
		assert(kernel_end - kernel_start <= num_kernel);
		double value = 0.0;
		size_t num_convolve = kernel_end - kernel_start;
		for (size_t j = 0; j < num_convolve; ++j) {
			value += kernel[kernel_start + j] * input_data[lower_index + j];
		}
		output_data[i] = value;
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

inline void CreateConvolve1DContextFloat(size_t num_data, size_t num_kernel,
		float const kernel[],
		bool use_fft, LIBSAKURA_SYMBOL(Convolve1DContextFloat)** context) {

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
	work_context->use_fft = use_fft;
	work_context->num_data = num_data;

	if (use_fft) {
		assert(num_data == num_kernel);
		size_t num_fft_kernel = num_data / 2 + 1;
		float *real_kernel_array = nullptr;
		std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> real_kernel_array_work(
				LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
						sizeof(float) * num_data, &real_kernel_array));
		assert(LIBSAKURA_SYMBOL(IsAligned)(real_kernel_array));
		FlipData1D(num_data, kernel, real_kernel_array);
		float *real_array = nullptr;
		std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> real_array_work(
				LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
						sizeof(float) * num_data, &real_array));
		assert(LIBSAKURA_SYMBOL(IsAligned)(real_array));
		std::unique_ptr<fftwf_complex[], decltype(&FreeFFTArray)> fft_applied_complex_kernel(
				AllocateFFTArray(num_fft_kernel), FreeFFTArray);
		if (fft_applied_complex_kernel == nullptr) {
			throw std::bad_alloc();
		}
		std::unique_ptr<fftwf_complex[], decltype(&FreeFFTArray)> fft_applied_complex_input_data(
				AllocateFFTArray(num_fft_kernel), FreeFFTArray);
		if (fft_applied_complex_input_data == nullptr) {
			throw std::bad_alloc();
		}
		std::unique_ptr<fftwf_complex[], decltype(&FreeFFTArray)> multiplied_complex_data(
				AllocateFFTArray(num_fft_kernel), FreeFFTArray);
		if (multiplied_complex_data == nullptr) {
			throw std::bad_alloc();
		}
		fftwf_plan plan_real_to_complex_float_kernel = fftwf_plan_dft_r2c_1d(
				num_data, real_kernel_array, fft_applied_complex_kernel.get(),
				FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
		ScopeGuard guard_for_fft_plan_kernel([&]() {
			DestroyFFTPlan(plan_real_to_complex_float_kernel);
		});
		if (plan_real_to_complex_float_kernel == nullptr) {
			guard_for_fft_plan_kernel.Disable();
			throw std::bad_alloc();
		}
		fftwf_plan plan_real_to_complex_float = fftwf_plan_dft_r2c_1d(num_data,
				real_array, fft_applied_complex_input_data.get(),
				FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
		ScopeGuard guard_for_fft_plan([&]() {
			DestroyFFTPlan(plan_real_to_complex_float);
		});
		if (plan_real_to_complex_float == nullptr) {
			guard_for_fft_plan.Disable();
			throw std::bad_alloc();
		}
		fftwf_plan plan_complex_to_real_float = fftwf_plan_dft_c2r_1d(num_data,
				multiplied_complex_data.get(), real_array,
				FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
		ScopeGuard guard_for_ifft_plan([&]() {
			DestroyFFTPlan(plan_complex_to_real_float);
		});
		if (plan_complex_to_real_float == nullptr) {
			guard_for_ifft_plan.Disable();
			throw std::bad_alloc();
		}
		fftwf_execute(plan_real_to_complex_float_kernel);
		guard_for_fft_plan_kernel.CleanUpNow();

		work_context->kernel_width = num_kernel; //kernel_width;
		work_context->plan_real_to_complex_float = plan_real_to_complex_float;
		work_context->plan_complex_to_real_float = plan_complex_to_real_float;
		work_context->fft_applied_complex_kernel =
				fft_applied_complex_kernel.release();
		work_context->real_array = real_array;
		work_context->real_array_work = real_array_work.release();
		work_context->real_kernel_array = nullptr;
		work_context->real_kernel_array_work = nullptr;
		guard_for_fft_plan.Disable();
		guard_for_ifft_plan.Disable();

	} else {
		float *real_kernel_array = nullptr;
		std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> real_kernel_array_work(
				LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
						sizeof(float) * num_kernel, &real_kernel_array));
		assert(LIBSAKURA_SYMBOL(IsAligned)(real_kernel_array));
		for (size_t i = 0; i < num_kernel; ++i) {
			real_kernel_array[i] = kernel[i];
		}

		work_context->kernel_width = num_kernel;
		work_context->plan_real_to_complex_float = nullptr;
		work_context->plan_complex_to_real_float = nullptr;
		work_context->fft_applied_complex_kernel = nullptr;
		work_context->real_array = nullptr;
		work_context->real_array_work = nullptr;
		work_context->real_kernel_array = real_kernel_array;
		work_context->real_kernel_array_work = real_kernel_array_work.release();
	}
	*context = work_context.release();
}

inline void Convolve1DFloat(
LIBSAKURA_SYMBOL(Convolve1DContextFloat) const *context, size_t num_data,
		float const input_data_arg[/*num_data*/],
		float output_data_arg[/*num_data*/]) {
	assert(context->num_data == num_data);
	assert(LIBSAKURA_SYMBOL(IsAligned)(input_data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(output_data_arg));
	auto input_data = AssumeAligned(input_data_arg);
	auto output_data = AssumeAligned(output_data_arg);
	if (context->use_fft) {
		size_t num_fft_data = num_data / 2 + 1;
		std::unique_ptr<fftwf_complex[], decltype(&FreeFFTArray)> fft_applied_complex_input_data(
				AllocateFFTArray(num_fft_data), FreeFFTArray);
		if (fft_applied_complex_input_data == nullptr) {
			throw std::bad_alloc();
		}
		std::unique_ptr<fftwf_complex[], decltype(&FreeFFTArray)> multiplied_complex_data(
				AllocateFFTArray(num_fft_data), FreeFFTArray);
		if (multiplied_complex_data == nullptr) {
			throw std::bad_alloc();
		}
		fftwf_execute_dft_r2c(context->plan_real_to_complex_float,
				const_cast<float*>(input_data),
				fft_applied_complex_input_data.get());
		float scale = 1.0f / static_cast<float>(num_data);
		for (size_t i = 0; i < num_fft_data; ++i) {
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
	} else {
		ConvolutionWithoutFFT(num_data, input_data, context->kernel_width,
				context->real_kernel_array, output_data);
	}
}

inline void DestroyConvolve1DContextFloat(
LIBSAKURA_SYMBOL(Convolve1DContextFloat)* context) {
	if (context != nullptr) {
		DestroyFFTPlan(context->plan_real_to_complex_float);
		DestroyFFTPlan(context->plan_complex_to_real_float);
		FreeFFTArray(context->fft_applied_complex_kernel);
		if (context->real_kernel_array_work != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(context->real_kernel_array_work);
		}
		if (context->real_array_work != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(context->real_array_work);
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
	CHECK_ARGS(1 < num_kernel && num_kernel <= INT_MAX);
	CHECK_ARGS(kernel_width > 0);
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

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateConvolve1DContextFloat)(
		size_t num_data, size_t num_kernel, float const kernel[], bool use_fft,
		LIBSAKURA_SYMBOL(Convolve1DContextFloat) **context) noexcept {
	CHECK_ARGS_WITH_MESSAGE(0 < num_data && num_data <= INT_MAX,
			"num_data must satify '0 < num_data <= INT_MAX'");
	CHECK_ARGS_WITH_MESSAGE(0 < num_kernel && num_kernel <= INT_MAX,
			"num_kernel must satisfy '0 < num_data <= INT_MAX'");
	CHECK_ARGS_WITH_MESSAGE(kernel != nullptr, "kernel should not be NULL");
	CHECK_ARGS_WITH_MESSAGE(LIBSAKURA_SYMBOL(IsAligned)(kernel),
			"kernel should be aligned");
	if (use_fft) {
		CHECK_ARGS_WITH_MESSAGE(num_kernel == num_data,
				"num_kernel must be equal to num_data if use_fft is true");
	}
	CHECK_ARGS_WITH_MESSAGE(context != nullptr, "context should not be NULL");
	try {
		CreateConvolve1DContextFloat(num_data, num_kernel, kernel, use_fft,
				context);
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
LIBSAKURA_SYMBOL(Convolve1DContextFloat) const *context, size_t num_data,
		float const input_data[/*num_data*/], float output_data[/*num_data*/])
				noexcept {
	//TODO
	//CHECK_ARGS(context == nullptr);
	if (context == nullptr) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (0 >= num_data || num_data > INT_MAX) {
		LOG4CXX_ERROR(logger, "num_data must be '0 < num_data <= INT_MAX'");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (LIBSAKURA_SYMBOL(IsAligned)(input_data) == false) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (LIBSAKURA_SYMBOL(IsAligned)(output_data) == false) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	//assert(fftw_alignment_of((double *)input_data) == 0);
	//assert(fftw_alignment_of((double *)output_data) == 0);
	try {
		Convolve1DFloat(context, num_data, input_data, output_data);
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
