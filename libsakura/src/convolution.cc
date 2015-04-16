/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2014
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

inline void Create1DGaussianKernel(size_t num_kernel, bool use_fft,
		size_t kernel_width, float* kernel) {
	assert((2*num_kernel-1)>=0);
	size_t num_data_for_gauss = (use_fft) ? num_kernel : 2 * num_kernel - 1;
	assert(kernel_width != 0);
	float const reciprocal_of_denominator = 1.66510922231539551270632928979040
			/ kernel_width; // sqrt(log(16)) / kernel_width
	float const height = .939437278699651333772340328410
			/ static_cast<float>(kernel_width); // sqrt(8*log(2)/(2*M_PI)) / kernel_width
	float center =
			(num_data_for_gauss % 2 != 0) ?
					static_cast<float>(num_data_for_gauss - 1) / 2.f :
					static_cast<float>(num_data_for_gauss) / 2.f;
	size_t middle = (num_data_for_gauss) / 2;
	size_t loop_max = middle;
	kernel[0] = height;
	kernel[middle] = height
			* exp(
					-(center * reciprocal_of_denominator)
							* (center * reciprocal_of_denominator));
	size_t plus_one_for_odd = 0;
	if (num_data_for_gauss % 2 != 0) {
		plus_one_for_odd = 1;
		kernel[middle + 1] = kernel[middle];
	}
	for (size_t i = 1; i < loop_max; ++i) {
		float value = (i - center) * reciprocal_of_denominator;
		kernel[middle - i] = height * exp(-(value * value));
	}
	if (use_fft) {
		for (size_t i = 1; i < loop_max; ++i) {
			kernel[middle + i + plus_one_for_odd] = kernel[middle - i];
		}
	}

}

inline void Create1DKernel(size_t num_kernel, bool use_fft,
LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type, size_t kernel_width,
		float* kernel) {
	switch (kernel_type) {
	case LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian):
		Create1DGaussianKernel(num_kernel, use_fft, kernel_width, kernel);
		break;
	case LIBSAKURA_SYMBOL(Convolve1DKernelType_kBoxcar):
		break;
	case LIBSAKURA_SYMBOL(Convolve1DKernelType_kHanning):
		break;
	case LIBSAKURA_SYMBOL(Convolve1DKernelType_kHamming):
		break;
	default:
		assert(false);
		break;
	}
}

inline void ConvolutionWithoutFFT(size_t num_data, float const *input_data_arg,
		size_t num_kernel, float const *kernel_arg, float *output_data_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(input_data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(kernel_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(output_data_arg));
	auto input_data = AssumeAligned(input_data_arg);
	auto kernel = AssumeAligned(kernel_arg);
	auto output_data = AssumeAligned(output_data_arg);

	for (size_t i = 1; i < (num_data + 1); ++i) {
		float value = 0.0;
		size_t jmax = std::min(num_kernel, num_data - (i - 1));
		for (size_t j = 1; j < jmax + 1; ++j) {
			value += input_data[i - 1 + j - 1] * kernel[j - 1];
		}
		jmax = std::min(num_kernel, i);
		for (size_t j = 1; j < jmax; ++j) {
			value += input_data[i - 1 - j] * kernel[j];
		}
		output_data[i - 1] = value;
	}

}

inline void CreateConvolve1DContextFloat(size_t num_data,
LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type, size_t kernel_width,
bool use_fft, LIBSAKURA_SYMBOL(Convolve1DContextFloat)** context) {

	assert(context != nullptr);
	if (use_fft) {
		float *real_kernel_array = nullptr;
		std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> real_kernel_array_work(
				LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
						sizeof(float) * num_data, &real_kernel_array));
		assert(LIBSAKURA_SYMBOL(IsAligned)(real_kernel_array));
		float *real_array = nullptr;
		std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> real_array_work(
				LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
						sizeof(float) * num_data, &real_array));
		assert(LIBSAKURA_SYMBOL(IsAligned)(real_array));
		std::unique_ptr<fftwf_complex[], decltype(&FreeFFTArray)> fft_applied_complex_kernel(
				AllocateFFTArray(num_data / 2 + 1), FreeFFTArray);
		if (fft_applied_complex_kernel == nullptr) {
			throw std::bad_alloc();
		}
		std::unique_ptr<fftwf_complex[], decltype(&FreeFFTArray)> fft_applied_complex_input_data(
				AllocateFFTArray(num_data / 2 + 1), FreeFFTArray);
		if (fft_applied_complex_input_data == nullptr) {
			throw std::bad_alloc();
		}
		std::unique_ptr<fftwf_complex[], decltype(&FreeFFTArray)> multiplied_complex_data(
				AllocateFFTArray(num_data / 2 + 1), FreeFFTArray);
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
		Create1DKernel(num_data, use_fft, kernel_type, kernel_width,
				real_kernel_array);
		fftwf_execute(plan_real_to_complex_float_kernel);
		guard_for_fft_plan_kernel.CleanUpNow();
		std::unique_ptr<LIBSAKURA_SYMBOL(Convolve1DContextFloat),
		LIBSAKURA_PREFIX::Memory> work_context(
				static_cast<LIBSAKURA_SYMBOL(Convolve1DContextFloat)*>(LIBSAKURA_PREFIX::Memory::Allocate(
						sizeof(LIBSAKURA_SYMBOL(Convolve1DContextFloat)))),
				LIBSAKURA_PREFIX::Memory());
		if (work_context == nullptr) {
			throw std::bad_alloc();
		}
		work_context->fft_applied_complex_kernel =
				fft_applied_complex_kernel.release();
		work_context->real_array = real_array;
		work_context->real_array_work = real_array_work.release();
		work_context->real_kernel_array = nullptr;
		work_context->real_kernel_array_work = nullptr;
		work_context->plan_real_to_complex_float = plan_real_to_complex_float;
		guard_for_fft_plan.Disable();
		work_context->plan_complex_to_real_float = plan_complex_to_real_float;
		guard_for_ifft_plan.Disable();
		work_context->num_data = num_data;
		work_context->kernel_width = kernel_width;
		work_context->use_fft = use_fft;
		*context = work_context.release();
	} else {
		size_t threshold = kernel_width / 2;
		if (kernel_type == LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian)) {
			constexpr double six_sigma = 6.0 / sqrt(8.0 * log(2.0));
			threshold = kernel_width * six_sigma;
		}
		size_t num_kernel = threshold + 1;
		if (2 * num_kernel - 1 > num_data) {
			num_kernel = num_data;
		}
		float *real_kernel_array = nullptr;
		std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> real_kernel_array_work(
				LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
						sizeof(float) * num_kernel, &real_kernel_array));
		assert(LIBSAKURA_SYMBOL(IsAligned)(real_kernel_array));
		Create1DKernel(num_kernel, use_fft, kernel_type, kernel_width,
				real_kernel_array);
		std::unique_ptr<LIBSAKURA_SYMBOL(Convolve1DContextFloat),
		LIBSAKURA_PREFIX::Memory> work_context(
				static_cast<LIBSAKURA_SYMBOL(Convolve1DContextFloat)*>(LIBSAKURA_PREFIX::Memory::Allocate(
						sizeof(LIBSAKURA_SYMBOL(Convolve1DContextFloat)))),
				LIBSAKURA_PREFIX::Memory());
		if (work_context == nullptr) {
			throw std::bad_alloc();
		}
		work_context->fft_applied_complex_kernel = nullptr;
		work_context->plan_complex_to_real_float = nullptr;
		work_context->plan_real_to_complex_float = nullptr;
		work_context->real_array = nullptr;
		work_context->real_array_work = nullptr;
		work_context->real_kernel_array = real_kernel_array;
		work_context->real_kernel_array_work = real_kernel_array_work.release();
		work_context->num_data = num_data;
		work_context->kernel_width = num_kernel;
		work_context->use_fft = use_fft;
		*context = work_context.release();
	}
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
		std::unique_ptr<fftwf_complex[], decltype(&FreeFFTArray)> fft_applied_complex_input_data(
				AllocateFFTArray(num_data / 2 + 1), FreeFFTArray);
		if (fft_applied_complex_input_data == nullptr) {
			throw std::bad_alloc();
		}
		std::unique_ptr<fftwf_complex[], decltype(&FreeFFTArray)> multiplied_complex_data(
				AllocateFFTArray(num_data / 2 + 1), FreeFFTArray);
		if (multiplied_complex_data == nullptr) {
			throw std::bad_alloc();
		}
		fftwf_execute_dft_r2c(context->plan_real_to_complex_float,
				const_cast<float*>(input_data),
				fft_applied_complex_input_data.get());
		float scale = 1.0f / static_cast<float>(num_data);
		for (size_t i = 1; i < num_data / 2 + 1 + 1; ++i) {
			multiplied_complex_data[i - 1][0] =
					(context->fft_applied_complex_kernel[i - 1][0]
							* fft_applied_complex_input_data[i - 1][0]
							- context->fft_applied_complex_kernel[i - 1][1]
									* fft_applied_complex_input_data[i - 1][1])
							* scale;
			multiplied_complex_data[i - 1][1] =
					(context->fft_applied_complex_kernel[i - 1][0]
							* fft_applied_complex_input_data[i - 1][1]
							+ context->fft_applied_complex_kernel[i - 1][1]
									* fft_applied_complex_input_data[i - 1][0])
							* scale;
		}
		fftwf_execute_dft_c2r(context->plan_complex_to_real_float,
				multiplied_complex_data.get(), output_data);
	} else {
		ConvolutionWithoutFFT(num_data, const_cast<float*>(input_data),
				context->kernel_width, context->real_kernel_array, output_data);
	}
}

inline void DestroyConvolve1DContextFloat(
LIBSAKURA_SYMBOL(Convolve1DContextFloat)* context) {
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

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateConvolve1DContextFloat)(
		size_t num_data, LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
		size_t kernel_width, bool use_fft,
		LIBSAKURA_SYMBOL(Convolve1DContextFloat) **context) noexcept {
	if (context == nullptr) {
		LOG4CXX_ERROR(logger, "context should not be NULL");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (!(0 < num_data && num_data <= INT_MAX)) {
		LOG4CXX_ERROR(logger, "num_data must be '0 < num_data <= INT_MAX'");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (!(0 <= kernel_type
			&& kernel_type < LIBSAKURA_SYMBOL(Convolve1DKernelType_kNumElements))) {
		LOG4CXX_ERROR(logger, "Invalid Kernel Type");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (!(0 < kernel_width)) {
		LOG4CXX_ERROR(logger, "kernel_width must be >0");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	try {
		CreateConvolve1DContextFloat(num_data, kernel_type, kernel_width,
				use_fft, context);
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
