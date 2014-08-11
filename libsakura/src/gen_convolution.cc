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
#include <cassert>
#include <cstdlib>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"
#include "libsakura/logger.h"

namespace {
auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("Convolution");
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateConvolve1DContext)(
		size_t num_data, LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
		size_t kernel_width, bool use_fft,
		LIBSAKURA_SYMBOL(Convolve1DContext) **context) {
	if (context == nullptr) {
		LOG4CXX_ERROR(logger, "context should not be NULL");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (!(0 < num_data && num_data <= INT_MAX)) {
		LOG4CXX_ERROR(logger, "num_data must be '0 < num_data <= INT_MAX'");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (!(0 <= kernel_type
			&& kernel_type < LIBSAKURA_SYMBOL(Convolve1DKernelType_kNumType))) {
		LOG4CXX_ERROR(logger, "Invalid Kernel Type");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (!(0 < kernel_width)) {
		LOG4CXX_ERROR(logger, "kernel_width must be >0");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	try {
		auto convolutionop =
				::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetConvolutionImpl();
		convolutionop->CreateConvolve1DContext(num_data, kernel_type,
				kernel_width, use_fft, context);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const std::invalid_argument &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	} catch (const std::runtime_error &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kNG);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Convolve1D)(
LIBSAKURA_SYMBOL(Convolve1DContext) const *context, size_t num_data,
		float const input_data[/*num_data*/], float output_data[/*num_data*/]) {
	if (context == nullptr) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (!(0 < num_data && num_data <= INT_MAX)) {
		LOG4CXX_ERROR(logger, "num_data must be '0 < num_data <= INT_MAX'");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (!( LIBSAKURA_SYMBOL(IsAligned)(input_data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(output_data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	try {
		auto convolutionop =
				::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetConvolutionImpl();
		convolutionop->Convolve1D(context, num_data, input_data, output_data);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const std::invalid_argument &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(
LIBSAKURA_SYMBOL(Convolve1DContext) *context) {
	if (context == nullptr) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	try {
		auto convolutionop =
				::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetConvolutionImpl();
		convolutionop->DestroyConvolve1DContext(context);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}
