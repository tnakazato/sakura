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
#include <sstream>
#include <memory>
#include <utility>
#include <vector>

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
#include <libsakura/logger.h>
#include <libsakura/memory_manager.h>

#include "locator.tcc"
#include "interpolation_utils.tcc"
#include "interpolation.tcc"

namespace {
// a logger for this module
auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("interpolation");

template<class Interpolator, class XDataType, class YDataType>
void Interpolate1D(uint8_t polynomial_order, size_t num_base,
		XDataType const base_position[/*num_x_base*/], size_t num_array,
		YDataType const base_data[/*num_x_base*num_y*/],
		size_t num_interpolated,
		XDataType const interpolated_position[/*num_x_interpolated*/],
		YDataType interpolated_data[/*num_x_interpolated*num_y*/]) {
	assert(num_base > 0);
	assert(num_array > 0);
	assert(num_interpolated > 0);
	assert(base_position != nullptr);
	assert(base_data != nullptr);
	assert(interpolated_position != nullptr);
	assert(interpolated_data != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(base_position));
	assert(LIBSAKURA_SYMBOL(IsAligned)(base_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(interpolated_position));
	assert(LIBSAKURA_SYMBOL(IsAligned)(interpolated_data));

	if (num_base == 1) {
		// No need to interpolate, just substitute base_data
		// to all elements in y_interpolated
		Interpolator::SubstituteSingleBaseData(num_base, num_array,
				num_interpolated, base_data, interpolated_data);
		return;
	}

	std::vector<StorageAndAlignedPointer<XDataType> > xdatatype_holder(2);
	StorageAndAlignedPointer<YDataType> ydatatype_holder;
	GetAscendingArray<typename Interpolator::XDataReorderer, XDataType,
			XDataType>(num_base, base_position, 1, base_position,
			&xdatatype_holder[0]);
	GetAscendingArray<typename Interpolator::YDataReorderer, XDataType,
			YDataType>(num_base, base_position, num_array, base_data,
			&ydatatype_holder);
	GetAscendingArray<typename Interpolator::XDataReorderer, XDataType,
			XDataType>(num_interpolated, interpolated_position, 1,
			interpolated_position, &xdatatype_holder[1]);
	XDataType const *base_position_work = xdatatype_holder[0].pointer;
	YDataType const *base_data_work = ydatatype_holder.pointer;
	XDataType const *interpolated_position_work = xdatatype_holder[1].pointer;

	// Generate worker class
	Interpolator interpolator;

	// Perform 1-dimensional interpolation
	// Any preparation for interpolation should be done here
	interpolator.PrepareForInterpolation(polynomial_order, num_base, num_array,
			base_position_work, base_data_work);

	// Locate each element in x_base against x_interpolated
	StorageAndAlignedPointer<size_t> size_t_holder;
	AllocateAndAlign<size_t>(num_base, &size_t_holder);
	size_t *location_base = size_t_holder.pointer;
	size_t num_location_base = Locate<XDataType>(num_interpolated, num_base,
			interpolated_position_work, base_position_work, location_base);

	// Outside of x_base[0]
	Interpolator::SubstituteLeftMostData(location_base[0], num_base, num_array,
			num_interpolated, base_data_work, interpolated_data);

	// Between x_base[0] and x_base[num_x_base-1]
	size_t offset = 0;
	if (base_position_work[0] < interpolated_position_work[0]) {
		for (size_t i = 1; i < num_base; ++i) {
			if (base_position_work[offset + 1]
					< interpolated_position_work[0]) {
				offset++;
			} else {
				break;
			}
		}
	}
	interpolator.Interpolate1D(num_base, base_position_work, num_array,
			base_data_work, num_interpolated, interpolated_position_work,
			interpolated_data, num_location_base, location_base, offset);

	// Outside of x_base[num_x_base-1]
	Interpolator::SubstituteRightMostData(location_base[num_location_base - 1],
			num_base, num_array, num_interpolated, base_data_work,
			interpolated_data);

	// swap output array
	if (interpolated_position[0]
			> interpolated_position[num_interpolated - 1]) {
		Interpolator::SwapResult(num_array, num_interpolated,
				interpolated_data);
	}
}

template<class InterpolatorSet, class XDataType, class YDataType>
void ExecuteInterpolate1D(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_base,
		XDataType const base_position[], size_t num_array,
		YDataType const base_data[], size_t num_interpolated,
		XDataType const interpolated_position[],
		YDataType interpolated_data[]) {
	typedef void (*Interpolate1DFunc)(uint8_t, size_t, XDataType const *,
			size_t, YDataType const *, size_t, XDataType const *, YDataType *);
	Interpolate1DFunc func = [](uint8_t, size_t, XDataType const *,
			size_t, YDataType const *, size_t, XDataType const *, YDataType *) {
		throw LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	};
	switch (interpolation_method) {
	case LIBSAKURA_SYMBOL(InterpolationMethod_kNearest):
		func = Interpolate1D<typename InterpolatorSet::NearestInterpolator,
				XDataType, YDataType>;
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kLinear):
		func = Interpolate1D<typename InterpolatorSet::LinearInterpolator,
				XDataType, YDataType>;
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kPolynomial):
		if (polynomial_order == 0) {
			// This is special case: 0-th polynomial interpolation
			// acts like nearest interpolation
			func = Interpolate1D<typename InterpolatorSet::NearestInterpolator,
					XDataType, YDataType>;
		} else {
			func = Interpolate1D<
					typename InterpolatorSet::PolynomialInterpolator, XDataType,
					YDataType>;
		}
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kSpline):
		func = Interpolate1D<typename InterpolatorSet::SplineInterpolator,
				XDataType, YDataType>;
		break;
	default:
		// invalid interpolation method type
		break;
	}
	(*func)(polynomial_order, num_base, base_position, num_array, base_data,
			num_interpolated, interpolated_position, interpolated_data);
}

// basic check of arguments
bool CheckArguments(LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_interpolation_axis,
		double const base[], size_t num_array, float const data_base[],
		size_t num_interpolated, double const interpolated[],
		float const data_interpolated[], LIBSAKURA_SYMBOL(Status) *status) {

	bool process_data = true;

	// check interpolation_method
	if (interpolation_method != LIBSAKURA_SYMBOL(InterpolationMethod_kNearest)
			&& interpolation_method
					!= LIBSAKURA_SYMBOL(InterpolationMethod_kLinear)
			&& interpolation_method
					!= LIBSAKURA_SYMBOL(InterpolationMethod_kPolynomial)
			&& interpolation_method
					!= LIBSAKURA_SYMBOL(InterpolationMethod_kSpline)) {
		LOG4CXX_ERROR(logger, "Invalid interpolation method");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		process_data = false;
	}

	// num_base must be non-zero
	if (num_interpolation_axis == 0) {
		LOG4CXX_ERROR(logger, "num_base must be >0");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		process_data = false;
	}

	// no interpolation will be done
	if (num_interpolated == 0 || num_array == 0) {
		// Nothing to do
		LOG4CXX_INFO(logger, "Nothing has been done since num_interpolated is 0");
		*status = LIBSAKURA_SYMBOL(Status_kOK);
		process_data = false;
	}

	// input arrays are not aligned
	if (!LIBSAKURA_SYMBOL(IsAligned)(base)
			|| !LIBSAKURA_SYMBOL(IsAligned)(data_base)
			|| !LIBSAKURA_SYMBOL(IsAligned)(interpolated)
			|| !LIBSAKURA_SYMBOL(IsAligned)(data_interpolated)) {
		LOG4CXX_ERROR(logger, "input arrays are not aligned");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		process_data = false;
	}

	// input arrays are null
	if (base == nullptr || data_base == nullptr || interpolated == nullptr
			|| data_interpolated == nullptr) {
		LOG4CXX_ERROR(logger, "input arrays are null");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		process_data = false;
	}
	return process_data;
}

/**
 * Interpolate1DAlongColumn performs 1D interpolation along column based on base_position
 * and data_base. data_base is a serial array of column-major matrix data.
 * Its memory layout is assumed to be:
 *     data_base[0]     = data[0][0]
 *     data_base[1]     = data[1][0]
 *     data_base[2]     = data[2][0]
 *     ...
 *     data_base[N]     = data[N][0]
 *     data_base[N+1]   = data[1][0]
 *     ...
 *     data_base[N*M-1] = data[N][M]
 * where N and M correspond to num_x_base and num_y respectively.
 */
template<class XDataType, class YDataType>
void InterpolateXAxis(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_x_base,
		XDataType const x_base[/*num_x_base*/], size_t num_y,
		YDataType const data_base[/*num_x_base*num_y*/],
		size_t num_x_interpolated,
		XDataType const x_interpolated[/*num_x_interpolated*/],
		YDataType data_interpolated[/*num_x_interpolated*num_y*/]) {
	ExecuteInterpolate1D<XInterpolatorSet<XDataType, YDataType>, XDataType,
			YDataType>(interpolation_method, polynomial_order, num_x_base,
			x_base, num_y, data_base, num_x_interpolated, x_interpolated,
			data_interpolated);
}

/**
 * Interpolate1DAlongColumn performs 1D interpolation along row based on y_base
 * and data_base. data_base is a serial array of row-major matrix data.
 * Its memory layout is assumed to be:
 *     data_base[0]     = data[0][0]
 *     data_base[1]     = data[0][1]
 *     data_base[2]     = data[0][2]
 *     ...
 *     data_base[M]     = data[0][M]
 *     data_base[M+1]   = data[1][0]
 *     ...
 *     data_base[M*N-1] = data[N][M]
 * where N and M correspond to num_y_base and num_x respectively.
 */
template<class XDataType, class YDataType>
void InterpolateYAxis(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_y_base,
		XDataType const y_base[/*num_y_base*/], size_t num_x,
		YDataType const data_base[/*num_y_base*num_x*/],
		size_t num_y_interpolated,
		XDataType const y_interpolated[/*num_y_interpolated*/],
		YDataType data_interpolated[/*num_y_interpolated*num_x*/]) {
	ExecuteInterpolate1D<YInterpolatorSet<XDataType, YDataType>, XDataType,
			YDataType>(interpolation_method, polynomial_order, num_y_base,
			y_base, num_x, data_base, num_y_interpolated, y_interpolated,
			data_interpolated);
}

} /* anonymous namespace */

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InterpolateXAxisFloat)(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_x_base,
		double const x_base[/*num_x_base*/], size_t num_y,
		float const data_base[/*num_x_base*num_y*/], size_t num_x_interpolated,
		double const x_interpolated[/*num_x_interpolated*/],
		float data_interpolated[/*num_x_interpolated*num_y*/]) {

	// check arguments
	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Status_kOK);
	if (!CheckArguments(interpolation_method, polynomial_order, num_x_base,
			x_base, num_y, data_base, num_x_interpolated, x_interpolated,
			data_interpolated, &status)) {
		return status;
	}

	try {
		InterpolateXAxis(interpolation_method,
				polynomial_order, num_x_base, x_base, num_y, data_base,
				num_x_interpolated, x_interpolated, data_interpolated);
	} catch (const std::bad_alloc &e) {
		// failed to allocate memory
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const LIBSAKURA_SYMBOL(Status) &LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		// failed to allocate memory
		LOG4CXX_ERROR(logger, "Invalid interpolation type.");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	} catch (...) {
		// any exception is thrown during interpolation
		assert(false);
		LOG4CXX_ERROR(logger, "Aborted due to unknown error");
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return status;
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InterpolateYAxisFloat)(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_y_base,
		double const y_base[/*num_y_base*/], size_t num_x,
		float const data_base[/*num_y_base*num_x*/], size_t num_y_interpolated,
		double const y_interpolated[/*num_y_interpolated*/],
		float data_interpolated[/*num_y_interpolated*num_x*/]) {

	// check arguments
	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Status_kOK);
	if (!CheckArguments(interpolation_method, polynomial_order, num_y_base,
			y_base, num_x, data_base, num_y_interpolated, y_interpolated,
			data_interpolated, &status)) {
		return status;
	}

	try {
		InterpolateYAxis(interpolation_method,
				polynomial_order, num_y_base, y_base, num_x, data_base,
				num_y_interpolated, y_interpolated, data_interpolated);
	} catch (const std::bad_alloc &e) {
		// failed to allocate memory
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const LIBSAKURA_SYMBOL(Status) &LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		// failed to allocate memory
		LOG4CXX_ERROR(logger, "Invalid interpolation type.");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	} catch (...) {
		// any exception is thrown during interpolation
		assert(false);
		LOG4CXX_ERROR(logger, "Aborted due to unknown error");
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return status;
}
