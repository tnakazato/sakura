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
 * @file
 *  Sakura main header file.
 *
 * sakura.h
 *
 *  Created on: 2013/02/20
 *      Author: kohji
 */

#ifndef LIBSAKURA_LIBSAKURA_SAKURA_H_
#define LIBSAKURA_LIBSAKURA_SAKURA_H_

#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>
#include <sys/types.h>

#include <libsakura/config.h>

#if defined(__GNUC__) || defined(__GNUG__)
#	define LIBSAKURA_WARN_UNUSED_RESULT __attribute__((warn_unused_result))
#else
#	define LIBSAKURA_WARN_UNUSED_RESULT /* Don't ignore result value */
#endif

#define LIBSAKURA_NOEXCEPT /* noexcept */

#ifdef __cplusplus
extern "C" {

#if __cplusplus >= 201103L
# undef LIBSAKURA_NOEXCEPT
# define LIBSAKURA_NOEXCEPT noexcept
#endif

#endif

/**
 * @brief A result of function call.
 */
typedef enum {
	/**
	 * @brief Success or normal end
	 */LIBSAKURA_SYMBOL(Status_kOK) = 0,
	/**
	 * @brief Failure or abnormal end
	 */LIBSAKURA_SYMBOL(Status_kNG) = 1,
	/**
	 * @brief Illegal argument(s)
	 *
	 * This includes a violation of the must-be-aligned constraint.
	 */LIBSAKURA_SYMBOL(Status_kInvalidArgument) = 2,
	/**
	 * @brief No memory
	 */LIBSAKURA_SYMBOL(Status_kNoMemory) = 3,
	/**
	 * @brief Unknown error
	 */LIBSAKURA_SYMBOL(Status_kUnknownError) = 99
}LIBSAKURA_SYMBOL (Status);

/**
 * @brief A type of the allocator function used by Sakura Library.
 *
 * Implementation of the function of this type must be reentrant.
 *
 * @note
 * It must return a valid pointer to a memory region of size 0 if 0 is passed as @a size parameter.
 *
 * @param[in] size	Size of required memory in bytes
 * @return Allocated memory. NULL if failed to allocate.
 *
 * MT-safe
 */
typedef void *(*LIBSAKURA_SYMBOL(UserAllocator))(size_t size);

/**
 * @brief A type of the deallocator function used by Sakura Library.
 *
 * Implementation of the function of this type must be reentrant.
 *
 * @note
 * It must do nothing if NULL is passed as @a pointer parameter.
 *
 * @param[in] pointer	NULL or an address to be released which was allocated by the allocator of type @ref sakura_UserAllocator .
 *
 * MT-safe
 */
typedef void (*LIBSAKURA_SYMBOL(UserDeallocator))(void *pointer);

/**
 * @brief Initializes Sakura Library.
 *
 * You must initialize libsakura by calling this function before calling any other function of Sakura Library.
 *
 * Without calling @ref sakura_CleanUp() , don't call this function again.
 *
 * @param[in]	allocator	An allocator which is used when Sakura Library needs to allocate memory dynamically. posix_memalign(3) is used if NULL is provided. See @ref sakura_UserAllocator .
 * @param[in]	deallocator	A deallocator which is used when Sakura Library needs to free dynamically allocated memory. free(3) is used if NULL is provided. See @ref sakura_UserDeallocator .
 * @return Only when @ref sakura_Status_kOK is returned, you can use Sakura Library.
 *
 * MT-unsafe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Initialize)(
LIBSAKURA_SYMBOL(UserAllocator) allocator,
LIBSAKURA_SYMBOL(UserDeallocator) deallocator)
		LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @brief Cleans up Sakura Library.
 *
 * When you call this function, no function of Sakura Library must be running.
 *
 * MT-unsafe
 */
void LIBSAKURA_SYMBOL(CleanUp)() LIBSAKURA_NOEXCEPT;

/*
 * memory alignment(for SIMD)
 */
/**
 * @brief Checks if @a ptr points the aligned address Sakura Library requires.
 *
 * @param[in] ptr An address to be checked. NULL is allowed.
 * @return true if the address is aligned, otherwise false
 *
 * MT-safe
 */

bool LIBSAKURA_SYMBOL(IsAligned)(void const *ptr) LIBSAKURA_NOEXCEPT;

/**
 * @brief Returns an alignment that Sakura Library expects for arrays on which vector operations are performed.
 *
 * @return An expected alignment for arrays marked as must-be-aligned.
 *
 * MT-safe
 */
size_t LIBSAKURA_SYMBOL (GetAlignment)() LIBSAKURA_NOEXCEPT;

/**
 * @brief Returns an aligned address close to @a arena by adding 0 or minimum offset.
 *
 * It returns @a arena if @a arena is already aligned.
 *
 * @param[in] size_of_arena Size of the memory region pointed by @a arena
 * @param[in] arena Start address of a memory region
 * @param[in] size_required Required size after alignment
 * @return Aligned address if at least @a size_required bytes are available in @a arena after alignment,
 * otherwise NULL.
 *
 * MT-safe
 */
void *LIBSAKURA_SYMBOL(AlignAny)(size_t size_of_arena, void *arena,
		size_t size_required) LIBSAKURA_NOEXCEPT;
/**
 * @copybrief sakura_AlignAny
 *
 * It returns @a arena if @a arena is already aligned.
 *
 * @param[in] elements_in_arena The number of elements in @a arena , not a size in bytes.
 * @param[in] arena Start address of an array
 * @param[in] elements_required Required number of elements after alignment
 * @return Aligned address if at least @a elements_required are available in @a arena after alignment,
 * otherwise NULL.
 *
 * MT-safe
 */
float *LIBSAKURA_SYMBOL(AlignFloat)(size_t elements_in_arena, float *arena,
		size_t elements_required) LIBSAKURA_NOEXCEPT;

/**
 * @copydoc sakura_AlignFloat()
 */
double *LIBSAKURA_SYMBOL(AlignDouble)(size_t elements_in_arena, double *arena,
		size_t elements_required) LIBSAKURA_NOEXCEPT;

/**
 * @brief A structure in which the result of @ref sakura_ComputeStatisticsFloat and @ref sakura_ComputeAccurateStatisticsFloat is stored.
 *
 * You can also calculate following statistics from the members of this struct if count > 0:
 *  - mean = sum / count
 *  - rms = sqrt(square_sum / count)
 *  - variance = abs(square_sum / count - mean * mean)
 *  - stddev = sqrt(variance)
 */
typedef struct {
	/**
	 * number of valid data
	 */
	size_t count;
	/**
	 * sum of valid data
	 */
	double sum;
	/**
	 * sum of squared valid data
	 */
	double square_sum;
	/**
	 * min value of valid data. NaN if no valid data.
	 */
	float min;
	/**
	 * max value of valid data. NaN if no valid data.
	 */
	float max;
	/**
	 * index for one of min value. -1 if there is no valid data.
	 */
	ssize_t index_of_min;
	/**
	 * index for one of max value. -1 if there is no valid data.
	 */
	ssize_t index_of_max;
}LIBSAKURA_SYMBOL(StatisticsResultFloat);

/**
 * @brief Computes statistics. Refer to
 * @ref sakura_StatisticsResultFloat to see what kind of statistics are computed.
 *
 * @param[in] num_data The number of elements in @a data and @a is_valid . @a num_data <= INT32_MAX
 * @param[in] data Data. If corresponding element in @a is_valid is true, the element in @a data must not be Inf nor NaN.
 * <br/>must-be-aligned
 * @param[in] is_valid Masks of @a data. If a value of element is false,
 * the corresponding element in @a data is ignored.
 * <br/>must-be-aligned
 * @param[out] result An address where the result should be stored. Some fields may be set to NaN if it is impossible to figure out.
 * If there is more than one occurrences of min or max value, it is undefined which index of the occurrences is selected for @a index_of_min or @a index_of_max.
 * @return Status code
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ComputeStatisticsFloat)(
		size_t num_data, float const data[], bool const is_valid[],
		LIBSAKURA_SYMBOL(StatisticsResultFloat) *result) LIBSAKURA_NOEXCEPT;

/**
 * @copybrief sakura_ComputeStatisticsFloat
 * @copydetails sakura_ComputeStatisticsFloat
 *
 * The result of this function is more accurate than that of @ref sakura_ComputeStatisticsFloat if
 * num_data is large. This function is slower than @ref sakura_ComputeStatisticsFloat .
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ComputeAccurateStatisticsFloat)(
		size_t num_data, float const data[], bool const is_valid[],
		LIBSAKURA_SYMBOL(StatisticsResultFloat) *result) LIBSAKURA_NOEXCEPT;

/**
 * @brief Sorts only valid data in ascending order.
 *
 * @param[in] num_data The number of elements in @a data and @a is_valid .
 * @param[in] is_valid Masks of @a data. If a value of element is false,
 * the corresponding element in @a data is ignored.
 * @param[in,out] data Data to be sorted. Since data is sorted in place, contents of this array are not preserved.
 * If corresponding element in @a is_valid is true, the element in @a data must not be Inf nor NaN.
 * @param[out] new_num_data The number of sorted elements that don't include invalid data( <= @a num_data ) is stored here.
 * @return Status code
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SortValidValuesDenselyFloat)(
		size_t num_data, bool const is_valid[], float data[],
		size_t *new_num_data) LIBSAKURA_NOEXCEPT;

/**
 * @brief Computes median absolute deviation.
 *
 * This function applies abs(@a data[i] - median) for each element in @a data and
 * stores results to @a new_data. Then it sorts @a new_data in ascending order.
 * @a data must be sorted in advance using e.g. @ref sakura_SortValidValuesDenselyFloat .
 * The median is a center value if @a num_data is odd. Otherwise, the median is
 * an average of two center values.
 * The caller is responsible to take a median value from @a new_data as
 * a median absolute deviation as below.
 *
 *  mad = (@a new_data[@a num_data/2] + @a new_data[@a num_data/2 - @a num_data%2]) / 2 if @a num_data > 0
 *
 * @param[in] num_data	The number of elements in @a data and @a new_data .
 * @param[in] data	Each value in data must not be Inf nor NaN. @a data must be sorted.
 * <br/>must-be-aligned
 * @param[out] new_data	An array where the results are stored. @a new_data may points @a data.
 * <br/>must-be-aligned
 * @return Status code
 *
 * MT-safe
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ComputeMedianAbsoluteDeviationFloat)(
		size_t num_data, float const data[], float new_data[])
				LIBSAKURA_NOEXCEPT;

/**
 * @brief Grids data with convolution.
 *
 * This function plot @a value from @a start_spectrum to @a end_spectrum whose location is at (@a x , @a y)
 * onto @a grid as points which have width represented by @a convolution_table.
 * The results of 2 * @a support pixels from the edges of the grids are not reliable due to
 * implementation to gain speed.
 * Polarizations and channels are mapped by @a polarization_map and @a channel_map when gridding.
 * All floating-point values passed to this function must not be NaN nor +-Inf.
 *
 * @param[in] num_spectra	The number of spectra. It should be  0 <= @a start_spectrum <= @a end_spectrum <= @a num_spectra .
 * @param[in] start_spectrum	Starting index of spectrum to be processed.
 * @param[in] end_spectrum	An index next to the last index of spectrum to be processed.
 * @param[in] spectrum_mask	Masks to each spectrum. The number of elements must be @a num_spectra.
 * If a value of element is false, the corresponding spectrum is ignored.
 * <br/>must-be-aligned
 * @param[in] x	X position of the point projected to 2D plane of the grid.
 * The number of elements must be @a num_spectra . Each element must be INT32_MIN < @a x[i] < INT32_MAX .
 * <br/>must-be-aligned
 * @param[in] y	Y position of the point projected to 2D plane of the grid.
 * The number of elements must be @a num_spectra . Each element must be INT32_MIN < @a y[i] < INT32_MAX .
 * <br/>must-be-aligned
 * @param[in] support	A half width(the number of pixels from center to edge) of convolution kernel on the grid plane of size @a width x @a height. It must be 0 < @a support <= 256 .<br/>
 * Note that this function may cause stack overflow because it allocates memory of a size proportional to @a support * @a sampling on stack
 * when @a support * @a sampling is too large.
 * @param[in] sampling	A resolution of convolution kernel(samples/grid aka samples/pixel). It must be 0 < @a sampling <= INT32_MAX .
 * @param[in] num_polarizations	The number of polarizations. It must be 0 < @a num_polarizations <= INT32_MAX .
 * @param[in] polarization_map	The number of elements must be @a num_polarizations. Each element must be in range [0,@a num_polarizations_for_grid).
 * <br/>must-be-aligned
 * @param[in] num_channels	The number of channels. It must be 0 < @a num_channels <= INT32_MAX .
 * @param[in] channel_map	The number of elements must be @a num_channels. Each element must be in range [0,@a num_channels_for_grid).
 * <br/>must-be-aligned
 * @param[in] mask	Masks to each @a value . Its memory layout should be [@a num_spectra][@a num_polarizations][@a num_channels].
 * If a value of element is false, the corresponding combination of the spectrum, polarization and channel is ignored.
 * <br/>must-be-aligned
 * @param[in] value	Values to be gridded. Its memory layout should be [@a num_spectra][@a num_polarizations][@a num_channels].
 * <br/>must-be-aligned
 * @param[in] weight Weights for @a value. Its memory layout should be [@a num_spectra][@a num_channels].
 * <br/>must-be-aligned
 * @param[in] weight_only	True if you want to get a grid of weight itself rather than production of @a value and @a weight. Otherwise false.
 * @param[in] num_convolution_table	The number of elements of @a convolution_table. It should be ceil(sqrt(2.)*(@a support+1)*@a sampling) <= @a num_convolution_table <= INT32_MAX / 32 .
 * @param[in] convolution_table	An array which represents convolution kernel. The number of elements must be @a num_convolution_table. The first element corresponds to center of the point.
 * <br/>must-be-aligned
 * @param[in] num_polarizations_for_grid	The number of polarizations on the grid. It should be 0 < @a num_polarizations_for_grid <= INT32_MAX .
 * @param[in] num_channels_for_grid	The number of channels on the grid. It should be 0 < @a num_channels_for_grid <= INT32_MAX .
 * @param[in] width	Width of the grid . It should be 0 < @a width <= INT32_MAX .
 * @param[in] height Height of the grid . It should be 0 < @a height <= INT32_MAX .
 * @param[out] weight_sum	Sum of weights. Its memory layout should be [@a num_polarizations_for_grid][@a num_channels_for_grid].
 * <br/>must-be-aligned
 * @param[out] weight_of_grid	Weight for each grid. Its memory layout should be [@a height][@a width][@a num_polarizations_for_grid][@a num_channels_for_grid].
 * <br/>must-be-aligned
 * @param[out] grid	The resulting grid. Its memory layout should be [@a height][@a width][@a num_polarizations_for_grid][@a num_channels_for_grid].
 * <br/>must-be-aligned
 * @return Status code
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GridConvolvingFloat)(
		size_t num_spectra, size_t start_spectrum, size_t end_spectrum,
		bool const spectrum_mask[/*num_spectra*/],
		double const x[/*num_spectra*/], double const y[/*num_spectra*/],
		size_t support, size_t sampling, size_t num_polarizations,
		uint32_t const polarization_map[/*num_polarizations*/],
		size_t num_channels, uint32_t const channel_map[/*num_channels*/],
		bool const mask/*[num_spectra][num_polarizations]*/[/*num_channels*/],
		float const value/*[num_spectra][num_polarizations]*/[/*num_channels*/],
		float const weight/*[num_spectra]*/[/*num_channels*/], bool weight_only,
		size_t num_convolution_table/*= ceil(sqrt(2.)*(support+1)*sampling)*/,
		float const convolution_table[/*num_convolution_table*/],
		size_t num_polarizations_for_grid, size_t num_channels_for_grid,
		size_t width, size_t height,
		double weight_sum/*[num_polarizations_for_grid]*/[/*num_channels_for_grid*/],
		float weight_of_grid/*[height][width][num_polarizations_for_grid]*/[/*num_channels_for_grid*/],
		float grid/*[height][width][num_polarizations_for_grid]*/[/*num_channels_for_grid*/])
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @brief Sets true if the values in an input array (@a data ) are in any of specified range (inclusive).
 * @details Elements of the output array are set to true if the corresponding element in the input array
 * is in range of upper and lower boundary pairs,
 * @par
 * @a lower_bound[k] <= @a data[i] <= @a upper_bound[k] ,
 *
 * otherwise they are set to false.
 *
 * The function takes more than one upper and lower boundary pairs as arrays,
 * @a lower_bounds and @a upper_bounds.@n
 *
 * @note
 * No evaluation is done when the data array is zero length, i.e., @a num_data = 0.@n
 * All elements in @a result are set to false when no condition is given, i.e., @a num_condition = 0.
 *
 * @param[in] num_data The number of elements in the arrays, @a data and @a result .
 * @param[in] data An input array of size, @a num_data.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] num_condition The number of elements in the arrays, @a lower_bounds
 * and @a upper_bounds.
 * @param[in] lower_bounds The input array of size, @a num_condition.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] upper_bounds The input array of size, @a num_condition.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[out] result The output array of size, @a num_data.
 * @n must-be-aligned
 * @return Status code
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfInRangesInclusiveFloat)(
		size_t num_data, float const data[/*num_data*/], size_t num_condition,
		float const lower_bounds[/*num_condition*/],
		float const upper_bounds[/*num_condition*/], bool result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @copybrief sakura_SetTrueIfInRangesInclusiveFloat
 * @copydetails sakura_SetTrueIfInRangesInclusiveFloat
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfInRangesInclusiveInt)(
		size_t num_data, int const data[/*num_data*/], size_t num_condition,
		int const lower_bounds[/*num_condition*/],
		int const upper_bounds[/*num_condition*/], bool result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @brief Sets true if the values in an input array (@a data ) are in any of specified range (exclusive).
 * @details Elements of the output array are set to true if the corresponding element in the input array
 * is in range of upper and lower boundary pairs,
 * @par
 * @a lower_bound[k] < @a data[i] < @a upper_bound[k] ,
 *
 * otherwise they are set to false.
 *
 * The function takes more than one upper and lower boundary pairs as arrays,
 * @a lower_bounds and @a upper_bounds.@n
 *
 * @note
 * No evaluation is done when the data array is zero length, i.e., @a num_data = 0.@n
 * All elements in @a result are set to false when no condition is given, i.e., @a num_condition = 0.
 *
 * @param[in] num_data The number of elements in the arrays, @a data
 * and @a result
 * @param[in] data An input array of size, @a num_data.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] num_condition The number of elements in the arrays, @a lower_bounds
 * and @a upper_bounds.
 * @param[in] lower_bounds The input array of size, @a num_condition.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] upper_bounds The input array of size, @a num_condition.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[out] result The output array of size, @a num_data.
 * @n must-be-aligned
 * @return Status code
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfInRangesExclusiveFloat)(
		size_t num_data, float const data[/*num_data*/], size_t num_condition,
		float const lower_bounds[/*num_condition*/],
		float const upper_bounds[/*num_condition*/], bool result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @copybrief sakura_SetTrueIfInRangesExclusiveFloat
 * @copydetails sakura_SetTrueIfInRangesExclusiveFloat
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfInRangesExclusiveInt)(
		size_t num_data, int const data[/*num_data*/], size_t num_condition,
		int const lower_bounds[/*num_condition*/],
		int const upper_bounds[/*num_condition*/], bool result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @brief Sets true if the values in the input array (@a data ) are greater than a threshold.
 * @details Elements of the output array are set to true if the corresponding element in the input array
 * is greater than a @a threshold,
 * @par
 * @a data[i] > @a threshold ,
 *
 * otherwise they are set to false.
 *
 * @note
 * No evaluation is done when the data array is zero length, i.e., @a num_data = 0.@n
 *
 * @param[in] num_data The number of elements in the arrays, @a data
 * and @a result
 * @param[in] data An input array of size, @a num_data.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] threshold The threshold of evaluation.
 * In case the parameter is floating-point type, the value should not be Inf nor NaN.
 * @param[out] result The output array of size, @a num_data.
 * @n must-be-aligned
 * @return Status code
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanFloat)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

/**
 * @copybrief sakura_SetTrueIfGreaterThanFloat
 * @copydetails sakura_SetTrueIfGreaterThanFloat
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

/**
 * @brief Sets true if the values in the input array (@a data ) are greater than or equal to a threshold.
 * @details Elements of the output array are set to true if the corresponding element in the input array
 * is greater than or equals to a @a threshold,
 * @par
 * @a data[i] >= @a threshold ,
 *
 * otherwise they are set to false.
 *
 * @note
 * No evaluation is done when the data array is zero length, i.e., @a num_data = 0.@n
 *
 * @param[in] num_data The number of elements in the arrays, @a data
 * and @a result
 * @param[in] data An input array of size, @a num_data.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] threshold The threshold of evaluation.
 * In case the parameter is floating-point type, the value should not be Inf nor NaN.
 * @param[out] result The output array of size, @a num_data.
 * @n must-be-aligned
 * @return Status code
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanOrEqualsFloat)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

/**
 * @copybrief sakura_SetTrueIfGreaterThanOrEqualsFloat
 * @copydetails sakura_SetTrueIfGreaterThanOrEqualsFloat
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanOrEqualsInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

/**
 * @brief Sets true if the values in the input array (@a data ) are less than a threshold.
 * @details Elements of the output array are set to true if the corresponding element in the input array
 * is less than a @a threshold,
 * @par
 * @a data[i] < @a threshold ,
 *
 * otherwise they are set to false.
 *
 * @note
 * No evaluation is done when the data array is zero length, i.e., @a num_data = 0.@n
 *
 * @param[in] num_data The number of elements in the arrays, @a data
 * and @a result
 * @param[in] data An input array of size, @a num_data.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] threshold The threshold of evaluation.
 * In case the parameter is floating-point type, the value should not be Inf nor NaN.
 * @param[out] result The output array of size, @a num_data.
 * @n must-be-aligned
 * @return Status code
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanFloat)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

/**
 * @copybrief sakura_SetTrueIfLessThanFloat
 * @copydetails sakura_SetTrueIfLessThanFloat
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

/**
 * @brief Sets true if the values in the input array (@a data ) are less than or equal to a threshold.
 * @details Elements of the output array are set to true if the corresponding element in the input array
 * is less than or equals to a @a threshold,
 * @par
 * @a data[i] <= @a threshold ,
 *
 * otherwise they are set to false.
 *
 * @note
 * No evaluation is done when the data array is zero length, i.e., @a num_data = 0.@n
 *
 * @param[in] num_data The number of elements in the arrays, @a data
 * and @a result
 * @param[in] data An input array of size, @a num_data.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] threshold The threshold of evaluation.
 * In case the parameter is floating-point type, the value should not be Inf nor NaN.
 * @param[out] result The output array of size, @a num_data.
 * @n must-be-aligned
 * @return Status code
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanOrEqualsFloat)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

/**
 * @copybrief sakura_SetTrueIfLessThanOrEqualsFloat
 * @copydetails sakura_SetTrueIfLessThanOrEqualsFloat
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanOrEqualsInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

/**
 * @brief Sets false if the values in the input array (@a data ) are either NaN or infinity.
 * @details Elements of the output array are set to false if the corresponding element in the input array
 * is not a number (NaN) or infinity. Otherwise, they are set to true.
 *
 * @note No operation is done when the data array is zero length, i.e., @a num_data = 0.
 *
 * @param[in] num_data The number of elements in the arrays, @a data
 * and @a result
 * @param[in] data The input array of of size, @a num_data.
 * @n must-be-aligned
 * @param[out] result The output array of of size, @a num_data.
 * @n must-be-aligned
 * @return Status code
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetFalseIfNanOrInfFloat)(
		size_t num_data, float const data[/*num_data*/],
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

/**
 * @brief Convert an input array to a boolean array.
 * @details Returns true if the corresponding element in input array != 0.
 * @note No operation is done when the data array is zero length, i.e., @a num_data = 0.
 *
 * @param[in] num_data The number of elements in the arrays, @a data
 * and @a result
 * @param[in] data The input array of of size, @a num_data.
 * @n must-be-aligned
 * @param[out] result The output array of of size, @a num_data.
 * @n must-be-aligned
 * @return Status code
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint8ToBool)(size_t num_data,
		uint8_t const data[/*num_data*/], bool result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @copybrief sakura_Uint8ToBool
 * @copydetails sakura_Uint8ToBool
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint32ToBool)(size_t num_data,
		uint32_t const data[/*num_data*/], bool result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @brief Invert a boolean array
 * @note No operation is done when the data array is zero length, i.e., @a num_data = 0.
 *
 * @param[in] num_data The number of elements in the arrays, @a data
 * and @a result
 * @param[in] data The input array of of size, @a num_data.
 * @n must-be-aligned
 * @param[out] result The output array of of size, @a num_data.
 * The pointer of @a out is allowed to be equal to
 * that of @a in (@a result == @a data), indicating in-place operation.
 * @n must-be-aligned
 * @return Status code
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InvertBool)(size_t num_data,
		bool const data[/*num_data*/], bool result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @brief Invoke bit operation AND between a bit mask and an array.
 * @details Invokes the following bit operation to @a i- th element of @a result :
 * @code
 * result [i] = edit_mask[i] ? (data[i] & bit_mask) : data[i]
 * @endcode
 *
 * @note
 * No operation is done when the data array is zero length, i.e., @a num_data = 0.
 *
 * @note
 * The function can also be used to invoke material nonimplication of @a data and @a bit_mask .
 * Input the complement of @a bit_mask (~@a bit_mask ) to invoke material nonimplication.
 * For details of material nonimplication, see, e.g.,@n
 * http://en.wikipedia.org/wiki/Truth_table
 *
 * @param[in] bit_mask A bit mask. The bit operation is invoked
 * between this value and the array, @a data.
 * @param[in] num_data The number of elements in the arrays, @a data,
 * @a edit_mask, and @a result.
 * @param[in] data An input array of size, @a num_data. The bit operation
 * is invoked between this array and @a bit_mask.@n
 * must-be-aligned
 * @param[in] edit_mask A boolean mask array of size, @a num_data. The bit operation
 * is skipped for the elements with the value, false.@n
 * must-be-aligned
 * @param[out] result The output array of size, @a num_data. It stores the result
 * of the bit operation between @a bit_mask and @a data. The bit operation is skipped
 * and the value in array, @a data, is adopted for the elements where corresponding
 * elements in @a edit_mask is false. The pointer of @a out is allowed to be equal to
 * that of @a in (@a result == @a data), indicating in-place operation.@n
 * must-be-aligned
 * @return Status code
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseAndUint8)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint8_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;
/**
 * @copybrief sakura_OperateBitwiseAndUint8
 * @copydetails sakura_OperateBitwiseAndUint8
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseAndUint32)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint32_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @brief Invoke bit operation, Converse Nonimplication, between a bit mask and an array.
 * @details Invokes the following bit operation to the @a i- th element of @a result :
 * @code
 * result [i] = edit_mask[i] ? (~data[i] & bit_mask) : data[i]
 * @endcode
 *
 * @note
 * No operation is done when the data array is zero length, i.e., @a num_data = 0.
 *
 * @note
 * The function can also be used to invoke bitwise NOR operation of @a data and @a bit_mask .
 * Input the complement of @a bit_mask (~@a bit_mask ) to invoke bitwise NOR operation.
 * For details of bitwise NOR operation, see, e.g.,@n
 * http://en.wikipedia.org/wiki/Truth_table
 *
 * @param[in] bit_mask A bit mask. The bit operation is invoked
 * between this value and the array, @a data.
 * @param[in] num_data The number of elements in the arrays, @a data,
 * @a edit_mask, and @a result.
 * @param[in] data An input array of size, @a num_data. The bit operation
 * is invoked between this array and @a bit_mask.@n
 * must-be-aligned
 * @param[in] edit_mask A boolean mask array of size, @a num_data. The bit operation
 * is skipped for the elements with the value, false.@n
 * must-be-aligned
 * @param[out] result The output array of size, @a num_data. It stores the result
 * of the bit operation between @a bit_mask and @a data. The bit operation is skipped
 * and the value in array, @a data, is adopted for the elements where corresponding
 * elements in @a edit_mask is false. The pointer of @a out is allowed to be equal to
 * that of @a in (@a result == @a data), indicating in-place operation.@n
 * must-be-aligned
 * @return Status code
 *
 * MT-safe
 *
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseConverseNonImplicationUint8)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint8_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;
/**
 * @copybrief sakura_OperateBitwiseConverseNonImplicationUint8
 * @copydetails sakura_OperateBitwiseConverseNonImplicationUint8
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseConverseNonImplicationUint32)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint32_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @brief Invoke bit operation, Material Implication, between a bit mask and an array.
 * @details Invokes the following bit operation to the @a i- th element of @a result :
 * @code
 * result [i] = edit_mask[i] ? (~data[i] | bit_mask) : data[i]
 * @endcode
 *
 * @note
 * No operation is done when the data array is zero length, i.e., @a num_data = 0.
 *
 * @note
 * The function can also be used to invoke bitwise NAND operation of @a data and @a bit_mask .
 * Input the complement of @a bit_mask (~@a bit_mask ) to invoke bitwise NAND operation.
 * For details of bitwise NAND operation, see, e.g.,@n
 * http://en.wikipedia.org/wiki/Truth_table
 *
 * @param[in] bit_mask A bit mask. The bit operation is invoked
 * between this value and the array, @a data.
 * @param[in] num_data The number of elements in the arrays, @a data,
 * @a edit_mask, and @a result.
 * @param[in] data An input array of size, @a num_data. The bit operation
 * is invoked between this array and @a bit_mask.@n
 * must-be-aligned
 * @param[in] edit_mask A boolean mask array of size, @a num_data. The bit operation
 * is skipped for the elements with the value, false.@n
 * must-be-aligned
 * @param[out] result The output array of size, @a num_data. It stores the result
 * of the bit operation between @a bit_mask and @a data. The bit operation is skipped
 * and the value in array, @a data, is adopted for the elements where corresponding
 * elements in @a edit_mask is false. The pointer of @a out is allowed to be equal to
 * that of @a in (@a result == @a data), indicating in-place operation.@n
 * must-be-aligned
 * @return Status code
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseImplicationUint8)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint8_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;
/**
 * @copybrief sakura_OperateBitwiseImplicationUint8
 * @copydetails sakura_OperateBitwiseImplicationUint8
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseImplicationUint32)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint32_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @brief Invoke bitwise NOT operation of an array.
 * @details Invokes the following bit operation to @a i- th element of @a result :
 * @code
 * result [i] = edit_mask[i] ? ~data[i] : data[i]
 * @endcode
 *
 * @note
 * No operation is done when the data array is zero length, i.e., @a num_data = 0.
 *
 * @param[in] num_data The number of elements in the arrays, @a data,
 * @a edit_mask, and @a result.
 * @param[in] data An input array of size, @a num_data. The bit operation
 * is invoked to this array.@n
 * must-be-aligned
 * @param[in] edit_mask A boolean mask array of size, @a num_data. The bit operation
 * is skipped for the elements with the value, false.@n
 * must-be-aligned
 * @param[out] result The output array of size, @a num_data. It stores the result
 * of the bit operation to @a data. The bit operation is skipped
 * and the value in array, @a data, is adopted for the elements where corresponding
 * elements in @a edit_mask is false. The pointer of @a out is allowed to be equal to
 * that of @a in (@a result == @a data), indicating in-place operation.@n
 * must-be-aligned
 * @return Status code
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseNotUint8)(
		size_t num_data, uint8_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint8_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;
/**
 * @copybrief sakura_OperateBitwiseNotUint8
 * @copydetails sakura_OperateBitwiseNotUint8
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseNotUint32)(
		size_t num_data, uint32_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint32_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @brief Invoke bit operation OR between a bit mask and an array.
 * @details Invokes the following bit operation to @a i- th element of @a result :
 * @code
 * result [i] = edit_mask[i] ? (data[i] | bit_mask) : data[i]
 * @endcode
 *
 * @note
 * No operation is done when the data array is zero length, i.e., @a num_data = 0.
 *
 * @note
 * The function can also be used to invoke converse implication of @a data and @a bit_mask .
 * Input the complement of @a bit_mask (~@a bit_mask ) to invoke converse implication.
 * For details of converse implication, see, e.g.,@n
 * http://en.wikipedia.org/wiki/Truth_table
 *
 * @param[in] bit_mask A bit mask. The bit operation is invoked
 * between this value and the array, @a data.
 * @param[in] num_data The number of elements in the arrays, @a data,
 * @a edit_mask, and @a result.
 * @param[in] data An input array of size, @a num_data. The bit operation
 * is invoked between this array and @a bit_mask.@n
 * must-be-aligned
 * @param[in] edit_mask A boolean mask array of size, @a num_data. The bit operation
 * is skipped for the elements with the value, false.@n
 * must-be-aligned
 * @param[out] result The output array of size, @a num_data. It stores the result
 * of the bit operation between @a bit_mask and @a data. The bit operation is skipped
 * and the value in array, @a data, is adopted for the elements where corresponding
 * elements in @a edit_mask is false. The pointer of @a out is allowed to be equal to
 * that of @a in (@a result == @a data), indicating in-place operation.@n
 * must-be-aligned
 * @return Status code
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseOrUint8)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint8_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;
/**
 * @copybrief sakura_OperateBitwiseOrUint8
 * @copydetails sakura_OperateBitwiseOrUint8
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseOrUint32)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint32_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @brief Invoke bit operation XOR between a bit mask and an array.
 * @details Invokes the following bit operation to @a i- th element of @a result :
 * @code
 * result [i] = edit_mask[i] ? (data[i] ^ bit_mask) : data[i]
 * @endcode
 *
 * @note
 * No operation is done when the data array is zero length, i.e., @a num_data = 0.
 *
 * @note
 * The function can also be used to invoke bitwise XNOR operation of @a data and @a bit_mask .
 * Input the complement of @a bit_mask (~@a bit_mask ) to invoke bitwise XNOR operation.
 * For details of bitwise XNOR operation, see, e.g.,@n
 * http://en.wikipedia.org/wiki/Truth_table
 *
 * @param[in] bit_mask A bit mask. The bit operation is invoked
 * between this value and the array, @a data.
 * @param[in] num_data The number of elements in the arrays, @a data,
 * @a edit_mask, and @a result.
 * @param[in] data An input array of size, @a num_data. The bit operation
 * is invoked between this array and @a bit_mask.@n
 * must-be-aligned
 * @param[in] edit_mask A boolean mask array of size, @a num_data. The bit operation
 * is skipped for the elements with the value, false.@n
 * must-be-aligned
 * @param[out] result The output array of size, @a num_data. It stores the result
 * of the bit operation between @a bit_mask and @a data. The bit operation is skipped
 * and the value in array, @a data, is adopted for the elements where corresponding
 * elements in @a edit_mask is false. The pointer of @a out is allowed to be equal to
 * that of @a in (@a result == @a data), indicating in-place operation.@n
 * must-be-aligned
 * @return Status code
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseXorUint8)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint8_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;
/**
 * @copybrief sakura_OperateBitwiseXorUint8
 * @copydetails sakura_OperateBitwiseXorUint8
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseXorUint32)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint32_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @brief Enumerations to define interpolation types.
 */
typedef enum {
	/**
	 * @brief Nearest interpolation
	 */LIBSAKURA_SYMBOL(InterpolationMethod_kNearest),

	/**
	 * @brief Linear interpolation
	 */LIBSAKURA_SYMBOL(InterpolationMethod_kLinear),

	/**
	 * @brief Polynomial interpolation
	 */LIBSAKURA_SYMBOL(InterpolationMethod_kPolynomial),

	/**
	 * @brief Spline interpolation (Natural cubic spline)
	 */LIBSAKURA_SYMBOL(InterpolationMethod_kSpline),

	/**
	 * @brief Number of interpolation methods implemented
	 */LIBSAKURA_SYMBOL(InterpolationMethod_kNumElements)
}LIBSAKURA_SYMBOL(InterpolationMethod);

/**
 * @brief Perform one-dimensional interpolation.
 * @details It performs one-dimensional interpolation based on two input arrays @a base_position
 * and @a base_data. Size of input array is @a num_base for @a base_position and @a num_interpolated
 * times @a num_array for @a base_data, where @a num_array is a number of arrays that are passed to the
 * function so that interpolation on multiple arrays can be performed simultaneously.
 * One can set boolean mask for @a base_data using @a base_mask, which has same array shape as
 * @a base_data. Mask value is true if data is valid while the value is false if the data is invalid.
 * Invalid data will be excluded from interpolation.
 *
 * List of locations where interpolated value is evaluated have to be specified via @a interpolate_position
 * whose length is @a num_interpolated. Interpolation result on each point specified by
 * @a interpolate_position is stored to @a interpolated_data. No extrapolation will be performed. Instead,
 * out of range points are filled by the value of nearest points.
 * Output boolean mask is @a interpolated_mask. Data will be invalid if corresponding mask is false, i.e.,
 * the value is just a nominal one but a result of actual interpolation. The mask will be false when
 * interpolation is skipped due to insufficient number of valid data elements.
 *
 * The function returns result status. In the successful run, returned value is
 * @link sakura_Status::sakura_Status_kOK sakura_Status_kOK @endlink while appropriate
 * error status will be returned for failure:
 * @link sakura_Status::sakura_Status_kInvalidArgument sakura_Status_kInvalidArgument @endlink
 * for invalid input arguments,
 * @link sakura_Status::sakura_Status_kNoMemory sakura_Status_kNoMemory @endlink
 * for memory allocation error for internal variables, and
 * @link sakura_Status::sakura_Status_kUnknownError sakura_Status_kUnknownError @endlink
 * for other unknown error.
 * Possible reason for
 * @link sakura_Status::sakura_Status_kInvalidArgument sakura_Status_kInvalidArgument @endlink
 * is either of (1) invalid @a interpolation_method or (2) arrays are not aligned.
 *
 * @pre @a base_position and @a interpolate_position must be sorted. Also, these arrays should not
 * have duplicate values.
 *
 * @par Difference between sakura_InterpolateXAxisFloat and sakura_InterpolateYAxisFloat:
 * Difference between these two similar functions is a memory layout of @a base_data and
 * @a interpolated_data. The former assumes the layout like [num_array][num_base] so that
 * @a base_data should store the data in the following order,
 * @verbatim data0[0], data0[1], ..., data1[0], data1[1], ... @endverbatim
 * On the other hand, the latter requires [num_base][num_array], i.e,
 * @verbatim data0[0], data1[0], ..., data0[1], data1[1], ... @endverbatim
 * Result array, @a interpolated_data, follows the memory layout required for @a base_data.
 *
 * @par Impact of sort order on performance:
 * When input arrays, @a base_position and/or @a base_data, are sorted in descending order,
 * the arrays are internally reversed in ascending order and then stored in working arrays.
 * Therefore, descending inputs may cause degradation of performance compared with ascending inputs.
 *
 * @par Note on polynomial interpolation:
 * Note that @a polynomial_order defines maximum order for polynomial interpolation.
 * In other words, it doesn't assure the interpolation to be specified order.
 * For example, suppose that @a polynomial_order is 2 and @a num_base is also 2.
 * In this case, effective polynomial order is 1 since we can obtain unique polynomial
 * with order 1, not 2, that passes through given two points.
 * Note also that @a polynomial_order 0 is equivalent to nearest interpolation.
 *
 * @param[in] interpolation_method Interpolation method.
 * @param[in] polynomial_order Maximum polynomial order for polynomial interpolation.
 * Actual order will be determined by a balance
 * between @a polynomial_order and @a num_base.
 * This parameter is effective only when @a interpolation_method is
 * @link sakura_InterpolationMethod::sakura_InterpolationMethod_kPolynomial sakura_InterpolationMethod_kPolynomial @endlink.
 * In other interpolation methods, it is ignored.
 * @param[in] num_base Number of elements for data points. Its value must be greater than 0.
 * @param[in] base_position Position of data points. Its length must be @a num_base.
 * It must be sorted either ascending or descending.
 * must-be-aligned
 * @param[in] num_array Number of arrays given in @a base_data.
 * @param[in] base_data Value of data points. Its length must be @a num_base times @a num_array.
 * must-be-aligned
 * @param[in] base_mask Boolean mask for data. Its length must be @a num_base times @a num_array.
 * False points will be excluded from the interpolation
 * must-be-aligned
 * @param[in] num_interpolated Number of elements for points that wants to get
 * interpolated value.
 * @param[in] interpolated_position Location of points that wants to get interpolated
 * value. Its length must be @a num_interpolated.
 * must-be-aligned
 * @param[out] interpolated_data Storage for interpolation result. Its length must be
 * @a num_interpolated times @a num_array.
 * must-be-aligned
 * @param[out] interpolated_mask Boolean mask for interpolation result. Its length must be
 * @a num_interpolated times @a num_array.
 * must-be-aligned
 * @return Status code.
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InterpolateXAxisFloat)(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_base,
		double const base_position[/*num_base*/], size_t num_array,
		float const base_data[/*num_base*num_array*/],
		bool const base_mask[/*num_base*num_array*/], size_t num_interpolated,
		double const interpolated_position[/*num_interpolated*/],
		float interpolated_data[/*num_interpolated*num_array*/],
		bool interpolated_mask[/*num_interpolated*num_array*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @copybrief sakura_InterpolateXAxisFloat
 * @copydetails sakura_InterpolateXAxisFloat
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InterpolateYAxisFloat)(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_base,
		double const base_position[/*num_base*/], size_t num_array,
		float const base_data[/*num_base*num_array*/],
		bool const base_mask[/*num_base*num_array*/], size_t num_interpolated,
		double const interpolated_position[/*num_interpolated*/],
		float interpolated_data[/*num_interpolated*num_array*/],
		bool interpolated_mask[/*num_interpolated*num_array*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @brief Normalize data against reference value with scaling factor.
 * @details
 * Normalize the data against reference value with scaling factor. The function normalizes the @a data
 * by the assumption that
 * the value in @a reference is normalized to @a scaling_factor. Specifically, it will calculate,
 * @verbatim result = scaling_factor * (data - reference) / reference @endverbatim
 *
 * @n
 * The function returns result status. For successful run, return value will be
 * @link sakura_Status::sakura_Status_kOK sakura_Status_kOK @endlink.
 * On the other hand, it will return
 * @link sakura_Status::sakura_Status_kInvalidArgument sakura_Status_kInvalidArgument @endlink
 * if any invalid values are passed to arguments.
 * @n
 * @n
 * The function allows in-place calculation, i.e., @a result array can be either @a data or
 * @a reference. In this case, @a data or @a reference will be overwritten.
 *
 * Only difference from @ref sakura_CalibrateDataWithConstScalingFloat is whether
 * @a scaling_factor is array or scalar. Note that, due to this difference,
 * order of the arguments is slightly different between the two.
 *
 * @param[in] num_data Number of data
 * @param[in] scaling_factor Scaling factor. This is an array and number of
 * elements must be @a num_data.
 * must-be-aligned
 * @param[in] data Data to be normalized. Number of elements must be @a num_data.
 * must-be-aligned
 * @param[in] reference Reference data. Number of elements must be @a num_data.
 * must-be-aligned
 * @param[out] result Resulting normalized data. One can give same array with either
 * @a data or @a reference for in-place calculation. Number of elements must be @a num_data.
 * must-be-aligned
 *
 * @return Status code.
 *
 * MT-safe
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CalibrateDataWithArrayScalingFloat)(
		size_t num_data, float const scaling_factor[/*num_data*/],
		float const data[/*num_data*/], float const reference[/*num_data*/],
		float result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

/**
 * @brief Normalize data against reference value with scaling factor.
 * @details
 * Normalize the data against reference value with scaling factor. The function normalizes the @a data
 * by the assumption that
 * the value in @a reference is normalized to @a scaling_factor. Specifically, it will calculate,
 * @verbatim result = scaling_factor * (data - reference) / reference @endverbatim
 *
 * @n
 * The function returns result status. For successful run, return value will be
 * @link sakura_Status::sakura_Status_kOK sakura_Status_kOK @endlink.
 * On the other hand, it will return
 * @link sakura_Status::sakura_Status_kInvalidArgument sakura_Status_kInvalidArgument @endlink
 * if any invalid values are passed to arguments.
 * @n
 * @n
 * The function allows in-place calculation, i.e., @a result array can be either @a data or
 * @a reference. In this case, @a data or @a reference will be overwritten.
 *
 * See @ref sakura_CalibrateDataWithArrayScalingFloat about the difference between those
 * two similar functions.
 *
 * @param[in] scaling_factor Scaling factor. This is a scalar.
 * @param[in] num_data Number of data
 * @param[in] data Data to be normalized. Number of elements must be @a num_data.
 * must-be-aligned
 * @param[in] reference Reference data. Number of elements must be @a num_data.
 * must-be-aligned
 * @param[out] result Resulting normalized data. One can give same array with either
 * @a data or @a reference for in-place calculation. Number of elements must be @a num_data.
 * must-be-aligned
 *
 * @return Status code.
 *
 * MT-safe
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CalibrateDataWithConstScalingFloat)(
		float scaling_factor, size_t num_data, float const data[/*num_data*/],
		float const reference[/*num_data*/], float result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @brief Context struct for convolution
 */
struct LIBSAKURA_SYMBOL(Convolve1DContextFloat);

/**
 * @brief Create Gaussian kernel.
 *
 * @details
 * Create 1 dimensional Gaussian kernel according to a peak location and FWHM given by @a peak_location
 * and @a kernel_width.
 *
 * The value of the kernel for i-th element is calculated by integrating Gaussian
 * from i-1/2 to i+1/2.
 * Since definite integral of Gaussian doesn't have analytic
 * form, this function is implemented using error function.
 * Resulting @a kernel is normalized so that its sum to be 1.0.
 * Actual formula for @a kernel is as follows.
 *
 * @verbatim
s = kernel_width / sqrt(log(16))
l = (i - 1/2 - peak_location) / s
r = (i + 1/2 - peak_location) / s
kernel[i] = 1/2 * {erf(r) - erf(l)}
@endverbatim
 *
 * where sqrt and log are mathematical functions of square root and logarithm to natural base.
 * The function erf is an error function.
 * The variable s corresponds to sqrt(2) times standard deviation of the Gaussian.
 * The above calculation is implemented using std::erf. See the link below
 * for details about std::erf:
 *
 * - http://www.cplusplus.com/reference/cmath/erf/
 * - http://en.cppreference.com/w/cpp/numeric/math/erf
 *
 * @param[in] peak_location the peak location of Gaussian
 * @param[in] kernel_width FWHM (Full Width of Half Maximum) of Gaussian.
 * @a kernel_width must be greater than 0.
 * @param[in] num_kernel The number of elements in the @a kernel
 * @param[out] kernel Output kernel array
 * @n must-be-aligned
 *
 * @return Status code.
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateGaussianKernelFloat)(
		float peak_location, float kernel_width, size_t num_kernel, float kernel[])
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @brief Create context for convolution using Fourier transformation.
 * @details
 * @param[in] num_kernel The number of elements in the @a kernel.
 * @a num_kernel must be positive.  0 < num_kernel <= INT_MAX
 * @param[in] kernel Convolution kernel. All elements in @a kernel must not be Inf nor NaN.
 * @n must-be-aligned
 * @param[out] context Context for convolution. The context can be shared between threads.
 * It has to be destroyed by @ref sakura_DestroyConvolve1DContextFloat after use by
 * @ref sakura_Convolve1DFFTFloat .
 * Note also that null pointer will be set to @a *context
 * in case this function fails.
 *
 * @return Status code.
 *
 * MT-unsafe
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateConvolve1DContextFFTFloat)(
		size_t num_kernel, float const kernel[/*num_kernel*/],
		struct LIBSAKURA_SYMBOL(Convolve1DContextFloat) **context)
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @brief Convolution is performed.
 * @details Convolution operation is performed by shifting a kernel over the input data.
 * The operation is carried out only for elements of data whose mask values are true.
 * @param[in] num_kernel  The number of elements in the @a kernel.
 * @a num_kernel must be positive.  0 < @a num_kernel <= INT_MAX
 * @param[in] kernel Convolution kernel. All elements in @a kernel must not be Inf nor NaN.
 * @n must-be-aligned
 * @param[in] num_data
 * The number of elements in @a input_data, @a input_mask, @a output_data,
 * and @a output_weight. 0 < @a num_data <= INT_MAX
 * @param[in] input_data Input data. If corresponding element in @a input_mask is true,
 * the element in @a input_data must not be Inf nor NaN.
 * @n must-be-aligned
 * @param[in] input_mask Mask of input data. @a input_data is used in operation
 * if corresponding element of @a input_mask is true. If not, the corresponding
 * elements in @a input_data is ignored.
 * @n must-be-aligned
 * @param[out] output_data Output data.
 * @n must-be-aligned
 * @param[out] output_weight Total weight of kernel summed up to corresponding elements
 * of @a output_data in the convolution operation.
 * @n must-be-aligned
 * @return Status code.
 *
 * @note
 * The function does not define mask of convolved data, @a output_data .
 * User is responsible for defining it if necessary.
 * Instead, the function gives @a output_weight to support user define
 * mask after convolution. @n
 * In general, weight value is expected to be smaller if some of elements,
 * that are supposed to contribute in convolution, are ignored due to
 * @a input_mask being false. Therefore values of @a output_weight can be
 * used as an indicator of degree of mask affecting the corresponding
 * elements in @a output_data. @n
 * Note, however, elements near the edge of array also have small weights
 * because of boundary effect.
 *
 * MT-safe
 *
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Convolve1DFloat)(size_t num_kernel,
		float const kernel[/*num_kernel*/], size_t num_data,
		float const input_data[/*num_data*/],
		bool const input_mask[/*num_data*/], float output_data[/*num_data*/],
		float output_weight[/*num_data*/]) LIBSAKURA_NOEXCEPT;
/**
 * @brief Convolution is performed using Fourier transformation.
 * @details Convolution operation is performed using Fourier transformation of data and kernel arrays.
 * The kernel is stored in a context in an internal format.
 * The operation is carried out for all elements of data.
 * @param[in] context
 * The context created by @ref sakura_CreateConvolve1DContextFFTFloat.
 * @param[in] num_data
 * The number of elements in @a input_data and @a output_data. @a num_data must be equal to @a num_kernel in @a context .
 * 0 < @a num_data <= INT_MAX
 * @param[in] input_data Input data. All elements in @a input_data must not be Inf nor NaN.
 * @n must-be-aligned
 * @param[out] output_data Output data. The pointer of @a out is allowed to be equal to
 *  that of @a in (@a input_data == @a output_data), indicating in-place operation.
 * @n must-be-aligned
 * @return Status code.
 *
 * MT-safe
 *
 * (But see @ref sakura_CreateConvolve1DContextFFTFloat for detail about the thread-safety)
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Convolve1DFFTFloat)(
		struct LIBSAKURA_SYMBOL(Convolve1DContextFloat) const *context,
		size_t num_data, float const input_data[/*num_data*/],
		float output_data[/*num_data*/]) LIBSAKURA_NOEXCEPT;
/**
 * @brief Destroy context for convolution.
 * @details
 * @param[in] context
 * The context created by @ref sakura_CreateConvolve1DContextFFTFloat
 * and to be destroyed.
 * @return Status code.
 *
 * MT-unsafe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DestroyConvolve1DContextFloat)(
		struct LIBSAKURA_SYMBOL(Convolve1DContextFloat) *context)
				LIBSAKURA_NOEXCEPT;

/**
 * @brief Compute coefficients of simultaneous equations used for Least-Square fitting.
 * @details
 * Suppose fitting ( @a num_data ) discrete data points yi with a linear
 * combination of ( @a num_model_bases ) bases (ai, bi, ..., ni), which are
 * also given as ( @a num_data ) discrete points. Assuming the best-fit model
 * is given as (A * ai + B * bi + ... + N * ni), where (A, B, C, ...) are
 * the coefficients to be solved, these values are connected via the following
 * simultaneous equations known as normal equation:
 *
 * @image html GetCoefficientsForLeastSquareFitting.png
 *
 * Note that the summation means all the data points except masked ones
 * are to be added.
 * This function computes the coefficients of the above simultaneous equations, i.e.,
 * the elements of the matrix at the left side and of the vector at the right side.
 * @par
 * @param[in] num_data The number of elements in the arrays @a data
 * and @a mask, and also the number of elements in each model data
 * (i.e., discrete values of basis function) consisting the entire model.
 * It must be a positive number.
 * @param[in] data Input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] mask The input mask data with length of @a num_data .
 * The @a i th element of @a data is included in input data if the
 * @a i th element of @a mask is true, while it is excluded from input
 * data if the corresponding element of @a mask is false.
 * @n must-be-aligned
 * @param[in] num_model_bases Number of model basis functions. It must be
 * in range 0 < @a num_model_bases <= @a num_data .
 * @param[in] basis_data A 1D array containing values of all basis functions
 * concatenated. Loop for basis index must be inside of that for data index,
 * i.e., the @a n -th data of the @a m -th model should be stored at
 * @a basis_data [ @a num_data * @a (n-1) + @a (m-1) ]. Its length must be
 * equal to ( @a num_model_bases * @a num_data ).
 * @n must-be-aligned
 * @param[in] num_lsq_bases The number of model basis functions to be used
 * in actual fitting. It must be a positive number and must not exceed
 * @a num_model_bases .
 * @param[in] use_bases_idx A 1D array containing indices of basis model
 * that are to be used for fitting. As for fitting types other than
 * sinusoidal, it should be always [0, 1, 2, ..., (num_lsq_bases-1)].
 * Element values must not be duplicated, and must be in ascending order.
 * @n must-be-aligned
 * @param[out] lsq_matrix A 1D array containing the values of a matrix
 * at the left side of simultaneous equations for least-square fitting.
 * Its length should therefore be equal to ( @a num_lsq_bases * @a num_lsq_bases ).
 * Loop for columns comes inside that for rows, i.e., the value at the
 * @a m -th row and @a n -th column is stored at @a out [ @a
 * num_lsq_bases * ( @a m -1) + ( @a n -1)], though @a out is actually
 * symmetric.
 * @n must-be-aligned
 * @param[out] lsq_vector The values of a vector at the right side of
 * simultaneous equations for least-square fitting. Its length should be
 * equal to @a num_model_bases .
 * @n must-be-aligned
 * @return Status code.
 * @link sakura_Status::sakura_Status_kOK sakura_Status_kOK @endlink
 * if finished successfully,
 * @link sakura_Status::sakura_Status_kInvalidArgument sakura_Status_kInvalidArgument @endlink
 * in case parameters does not meet the above criteria,
 * @link sakura_Status::sakura_Status_kNG sakura_Status_kNG @endlink
 * in case the number of unmasked data (for which mask is @a false ) is
 * less than the number of simultaneous equations, and
 * @link sakura_Status::sakura_Status_kUnknownError sakura_Status_kUnknownError @endlink
 * in case other exceptions emitted internally.
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetLSQCoefficientsDouble)(
		size_t const num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/], size_t const num_model_bases,
		double const basis_data[/*num_model_bases*num_data*/],
		size_t const num_lsq_bases,
		size_t const use_bases_idx[/*num_lsq_bases*/],
		double lsq_matrix[/*num_lsq_bases*num_lsq_bases*/],
		double lsq_vector[/*num_lsq_bases*/])
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @copybrief sakura_GetLSQCoefficientsDouble
 * @details
 * This function updates the coefficients of normal equation for LSQ fitting
 * created by @ref sakura_GetLSQCoefficientsDouble , by subtracting values
 * corresponding to data points which have been used in the previous
 * calculation but not this time. This is faster than newly calculating
 * coefficients if the number of points to be excluded this time is less
 * than half of those previously used.
 * @par
 * @param[in] num_data The number of elements in the array @a data and the
 * number of elements in each model data (i.e., discrete values of basis
 * function) consisting the entire model.
 * It must be a positive number.
 * @param[in] data Input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] mask The input mask data with length of @a num_data .
 * The @a i th element of @a data is included in input data if the
 * @a i th element of @a mask is true, while it is excluded from input
 * data if the corresponding element of @a mask is false.
 * @n must-be-aligned
 * @param[in] num_exclude_indices The number of data points to be excluded
 * this time. The range of allowed value is between 0 and @a num_data .
 * @param[in] exclude_indices An array containing indices of data points
 * (the row index of @a basis_data ) to be excluded this time. The indices
 * must be stored in the first @a num_exclude_indices elements. Its length
 * should be @a num_exclude_indices .
 * @n must-be-aligned
 * @param[in] num_model_bases Number of model basis functions. It must be
 * in range 0 < @a num_model_bases <= @a num_data .
 * @param[in] basis_data A 1D array containing values of all basis functions
 * concatenated. Loop for basis index must be inside of that for data index,
 * i.e., the @a n -th data of the @a m -th model should be stored at
 * @a basis_data [ @a num_data * @a (n-1) + @a (m-1) ]. Its length must be
 * equal to ( @a num_model_bases * @a num_data ).
 * @n must-be-aligned
 * @param[in] num_lsq_bases The number of model basis functions to be used
 * in actual fitting. It must be in range 0 < @a num_lsq_bases <=
 * @a num_model_bases.
 * @param[in] use_bases_idx A 1D array containing indices of basis model
 * that are to be used for fitting. As for fitting types other than
 * sinusoidal, it should be always [0, 1, 2, ..., (num_lsq_bases-1)].
 * Element values must be in ascending order.
 * @n must-be-aligned
 * @param[in,out] lsq_matrix A 1D array containing the values of a matrix
 * at the left side of simultaneous equations for least-square fitting.
 * Its length should therefore be equal to ( @a num_lsq_bases * @a num_lsq_bases ).
 * Loop for columns comes inside that for rows, i.e., the value at the
 * @a m -th row and @a n -th column is stored at @a out [ @a
 * num_lsq_bases * ( @a m -1) + ( @a n -1)], though @a out is actually
 * symmetric.
 * @param[in,out] lsq_vector The values of a vector at the right side of
 * simultaneous equations for least-square fitting. Its length should be
 * equal to @a num_lsq_bases .
 * @n must-be-aligned
 * @return Status code.
 * @link sakura_Status::sakura_Status_kOK sakura_Status_kOK @endlink
 * if finished successfully,
 * @link sakura_Status::sakura_Status_kInvalidArgument sakura_Status_kInvalidArgument @endlink
 * in case parameters does not meet the above criteria, and
 * @link sakura_Status::sakura_Status_kUnknownError sakura_Status_kUnknownError @endlink
 * in case other exceptions emitted internally.
 * @par Caution:
 * Users must be careful in using this function about which and how many
 * data are to be excluded not to fall into destructive cases that the
 * number of used data becomes less than @a num_model_bases or not to
 * exclude the same data in duplicate.
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UpdateLSQCoefficientsDouble)(
		size_t const num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/], size_t const num_exclude_indices,
		size_t const exclude_indices[/*num_data*/],
		size_t const num_model_bases,
		double const basis_data[/*num_model_bases*num_data*/],
		size_t const num_lsq_bases,
		size_t const use_bases_idx[/*num_lsq_bases*/],
		double lsq_matrix[/*num_lsq_bases*num_lsq_bases*/],
		double lsq_vector[/*num_lsq_bases*/])
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @brief Solve simultaneous equations via LU decomposition.
 * @details
 * Suppose solving simultaneous equations A x = y to derive x, where A
 * is a square matrix of @a num_equations rows and columns, and x and y
 * are vectors with length of @a num_equations . Given A and y values,
 * this function computes x values using LU decomposition of A.
 * @par
 * @param[in] num_equations Number of equations.
 * @param[in] in_matrix A 1D array containing values of the matrix A in
 * the left side of the above simultaneous equations. Loop for columns
 * comes inside that for rows, i.e., the value at the @a m -th row and
 * @a n -th column is stored at
 * @a in_matrix [ @a num_equations * ( @a m -1) + ( @a n -1)].
 * Its length must be (@a num_equations * @a num_equations).
 * @n must-be-aligned
 * @param[in] in_vector A 1D array containing values of the vector y in
 * the right side of the above simultaneous equations. Its length must be
 * @a num_equations .
 * @n must-be-aligned
 * @param[out] out The solution (x in the above equations). Its length
 * must be @a num_equations .
 * @n must-be-aligned
 * @return Status code.
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLUDouble)(
		size_t num_equations,
		double const in_matrix[/*num_equations*num_equations*/],
		double const in_vector[/*num_equations*/],
		double out[/*num_equations*/])
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

//LM part----------------------------------------------------------------
/**
 * @brief Fit Gaussian to the input data using Levenberg-Marquardt method.
 * @details
 * Fit Gaussian to the input data using Levenberg-Marquardt method.
 * The input data can contain multiple Gaussian profiles, however,
 * note that fitting multiple Gaussians at once is quite difficult
 * and thus the results may not be correct.
 * @param[in] num_data Number of input data. It must be equal to or
 * larger than the number of Gaussian parameters (3 * @a num_peaks ).
 * @param[in] data The input data with length of @a num_data .
 * @param[in] mask The input mask data with length of @a num_data .
 * @param[in] num_peaks Number of Gaussians to be fitted.
 * @param[in,out] height An array to store initial guess of Gaussian
 * heights. Its length must be @a num_peaks and will be overwritten
 * with the fitting results.
 * @param[in,out] center An array to store initial guess of Gaussian
 * centers. Its length must be @a num_peaks and will be overwritten
 * with the fitting results.
 * @param[in,out] sigma An array to store initial guess of Gaussian
 * sigmas. Its length must be @a num_peaks and will be overwritten
 * with the fitting results.
 * @param[out] err_height An array to store errors of fitted
 * Gaussian heights. Its length must be @a num_peaks .
 * @param[out] err_center An array to store errors of fitted
 * Gaussian centers. Its length must be @a num_peaks .
 * @param[out] err_sigma An array to store errors of fitted
 * Gaussian sigmas. Its length must be @a num_peaks .
 * @return Status code.
 *
 * MT-safe
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(LMFitGaussianFloat)(
		size_t const num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/], size_t const num_peaks,
		double height[/*num_peaks*/], double err_height[/*num_peaks*/],
		double center[/*num_peaks*/], double err_center[/*num_peaks*/],
		double sigma[/*num_peaks*/], double err_sigma[/*num_peaks*/])
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @brief Enumerations to define least-square fitting specific error code.
 */
typedef enum {
	/**
	 * @brief OK
	 */LIBSAKURA_SYMBOL(LSQFitStatus_kOK) = 0,

	/**
	 * @brief NG
	 */LIBSAKURA_SYMBOL(LSQFitStatus_kNG) = 1,

	/**
	 * @brief Not enough data for least-square fitting
	 */LIBSAKURA_SYMBOL(LSQFitStatus_kNotEnoughData) = 2,

	/**
	 * @brief Number of error codes implemented
	 */LIBSAKURA_SYMBOL(LSQFitStatus_kNumElements)
}LIBSAKURA_SYMBOL(LSQFitStatus);

/**
 * @brief Enumerations to define type of least-square fitting models.
 */
typedef enum {
	/**
	 * @brief Polynomial
	 */LIBSAKURA_SYMBOL(LSQFitType_kPolynomial),

	/**
	 * @brief Chebyshev Polynomial
	 */LIBSAKURA_SYMBOL(LSQFitType_kChebyshev),

	/**
	 * @brief Number of fitting models
	 */LIBSAKURA_SYMBOL(LSQFitType_kNumElements)
}LIBSAKURA_SYMBOL(LSQFitType);

/**
 * @brief Context struct for least-square fitting
 */
struct LIBSAKURA_SYMBOL(LSQFitContextFloat);

/**
 * @brief Create an object containing model data for least-square fitting.
 * @details
 * @note A lsqfit context object can not be shared between
 * threads, as it contains working areas exclusive for a specific
 * thread.
 * @param[in] lsqfit_type Type of basis function. It should
 * be either of
 * @link sakura_LSQFitType::sakura_LSQFitType_kPolynomial sakura_LSQFitType_kPolynomial @endlink
 * or
 * @link sakura_LSQFitType::sakura_LSQFitType_kChebyshev sakura_LSQFitType_kChebyshev @endlink .
 * @param[in] order Polynomial order. It must be positive or zero.
 * @param[in] num_data Number of data to fit. It must be equal
 * to or larger than the number of model bases, which is
 * ( @a order+1 ), thus the smallest allowed value of
 * @a num_data is 1.
 * @param[out] context Pointer to pointer to an object containing
 * model data. When no longer used, the object must be destroyed by
 * @ref sakura_DestroyLSQFitContextFloat .
 * Note also that null pointer will be set to @a *context
 * in case this function fails.
 * @return Status code.
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateLSQFitContextPolynomialFloat)(
LIBSAKURA_SYMBOL(LSQFitType) const lsqfit_type, uint16_t order, size_t num_data,
		struct LIBSAKURA_SYMBOL(LSQFitContextFloat) **context)
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @copybrief sakura_CreateLSQFitContextPolynomialFloat
 * @details
 * @note A lsqfit context object can not be shared between
 * threads, as it contains working areas exclusive for a specific
 * thread.
 * @param[in] npiece Number of spline pieces. It must be a
 * positive value.
 * @param[in] num_data Number of data to fit. It must be equal
 * to or larger than the number of model bases ( @a npiece+3 ),
 * thus the smallest allowed value of @a num_data is 4.
 * @param[out] context Pointer to pointer to an object containing
 * model data. When no longer used, the object must be destroyed by
 * @ref sakura_DestroyLSQFitContextFloat .
 * Note also that null pointer will be set to @a *context
 * in case this function fails.
 * @return Status code.
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateLSQFitContextCubicSplineFloat)(
		uint16_t npiece, size_t num_data,
		struct LIBSAKURA_SYMBOL(LSQFitContextFloat) **context)
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @copybrief sakura_CreateLSQFitContextPolynomialFloat
 * @details
 * @note A lsqfit context object can not be shared between
 * threads, as it contains working areas exclusive for a specific
 * thread.
 * @param[in] nwave Maximum wave number of sinusoids. It must be
 * positive or zero.
 * @param[in] num_data Number of data to fit. It must be equal
 * to or larger than the number of model bases plus one,
 * which is ( @a nwave*2+2 ), thus the smallest allowed value of
 * @a num_data is 2.
 * @param[out] context Pointer to pointer to an object containing
 * model data. When no longer used, the object must be destroyed by
 * @ref sakura_DestroyLSQFitContextFloat .
 * Note also that null pointer will be set to @a *context
 * in case this function fails.
 * @return Status code.
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateLSQFitContextSinusoidFloat)(
		uint16_t nwave, size_t num_data,
		struct LIBSAKURA_SYMBOL(LSQFitContextFloat) **context)
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @brief Destroy an object containing model data for least-square fitting.
 * @details
 * @param[in] context A context created by
 * @ref sakura_CreateLSQFitContextPolynomialFloat or
 * @ref sakura_CreateLSQFitContextCubicSplineFloat or
 * @ref sakura_CreateLSQFitContextSinusoidFloat .
 * @return Status code.
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DestroyLSQFitContextFloat)(
		struct LIBSAKURA_SYMBOL(LSQFitContextFloat) *context)
				LIBSAKURA_NOEXCEPT;

/**
 * @brief Fit a polynomial curve to input data.
 * @details A polynomial or Chebyshev polynomial is fitted to input data
 * based on Least-Square method. As output, coefficients of bases of the
 * best-fit curve, the best-fit curve value itself and the residuals
 * (i.e., input - best-fit curve) are stored in @a coeff , @a best_fit
 * and @a residual , respectively. If @a num_fitting_max greater than 1
 * is given, fitting is executed recursively. Once best-fit curve is
 * subtracted from input data, mean and standard deviation of the
 * residual is computed to update @a mask so that @a mask has false
 * value at data points where residual value exceeds threshold defined
 * by @a clip_threshold_sigma , then Least-Square fitting is again
 * performed to the input data with updated @a mask to update @a coeff
 * values. This procedure is repeatedly done for ( @a num_fitting_max-1 )
 * times or until the fitting result converges. Once fitting is done,
 * the updated mask information is stored in @a final_mask .
 * @param[in] context A context created by @ref sakura_CreateLSQFitContextPolynomialFloat .
 * @param[in] order Polynomial order. It should not exceed the @a order
 * specified in creation of @a context .
 * @param[in] num_data Number of elements in the arrays @a data, @a mask,
 * @a final_mask, and @a out. It must be equal to @a num_data which was
 * given to @ref sakura_CreateLSQFitContextPolynomialFloat to create
 * @a context .
 * @param[in] data Input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] mask Input mask data with length of @a num_data .
 * The @a i th element of @a data is included in input spectrum if the
 * @a i th element of @a mask is true, while it is excluded from input
 * spectrum if the corresponding element of @a mask is false.
 * @n must-be-aligned
 * @param[in] clip_threshold_sigma Threshold of clipping in unit of
 * sigma. It must be a positive value.
 * @param[in] num_fitting_max Upper limit of how many times fitting is
 * performed recursively. Before executing the second or later fitting,
 * outlier in @a data is masked via clipping to be not used. If 1 is
 * given, fitting is done just once and no clipping will be applied.
 * If 0 is given, fitting is not executed and the values of @a residual
 * should be identical with those of @a data .
 * @param[in] num_coeff Number of elements in the array @a coeff .
 * In case @a coeff is not null pointer, it must be equal to ( @a order
 * +1), while the value is not checked when @a coeff is null pointer.
 * @param[out] coeff Coefficients of the polynomial bases. Its length
 * must be @a num_coeff . Coefficient values are stored in ascending
 * order of polynomial order: coefficient of constant term comes first,
 * followed by those of first order, second order, and so on. Null
 * pointer can be given in case users do not need this value.
 * @n must-be-aligned
 * @param[out] best_fit The best-fit curve data, i.e., the result of
 * least-square fitting itself. Its length must be @a num_data .
 * @a data can be set to @a best_fit if users want to overwrite it
 * in-place. Null pointer can be given in case users do not need
 * this value.
 * @n must-be-aligned
 * @param[out] residual Residual (input - best-fit) data. Its length
 * must be @a num_data . @a data can be set to @a residual if users
 * want to overwrite it in-place. Null pointer can be given in case
 * users do not need this value.
 * @n must-be-aligned
 * @param[out] final_mask The final status of mask data after
 * recursive clipping procedure finish. Its length must be
 * @a num_data . @a mask can be set to @a final_mask if users want
 * to overwrite it in-place.
 * @n must-be-aligned
 * @param[out] rms The root-mean-square of @a residual .
 * @param[out] lsqfit_status LSQFit-specific error code.
 * @return Status code.
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(
		struct LIBSAKURA_SYMBOL(LSQFitContextFloat) const *context,
		uint16_t order, size_t num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/], float clip_threshold_sigma,
		uint16_t num_fitting_max, size_t num_coeff, double coeff[/*num_coeff*/],
		float best_fit[/*num_data*/], float residual[/*num_data*/],
		bool final_mask[/*num_data*/], float *rms,
		LIBSAKURA_SYMBOL(LSQFitStatus) *lsqfit_status)
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @brief Fit a cubic spline curve to input data.
 * @details A cubic spline curve is fitted to input data based on
 * Least-Square method. As output, coefficients of bases of the
 * best-fit curve, the best-fit curve value itself and the residuals
 * (i.e., input - best-fit curve) are stored in @a coeff ,
 * @a best_fit and @a residual , respectively.
 * If @a num_fitting_max greater than 1 is given, fitting is executed
 * recursively. Once best-fit curve is subtracted from input data, mean
 * and standard deviation of the residual is computed to update @a mask
 * so that @a mask has false value at data points where residual value
 * exceeds threshold defined by @a clip_threshold_sigma , then
 * Least-Square fitting is again performed to the input data with updated
 * @a mask to update @a coeff values. This procedure is repeatedly done
 * for ( @a num_fitting_max-1 ) times or until the fitting result converges.
 * Once fitting is done, the updated mask information is stored in
 * @a final_mask .
 * @param[in] context A context created by @ref sakura_CreateLSQFitContextCubicSplineFloat .
 * @param[in] num_pieces Number of spline pieces. It must be positive
 * and also must not exceed the number of spline pieces specified in
 * creation of @a context .
 * @param[in] num_data Number of elements in the arrays @a data, @a mask,
 * @a final_mask, and @a out. It must be equal to @a num_data which was
 * given to @ref sakura_CreateLSQFitContextCubicSplineFloat to create
 * @a context .
 * @param[in] data Input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] mask Input mask data with length of @a num_data .
 * The @a i th element of @a data is included in input spectrum if the
 * @a i th element of @a mask is true, while it is excluded from input
 * spectrum if the corresponding element of @a mask is false.
 * @n must-be-aligned
 * @param[in] clip_threshold_sigma Threshold of clipping in unit of
 * sigma. It must be a positive value.
 * @param[in] num_fitting_max Upper limit of how many times fitting is
 * performed recursively. Before executing the second or later fitting,
 * outlier in @a data is masked via clipping to be not used. If 1 is
 * given, fitting is done just once and no clipping will be applied.
 * If 0 is given, fitting is not executed and the values of @a residual
 * should be identical with those of @a data .
 * @param[out] coeff The coefficients of the best-fit cubic spline curve.
 * It must be a 2D array with type of double[ @a num_pieces ][4].
 * @a coeff[i] is for the @a i th spline piece from the left side.
 * The first element in each @a coeff[i] is of constant term, followed
 * by those of first, second, and third orders. Null pointer can be
 * given in case users do not need this value.
 * @n must-be-aligned
 * @param[out] best_fit The best-fit curve data, i.e., the result of
 * least-square fitting itself. Its length must be @a num_data .
 * @a data can be set to @a best_fit if users want to overwrite it
 * in-place. Null pointer can be given in case users do not need
 * this value.
 * @n must-be-aligned
 * @param[out] residual Residual (input - best-fit) data. Its length
 * must be @a num_data . @a data can be set to @a residual if users
 * want to overwrite it in-place. Null pointer can be given in case
 * users do not need this value.
 * @n must-be-aligned
 * @param[out] final_mask The final status of mask data after
 * recursive clipping procedure finish. Its length must be
 * @a num_data . @a mask can be set to @a final_mask if users want
 * to overwrite it in-place.
 * @n must-be-aligned
 * @param[out] rms The root-mean-square of @a residual .
 * @param[out] boundary A 1D array containing the boundary positions
 * of spline pieces, which are indices of @a data .
 * Its length must be ( @a num_pieces +1). The element values will be
 * stored in ascending order. The first element will always be zero,
 * the left edge of the first (left-most) spline piece, while the last
 * element will be @a num_data , which is the next of the right edge
 * of the last (right-most) spline piece.
 * @n must-be-aligned
 * @param[out] lsqfit_status LSQFit-specific error code.
 * @return Status code.
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(LSQFitCubicSplineFloat)(
		struct LIBSAKURA_SYMBOL(LSQFitContextFloat) const *context,
		size_t num_pieces, size_t num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/], float clip_threshold_sigma,
		uint16_t num_fitting_max, double coeff[/*num_pieces*/][4],
		float best_fit[/*num_data*/], float residual[/*num_data*/],
		bool final_mask[/*num_data*/], float *rms,
		size_t boundary[/*num_pieces+1*/],
		LIBSAKURA_SYMBOL(LSQFitStatus) *lsqfit_status)
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @brief Fit a sinusoidal curve to input data.
 * @details A sinusoidal curve is fitted to input data based on
 * Least-Square method. As output, coefficients of bases of the
 * best-fit curve, the best-fit curve value itself and the residuals
 * (i.e., input - best-fit curve) are stored in @a coeff , @a best_fit
 * and @a residual , respectively. If @a num_fitting_max greater than 1
 * is given, fitting is executed recursively. Once best-fit curve is
 * subtracted from input data, mean and standard deviation of the
 * residual is computed to update @a mask so that @a mask has false
 * value at data points where residual value exceeds threshold defined
 * by @a clip_threshold_sigma , then Least-Square fitting is again
 * performed to the input data with updated @a mask to update @a coeff
 * values. This procedure is repeatedly done for ( @a num_fitting_max-1 )
 * times or until the fitting result converges. Once fitting is done,
 * the updated mask information is stored in @a final_mask .
 * @param[in] context A context created by
 * @ref sakura_CreateLSQFitContextSinusoidFloat .
 * @param[in] num_nwave The number of elements in the array @a nwave .
 * @param[in] nwave Wave numbers within the index range of @a data
 * to be used for sinusoidal fitting. The values must be positive or
 * zero (for constant term), but not exceed the @a maximum wave number
 * specified in creation of @a context . The values must be stored in
 * ascending order and must not be duplicate.
 * @param[in] num_data Number of elements in the arrays @a data, @a mask,
 * @a final_mask, and @a out. It must be equal to @a num_data which was
 * given to @ref sakura_CreateLSQFitContextSinusoidFloat to create
 * @a context .
 * @param[in] data Input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] mask Input mask data with length of @a num_data .
 * The @a i th element of @a data is included in input spectrum if the
 * @a i th element of @a mask is true, while it is excluded from input
 * spectrum if the corresponding element of @a mask is false.
 * @n must-be-aligned
 * @param[in] clip_threshold_sigma Threshold of clipping in unit of
 * sigma. It must be a positive value.
 * @param[in] num_fitting_max Upper limit of how many times fitting is
 * performed recursively. Before executing the second or later fitting,
 * outlier in @a data is masked via clipping to be not used. If 1 is
 * given, fitting is done just once and no clipping will be applied.
 * If 0 is given, fitting is not executed and the values of @a residual
 * should be identical with those of @a data .
 * @param[in] num_coeff The number of elements in the array @a coeff.
 * If @a coeff is not null pointer, it must be ( @a num_nwave*2-1 )
 * or ( @a num_nwave*2 ) in cases @a nwave contains zero or not,
 * respectively, and must not exceed @a num_data, while the value is
 * not checked when @a coeff is null pointer.
 * @param[out] coeff Coefficients of the sinusoidal fit. Its length
 * must be @a num_coeff . If @a nwave contains zero, the first element
 * is of the constant term. If not, the first and the second elements
 * are for sine and cosine terms for the smallest wave number in
 * @a nwave, respectively. Coefficients of sine and cosine terms for
 * the second-smallest wave number come next, and so on. Null pointer
 * can be given in case users do not need this value.
 * @n must-be-aligned
 * @param[out] best_fit The best-fit curve data, i.e., the result of
 * least-square fitting itself. Its length must be @a num_data .
 * @a data can be set to @a best_fit if users want to overwrite it
 * in-place. Null pointer can be given in case users do not need
 * this value.
 * @n must-be-aligned
 * @param[out] residual Residual (input - best-fit) data. Its length
 * must be @a num_data . @a data can be set to @a residual if users
 * want to overwrite it in-place. Null pointer can be given in case
 * users do not need this value.
 * @n must-be-aligned
 * @param[out] final_mask The final status of mask data after
 * recursive clipping procedure finish. Its length must be
 * @a num_data . @a mask can be set to @a final_mask if users want
 * to overwrite it in-place.
 * @n must-be-aligned
 * @param[out] rms The root-mean-square of @a residual .
 * @param[out] lsqfit_status LSQFit-specific error code.
 * @return Status code.
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(
		struct LIBSAKURA_SYMBOL(LSQFitContextFloat) const *context,
		size_t num_nwave, size_t const nwave[/*num_nwave*/], size_t num_data,
		float const data[/*num_data*/], bool const mask[/*num_data*/],
		float clip_threshold_sigma, uint16_t num_fitting_max, size_t num_coeff,
		double coeff[/*num_coeff*/], float best_fit[/*num_data*/],
		float residual[/*num_data*/], bool final_mask[/*num_data*/], float *rms,
		LIBSAKURA_SYMBOL(LSQFitStatus) *lsqfit_status)
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @brief Subtract best-fit polynomial model from input data.
 * @details
 * @param[in] context A context created by
 * @ref sakura_CreateLSQFitContextPolynomialFloat .
 * @param[in] num_data The number of elements in @a data and @a out.
 * It must be equal to @a num_data which was given to
 * @ref sakura_CreateLSQFitContextPolynomialFloat to create @a context .
 * @param[in] data The input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] num_coeff The number of elements in @a coeff .
 * It must be in range (0 < @a num_coeff <= @a order+1 ), where
 * @a order is the polynomial or Chebyshev polynomial order
 * specified in creation of @a context .
 * @param[in] coeff Coefficients of model data. Its length must
 * be @a num_coeff. Its first element must be of constant term,
 * followed by those of first order, second order,
 * and so on.
 * @n must-be-aligned
 * @param[out] out The output data. Its length must be @a num_data .
 * @n must-be-aligned
 * @return Status code.
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractPolynomialFloat)(
		struct LIBSAKURA_SYMBOL(LSQFitContextFloat) const *context,
		size_t num_data, float const data[/*num_data*/], size_t num_coeff,
		double const coeff[/*num_coeff*/], float out[/*num_data*/])
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @brief Subtract best-fit cubic spline model from input data.
 * @details
 * @param[in] context A context created by
 * @ref sakura_CreateLSQFitContextCubicSplineFloat .
 * @param[in] num_data The number of elements in @a data and @a out .
 * It must be equal to @a num_data which was given to
 * @ref sakura_CreateLSQFitContextCubicSplineFloat to create @a context .
 * @param[in] data The input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] num_pieces The number of spline pieces. If zero is
 * given, no subtraction executed.
 * @param[in] coeff Coefficients of cubic spline curve. It must be
 * a 2D array with type of double[ @a num_pieces ][4]. @a coeff[i]
 * is for the @a i th spline piece from the left side. The first
 * element in each @a coeff[i] must be of constant term, followed
 * by the ones of first, second, and third orders.
 * @n must-be-aligned
 * @param[in] boundary A 1D array containing the boundary positions,
 * which are indices of parameters @a data and @a out , of spline
 * pieces. Its length must be ( @a num_pieces +1). The element
 * values should be stored in ascending order. The first element
 * must always be zero, the left edge of the first (left-most) spline
 * piece, while the last element must be @a num_data , which is the
 * next of the right edge of the last (right-most) spline piece.
 * @n must-be-aligned
 * @param[out] out The output data. Its length must be @a num_data .
 * @n must-be-aligned
 * @return Status code.
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractCubicSplineFloat)(
		struct LIBSAKURA_SYMBOL(LSQFitContextFloat) const *context,
		size_t num_data, float const data[/*num_data*/], size_t num_pieces,
		double const coeff[/*num_pieces*/][4],
		size_t const boundary[/*num_pieces+1*/], float out[/*num_data*/])
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @brief Subtract best-fit sinusoidal model from input data.
 * @details
 * @param[in] context A context created by
 * @ref sakura_CreateLSQFitContextSinusoidFloat .
 * @param[in] num_data The number of elements in @a data and @a out.
 * It must be equal to @a num_data which was given to
 * @ref sakura_CreateLSQFitContextSinusoidFloat to create @a context .
 * @param[in] data The input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] num_nwave The number of elements in the array @a nwave .
 * @param[in] nwave Wave numbers within the index range of @a data
 * to be used for sinusoidal fitting. The values must be positive
 * or zero (for constant term), but not exceed the maximum wave
 * number specified in creation of @a context . The values must
 * be stored in ascending order and must not be duplicate.
 * @param[in] num_coeff The number of elements in @a coeff .
 * It must be in range (0 < @a num_coeff <= @a num_model_bases ),
 * where @a num_model_bases is ( @a num_nwave*2-1 ) or
 * ( @a num_nwave*2 ) in cases @a nwave contains zero or not,
 * respectively. Also it must not exceed @a num_data .
 * @param[in] coeff Coefficients of model data. Its length must
 * be @a num_coeff. If @a nwave contains zero, the first element
 * must be of the constant term. If not, the first and the second
 * elements are for sine and cosine terms for the smallest wave
 * number in @a nwave, respectively. Coefficients of sine and
 * cosine terms for the second-smallest wave number come next,
 * and so on.
 * @n must-be-aligned
 * @param[out] out The output data. Its length must be @a num_data .
 * @n must-be-aligned
 * @return Status code.
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractSinusoidFloat)(
		struct LIBSAKURA_SYMBOL(LSQFitContextFloat) const *context,
		size_t num_data, float const data[/*num_data*/], size_t num_nwave,
		size_t const nwave[/*num_nwave*/], size_t num_coeff,
		double const coeff[/*num_coeff*/], float out[/*num_data*/])
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @brief Return the number of basis functions used for least-square fitting.
 * @details
 * @param[in] context A context created by
 * @ref sakura_CreateLSQFitContextPolynomialFloat or
 * @ref sakura_CreateLSQFitContextCubicSplineFloat .
 * @param[in] order Parameter for the specified function.
 * It is the order (for polynomial and Chebyshev polynomial), or
 * the number of pieces (for cubic spline). It must be positive
 * for cubic spline, while other models accept zero value. The
 * value should not exceed the @a order specified in creation of
 * @a context .
 * @param[out] num_coeff Number of basis functions to be used for
 * least-square fitting. This value should be the actual number
 * of normal equations.
 * @return Status code.
 *
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetNumberOfCoefficientsFloat)(
		struct LIBSAKURA_SYMBOL(LSQFitContextFloat) const *context,
		uint16_t order, size_t *num_coeff)
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @brief Copy elements in the @a src array into the @a dst array with flipping elements to reorder
 * as some FFT library expects.
 *
 * @details
 * When you provide @a inner_most_untouched = false, @a elements = {3, 4} and @a src = {
 * @code
 *   1,   2,   3,
 *   4,   5,   6,
 *   7,   8,   9,
 *  10,  11,  12,
 * @endcode
 * }, then you will get @a dst as below.
 * @code
 *   9,   7,   8,
 *  12,  10,  11,
 *   3,   1,   2,
 *   6,   4,   5,
 * @endcode
 * If @a inner_most_untouched = true, then you will get @a dst as below.
 * @code
 *  7,   8,   9,
 * 10,  11,  12,
 *  1,   2,   3,
 *  4,   5,   6,
 * @endcode
 *
 * @param[in] inner_most_untouched If true, the order of the inner most dimension is untouched.
 * @param[in] dims Dimensions of the array @a src and @a dst. In other words, the number of elements in @a elements.
 * @param[in] elements Numbers of elements of each dimension of @a src and @a dst with the inner-to-outer order. For example,
 * when you have @code{.cpp}
 * float matrix[4][3];
 * @endcode
 * this parameter should be @code
 * { 3, 4 }
 * @endcode
 * @param[in] src	Source multidimensional array.
 * @n must-be-aligned
 * @param[in] dst	Destination multidimensional array.
 * @n must-be-aligned
 * @return Status code.
 *
 * MT-safe
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(FlipArrayFloat)(
		bool inner_most_untouched, size_t dims, size_t const elements[],
		float const src[], float dst[]) LIBSAKURA_NOEXCEPT;

/**
 * @brief Copy elements in the @a src array into the @a dst array with unflipping elements to the original order.
 *
 * @details
 * When you provide @a inner_most_untouched = false, @a elements = {3, 4} and @a src = {
 * @code
 *   9,   7,   8,
 *  12,  10,  11,
 *   3,   1,   2,
 *   6,   4,   5,
 * @endcode
 * }, then you will get @a dst as below.
 * @code
 *   1,   2,   3,
 *   4,   5,   6,
 *   7,   8,   9,
 *  10,  11,  12,
 * @endcode
 *
 * @param[in] inner_most_untouched If true, the order of the inner most dimension is untouched.
 * @param[in] dims Dimensions of the array @a src and @a dst. In other words, the number of elements in @a elements.
 * @param[in] elements Numbers of elements of each dimension of @a src and @a dst with the inner-to-outer order. For example,
 * when you have @code{.cpp}
 * float matrix[4][3];
 * @endcode
 * this parameter should be @code
 * { 3, 4 }
 * @endcode
 * @param[in] src	Source multidimensional array.
 * @n must-be-aligned
 * @param[in] dst	Destination multidimensional array.
 * @n must-be-aligned
 * @return Status code.
 *
 * MT-safe
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UnflipArrayFloat)(
		bool inner_most_untouched, size_t dims, size_t const elements[],
		float const src[], float dst[]) LIBSAKURA_NOEXCEPT;

/**
 * @copydoc sakura_FlipArrayFloat
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(FlipArrayDouble)(
		bool inner_most_untouched, size_t dims, size_t const elements[],
		double const src[], double dst[]) LIBSAKURA_NOEXCEPT;

/**
 * @copydoc sakura_UnflipArrayFloat
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UnflipArrayDouble)(
		bool inner_most_untouched, size_t dims, size_t const elements[],
		double const src[], double dst[]) LIBSAKURA_NOEXCEPT;

/**
 * @brief Same as @ref sakura_FlipArrayFloat except the element type of the multidimensional arrays.
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(FlipArrayDouble2)(
		bool inner_most_untouched, size_t dims, size_t const elements[],
		double const src[][2], double dst[][2]) LIBSAKURA_NOEXCEPT;

/**
 * @brief Same as @ref sakura_UnflipArrayFloat except the element type of the multidimensional arrays.
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UnflipArrayDouble2)(
		bool inner_most_untouched, size_t dims, size_t const elements[],
		double const src[][2], double dst[][2]) LIBSAKURA_NOEXCEPT;

/**
 * @brief Create mask
 *
 * @details
 * Based on two dimensional position distribution given by @a x and @a y,
 * CreateMaskNearEdgeDouble creates mask for each data point. The algorithm
 * is as follows:
 *
 * -# Conversion stage\n
 * Position distribution given by @a x and @a y are converted to pixel coordinate.
 * User is able to specify pixel size via @a pixel_size. If @a pixel_size is zero,
 * it is set to a half of median separation between neighboring two positions.
 * If @a num_data is odd, there are several choices for median separation since
 * there are number of separation between two positions are even. In such case,
 * the median separation is evaluated as
 * @verbatim
 median_separation = sorted_separation[(num_data - 1) / 2];
 @endverbatim
 * It is recommended to set @a pixel_size to zero. Please consider to set nonzero value only if
 * zero @a pixel_size doesn't work on your data.
 * Number of pixels along horizontal and vertical axes is determined by
 * @verbatim
 (ceil(wx / pixel_size), ceil(wy / pixel_size))
 @endverbatim
 * where wx and wy are evaluated from @a blc_x, @a blc_y, @a trc_x, and @a trc_y if they are
 * all non-NULL. The formula is
 * @verbatim
 wx = (*trc_x - *blc_x) * 1.1
 wy = (*trc_y - *blc_y) * 1.1
 @endverbatim
 * Otherwise, they are calculated by
 * @verbatim
 wx = (max(x) - min(x)) * 1.1
 wy = (max(y) - min(y)) * 1.1
 @endverbatim
 *
 * -# Count stage\n
 * Count up data points in each pixel.
 *
 * -# Binalization stage\n
 * Binarize pixel data. Pixel value is set to one if it has non-zero value. Zero remains zero.
 *
 * -# Fill stage\n
 * Fill zero pixels bracketed by one with one. It is repeated until there is no update on
 * pixel values.
 *
 * -# Detection stage\n
 * Find the area that has non-zero value, detect edge of the area, and mask data points
 * that belong to edges, i.e., mask for data in the edge are set to true. The process is
 * repeated until number of masked position exceeds the threshold given by
 * @a num_data times @a fraction.
 *
 * @remark Note that threshold for the detection stage will become 0 if @a num_data times
 * @a fraction is less than 1. No data will be masked in that case.
 *
 * @remark Note also that @a x and @a y are intended to represent any trajectory in the two-dimensional
 * plane (e.g., regularly sampled two-dimensional position data of a certain object).
 * Thus, orders of @a x and @a y are important.
 *
 * @param[in] fraction Fraction of the data to be masked. Threshold for mask operation
 * is evaluated by @a fraction times @a num_data. 0 <= @a fraction <= 1.
 * @param[in] pixel_size Size of the pixel. If it is zero, pixel size is set to a half of
 * median of the distance between neighboring two positions that are provided by @a x and @a y.
 * @a pixel_size must be greater than or equal to 0. Setting 0 works most of the cases.
 * @param[in] num_data Number of data.
 * @param[in] x List of horizontal positions. The length must be @a num_data.\n
 * must-be-aligned
 * @param[in] y List of vertical positions. The length must be @a num_data.\n
 * must-be-aligned
 * @param[in] blc_x Horizontal limit for bottom-left corner of user-defined range.
 * If non-NULL value is provided, the data having @a x value less than @a *blc_x will be ignored.
 * NULL corresponds to no limit.
 * @param[in] blc_y Vertical limit for bottom-left corner of user-defined range.
 * If non-NULL value is provided, the data having @a y value less than @a *blc_y will be ignored.
 * NULL corresponds to no limit.
 * @param[in] trc_x Horizontal limit for top-right corner of user-defined range.
 * If non-NULL value is provided, the data having @a x value greater than @a *trc_x will be ignored.
 * NULL corresponds to no limit.
 * @param[in] trc_y Vertical limit for top-right corner of user-defined range.
 * If non-NULL value is provided, the data having @a y value greater than @a *trc_y will be ignored.
 * NULL corresponds to no limit.
 * @param[out] mask Output mask. The length must be @a num_data since @a mask is associating
 * with each data points provided by @a x and @a y. The points near edge will be masked,
 * i.e., mask value will be set to true for those points.\n
 * must-be-aligned
 *
 * @return Status code.
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateMaskNearEdgeDouble)(
		float fraction, double pixel_size, size_t num_data, double const x[],
		double const y[], double const *blc_x, double const *blc_y,
		double const *trc_x, double const *trc_y, bool mask[]);

#ifdef __cplusplus
}
/* extern "C" */
#endif

#endif /* LIBSAKURA_LIBSAKURA_SAKURA_H_ */
