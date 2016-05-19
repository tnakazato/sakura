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
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <memory>

#include <libsakura/config.h>
#include <libsakura/localdef.h>
#include <libsakura/logger.h>
#include <libsakura/sakura.h>
#include <libsakura/memory_manager.h>

namespace {
// a logger for this module
auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("mask_edge");

template<typename DataType>
void InitializeArray(size_t n, DataType value, DataType data[]) {
	for (size_t i = 0; i < n; ++i) {
		data[i] = value;
	}
}

/**
 * @brief Convert location of data points to the one in pixel coordinate.
 *
 * @param[in] pixel_size pixel size. If it is zero, pixel size is evaluated as
 * a half of median separation between two neighboring data points.
 * @param[in] num_data number of data points
 * @param[in] x x-axis coordinate value of data points. Its length must be @a num_data.
 * must-be-aligned
 * @param[in] y y-axis coordinate value of data points. Its length must be @a num_data.
 * must-be-aligned
 * @param[out] pixel_x x-axis coordinate value of data points in pixel coordinate.
 * Its length must be @a num_data. must-be-aligned.
 * @param[out] pixel_y y-axis coordinate value of data points in pixel coordinate.
 * Its length must be @a num_data. must-be-aligned.
 * @param[out] center_x x-axis coordinate value of central (reference) point for
 * conversion to pixel coordinate.
 * @param[out] center_y y-axis coordinate value of central (reference) point for
 * conversion to pixel coordinate.
 * @param[out] num_horizontal number of pixels along x-axis.
 * @param[out] num_vertical number of pixels along y-axis.
 * @param[out] pixel_width resulting width of the pixel. Pixels are always square.
 */
template<typename DataType>
inline LIBSAKURA_SYMBOL(Status) ConvertToPixel(DataType pixel_size,
		size_t num_data, DataType const x[], DataType const y[],
		DataType pixel_x[], DataType pixel_y[], DataType const *blc_x_in,
		DataType const *blc_y_in, DataType const *trc_x_in,
		DataType const *trc_y_in, DataType *center_x, DataType *center_y,
		size_t *num_horizontal, size_t *num_vertical, DataType *pixel_width) {
	// To derive median separation between two neighboring data points
	// use pixel_x and pixel_y as a working storage
	STATIC_ASSERT(sizeof(DataType) >= sizeof(bool));
	for (size_t i = 0; i < num_data - 1; ++i) {
		DataType separation_x = x[i + 1] - x[i];
		DataType separation_y = y[i + 1] - y[i];
		pixel_x[i] = separation_x * separation_x + separation_y * separation_y;
	}

	// 2015/11/19 TN
	// API change: pixel_width is set to pixel_size if pixel_size is nonzero.
	//             If pixel_size is zero, pixel_width is set based on median separation.
	if (pixel_size == 0.0) {
		// sort data to evaluate median
		std::qsort(pixel_x, num_data - 1, sizeof(DataType),
				[](const void *a, const void *b) {
					DataType aa = *static_cast<DataType const *>(a);
					DataType bb = *static_cast<DataType const *>(b);
					if (aa < bb) {
						return -1;
					}
					else if (aa > bb) {
						return 1;
					}
					else {
						return 0;
					}
				});
		DataType median_separation = std::sqrt(pixel_x[(num_data - 1) / 2]);

		if (median_separation == 0.0) {
			return LIBSAKURA_SYMBOL(Status_kNG);
		}

		// pixel width is half of median_separation
		*pixel_width = median_separation * 0.5;
	} else {
		*pixel_width = pixel_size;
	}

	if (*pixel_width <= 0.0) {
		// Invalid pixel_width is set due to unknown reason
		return LIBSAKURA_SYMBOL(Status_kNG);
	}

	// minimul and maximum position
	DataType blc_x = x[0];
	DataType blc_y = y[0];
	DataType trc_x = x[0];
	DataType trc_y = y[0];
	for (size_t i = 0; i < num_data - 1; ++i) {
		DataType local_x = x[i];
		DataType local_y = y[i];
		blc_x = std::min(blc_x, local_x);
		trc_x = std::max(trc_x, local_x);
		blc_y = std::min(blc_y, local_y);
		trc_y = std::max(trc_y, local_y);
	}
	if (blc_x_in != nullptr) {
		blc_x = *blc_x_in;
	}
	if (blc_y_in != nullptr) {
		blc_y = *blc_y_in;
	}
	if (trc_x_in != nullptr) {
		trc_x = *trc_x_in;
	}
	if (trc_y_in != nullptr) {
		trc_y = *trc_y_in;
	}
	//std::cout << "blc = [" << blc_x << "," << blc_y << "] trc = [" << trc_x << "," << trc_y << "]" << std::endl;
	if (blc_x >= trc_x || blc_y >= trc_y) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}

	// number of pixels
	constexpr DataType kMargin = DataType(1.1);
	DataType const wx = (trc_x - blc_x) * kMargin;
	DataType const wy = (trc_y - blc_y) * kMargin;
	*num_horizontal = static_cast<size_t>(std::ceil(wx / *pixel_width));
	*num_vertical = static_cast<size_t>(std::ceil(wy / *pixel_width));

	// center coordinate
	constexpr DataType kTwo = DataType(2);
	*center_x = (blc_x + trc_x) / kTwo;
	*center_y = (blc_y + trc_y) / kTwo;

	// compute pixel_x and pixel_y
	DataType pixel_center_x = static_cast<DataType>(*num_horizontal - 1) / kTwo;
	DataType pixel_center_y = static_cast<DataType>(*num_vertical - 1) / kTwo;
	constexpr DataType kOutOfRangeValue = -1.0;
	for (size_t i = 0; i < num_data; ++i) {
		double xi = x[i];
		double yi = y[i];
		if (blc_x <= xi && xi <= trc_x) {
			pixel_x[i] = pixel_center_x + (x[i] - *center_x) / *pixel_width;
		} else {
			pixel_x[i] = kOutOfRangeValue;
		}
		if (blc_y <= yi && yi <= trc_y) {
			pixel_y[i] = pixel_center_y + (y[i] - *center_y) / *pixel_width;
		} else {
			pixel_y[i] = kOutOfRangeValue;
		}
	}

	return LIBSAKURA_SYMBOL(Status_kOK);
}

/**
 * @brief Count up number of data in each pixel
 *
 * @param[in] num_data number of data points
 * @param[in] pixel_x x-axis coordinate value of data points in pixel coordinate.
 * Its length must be @a num_data. must-be-aligned.
 * @param[in] pixel_y y-axis coordinate value of data points in pixel coordinate.
 * Its length must be @a num_data. must-be-aligned.
 * @param[in] num_horizontal number of pixels along x-axis.
 * @param[in] num_vertical number of pixels along y-axis.
 * @param[out] counts data counts in each pixel. This is a flattened one-dimensional
 * array of two-dimensional array with shape of [@a num_vertical][@a num_horizontal].
 * must-be-aligned
 */
template<typename DataType>
inline LIBSAKURA_SYMBOL(Status) CountUp(size_t num_data,
		DataType const pixel_x[], DataType const pixel_y[],
		size_t num_horizontal, size_t num_vertical, size_t counts[]) {

	// initialization
	InitializeArray(num_horizontal * num_vertical, 0ul, counts);

	// data counts for each pixel
	for (size_t i = 0; i < num_data; ++i) {
		auto pxi = pixel_x[i];
		auto pyi = pixel_y[i];
		if (pxi >= 0.0 && pyi >= 0.0) {
			size_t ix = static_cast<size_t>(std::round(pxi));
			size_t iy = static_cast<size_t>(std::round(pyi));
			assert(ix < num_horizontal);
			assert(iy < num_vertical);
			size_t index = ix * num_vertical + iy;
			assert(index < num_horizontal * num_vertical);
			counts[index] += 1;
		}
	}

	return LIBSAKURA_SYMBOL(Status_kOK);
}

/**
 * @brief Binarize data counts
 *
 * @param[in] num_horizontal number of pixels along x-axis.
 * @param[in] num_vertical number of pixels along y-axis.
 * @param[in,out] counts data counts in each pixel. This is a flattened one-dimensional
 * array of two-dimensional array with shape of [@a num_vertical][@a num_horizontal].
 * On output, counts is binarized so that non-zero vaule is converted to 1.
 * must-be-aligned
 */
inline LIBSAKURA_SYMBOL(Status) Binarize(size_t num_horizontal,
		size_t num_vertical, size_t counts[]) {
	for (size_t i = 0; i < num_horizontal * num_vertical; ++i) {
		counts[i] = (counts[i] > 0) ? 1 : 0;
	}

	return LIBSAKURA_SYMBOL(Status_kOK);
}

template<typename DataType>
inline size_t SearchForward(size_t start, size_t increment, size_t n,
		DataType const mask[]) {
	for (size_t i = 0; i < n; ++i) {
		size_t index = start + i * increment;
		if (mask[index] > 0) {
			return index;
		}
	}
	return start + (n - 1) * increment;
}

template<typename DataType>
inline size_t SearchBackward(size_t start, size_t increment, size_t n,
		DataType const mask[]) {
	for (ssize_t i = 0; i < n; ++i) {
		assert(i * increment <= start);
		size_t index = start - i * increment;
		if (mask[index] > 0) {
			return index;
		}
	}
	return 0;
}

/**
 * @brief Fill closed or approximately closed region
 *
 * Fill (approximately) closed region. To overcome incomplete closure, it iteratively
 * fill the region until number of filled pixels converges.
 *
 * @param[in] num_horizontal number of pixels along x-axis.
 * @param[in] num_vertical number of pixels along y-axis.
 * @param[in,out] pixel_mask binarized data counts. On input, pixels that contains any number
 * of data points have the value 1. Otherwise, the value is 0. On output, iterative fill process
 * is applied. must-be-aligned
 */
inline LIBSAKURA_SYMBOL(Status) Fill(size_t num_max_iteration,
		size_t num_horizontal, size_t num_vertical, size_t pixel_mask[]) {

	size_t num_masked = 1;
	size_t num_iteration = 0;

	do {
		num_masked = 0;
		for (size_t i = 0; i < num_horizontal; ++i) {
			size_t search_start = num_vertical * i;
			size_t increment = 1;
			size_t start = SearchForward(search_start, increment, num_vertical,
					pixel_mask);
			search_start += num_vertical - 1;
			size_t end = SearchBackward(search_start, increment, num_vertical,
					pixel_mask);
			for (size_t j = start; j <= end; ++j) {
				if (pixel_mask[j] == 0) {
					pixel_mask[j] = 1;
					++num_masked;
				}
			}
		}
		for (size_t i = 0; i < num_vertical; ++i) {
			size_t search_start = i;
			size_t increment = num_vertical;
			size_t start = SearchForward(search_start, increment,
					num_horizontal, pixel_mask);
			search_start += (num_horizontal - 1) * increment;
			size_t end = SearchBackward(search_start, increment, num_horizontal,
					pixel_mask);
			for (size_t j = start; j <= end; j += increment) {
				if (pixel_mask[j] == 0) {
					pixel_mask[j] = 1;
					++num_masked;
				}
			}
		}
		++num_iteration;
	} while (0 < num_masked && num_iteration < num_max_iteration);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

inline void DetectEdge(size_t num_horizontal, size_t num_vertical,
		size_t const pixel_mask[], size_t edge[]) {
	for (size_t i = 0; i < num_vertical; ++i) {
		edge[i] = pixel_mask[i];
	}
	size_t offset = num_vertical * (num_horizontal - 1);
	for (size_t i = 0; i < num_vertical; ++i) {
		edge[offset + i] = pixel_mask[offset + i];
	}
	for (size_t i = 1; i < num_horizontal - 1; ++i) {
		edge[i * num_vertical] = pixel_mask[i * num_vertical];
		edge[(i + 1) * num_vertical - 1] =
				pixel_mask[(i + 1) * num_vertical - 1];
	}
	for (size_t i = 1; i < num_horizontal - 1; ++i) {
		for (size_t j = 1; j < num_vertical - 1; ++j) {
			size_t index = num_vertical * i + j;
			if (pixel_mask[index] == 1) {
				size_t index_below = index - num_vertical;
				size_t index_above = index + num_vertical;
				size_t surroundings = pixel_mask[index - 1]
						* pixel_mask[index + 1] * pixel_mask[index_below - 1]
						* pixel_mask[index_below] * pixel_mask[index_below + 1]
						* pixel_mask[index_above - 1] * pixel_mask[index_above]
						* pixel_mask[index_above + 1];
				edge[index] = (surroundings == 0) ? 1 : 0;
			}
		}
	}
}

inline void TrimEdge(size_t num_horizontal, size_t num_vertical,
		size_t const edge[], size_t pixel_mask[]) {
	for (size_t i = 0; i < num_horizontal * num_vertical; ++i) {
		if (edge[i] > 0) {
			pixel_mask[i] = 0;
		}
	}
}

/**
 * @brief Find data near edge using reference pixel image
 *
 * @param[in] fraction
 * @param[in] num_horizontal
 * @param[in] num_vertical
 * @param[in] pixel_mask pixel mask information. The value is 1 if any number of data
 * points are included in the pixel, or 0 otherwise. It will be destroyed in the function.
 * must-be-aligned
 * @param[in] num_data
 * @param[in] pixel_x
 * @param[in] pixel_y
 * @param[out] mask
 */
template<typename DataType>
inline LIBSAKURA_SYMBOL(Status) DetectDataNearEdge(float fraction,
		size_t num_horizontal, size_t num_vertical, size_t pixel_mask[],
		size_t num_data, DataType const pixel_x[], DataType const pixel_y[],
		bool mask[], size_t *num_masked) {

	assert(fraction <= 1.0f);
	size_t threshold = static_cast<size_t>(static_cast<float>(num_data)
			* fraction);

	if (threshold == 0) {
		// do nothing
		return LIBSAKURA_SYMBOL(Status_kOK);
	}

	size_t num_masked_local = 0;

	// working array to store edge information
	size_t *edge = nullptr;
	size_t num_pixel = num_horizontal * num_vertical;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					num_pixel * sizeof(size_t), &edge));
	InitializeArray(num_pixel, 0ul, edge);

	// iteration loop for detection of data points and edge trimming
	do {
		// detection
		DetectEdge(num_horizontal, num_vertical, pixel_mask, edge);

		// marking
		for (size_t i = 0; i < num_data; ++i) {
			size_t ix = static_cast<size_t>(std::round(pixel_x[i]));
			size_t iy = static_cast<size_t>(std::round(pixel_y[i]));
			if (ix < num_horizontal && iy < num_vertical) {
				assert(ix < num_horizontal);
				assert(iy < num_vertical);
				size_t index = ix * num_vertical + iy;
				assert(index < num_pixel);
				if (mask[i] == false && edge[index] == 1) {
					mask[i] = true;
					++num_masked_local;
				}
			}
		}

		// trimming
		TrimEdge(num_horizontal, num_vertical, edge, pixel_mask);

	} while (num_masked_local < threshold);

	*num_masked = num_masked_local;

	return LIBSAKURA_SYMBOL(Status_kOK);
}

/**
 * Edge detection algorithm consists of the following steps:
 *
 *    ConvertToPixel
 *        set up pixel coordinate based on median separation
 *        between two neighboring points and user-specified
 *        scaling factor, pixel_size, and convert data list
 *        to pixel coordinate.
 *    CountUp
 *        count up number of points contained in each pixel.
 *    Binarize
 *        binarization of countup result such that pixel value
 *        is 1 if there is a data point in this pixel, otherwise
 *        the value is 0.
 *    Fill
 *        set value 1 to pixels that is located between two pixels
 *        with value of 1.
 *    DetectDatanearEdge
 *        find map edge and mask points that are located in edge pixels.
 *        edge finding and masking process is repeated until the fraction
 *        of number of masked points exceeds user-specified fraction.
 */
template<typename DataType>
inline LIBSAKURA_SYMBOL(Status) CreateMaskNearEdge(float fraction,
		DataType pixel_size, size_t num_data, DataType const x[],
		DataType const y[], DataType const *blc_x, DataType const *blc_y,
		DataType const *trc_x, DataType const *trc_y, bool mask[]) {
	// initialize mask
	InitializeArray(num_data, false, mask);

	// do nothing if effective fraction is zero
	if (fraction * static_cast<float>(num_data) < 1.0f) {
		return sakura_Status_kOK;
	}

	// ConvertToPixel
	// allocate storage for two aligned arrays by one allocation call
	STATIC_ASSERT(sizeof(DataType) < LIBSAKURA_ALIGNMENT);
	constexpr size_t kStride = LIBSAKURA_ALIGNMENT / sizeof(DataType);
	size_t const reminder = num_data % kStride;
	size_t const margin = (reminder == 0) ? 0 : (kStride - reminder);
	size_t const size_of_storage = (2 * num_data + margin) * sizeof(DataType);
	DataType *pixel_x = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					size_of_storage, &pixel_x));
	DataType *pixel_y = pixel_x + (num_data + margin);
	DataType center_x = 0.0;
	DataType center_y = 0.0;
	DataType pixel_width = 0.0;
	size_t num_horizontal = 0;
	size_t num_vertical = 0;

	LIBSAKURA_SYMBOL(Status) status = ConvertToPixel(pixel_size, num_data, x, y,
			pixel_x, pixel_y, blc_x, blc_y, trc_x, trc_y, &center_x, &center_y,
			&num_horizontal, &num_vertical, &pixel_width);

	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		return status;
	}

	// CountUp
	size_t *counts = nullptr;
	size_t num_pixel = num_horizontal * num_vertical;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage2(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					num_pixel * sizeof(size_t), &counts));
	status = CountUp(num_data, pixel_x, pixel_y, num_horizontal, num_vertical,
			counts);

	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		return status;
	}

	// Binarize
	status = Binarize(num_horizontal, num_vertical, counts);

	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		return status;
	}

	// Fill
	constexpr size_t kMaxIteration = 10;
	status = Fill(kMaxIteration, num_horizontal, num_vertical, counts);

	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		return status;
	}

	// Detect
	size_t num_masked = 0;
	status = DetectDataNearEdge(fraction, num_horizontal, num_vertical, counts,
			num_data, pixel_x, pixel_y, mask, &num_masked);

	return status;
}

} // anonymous namespace

#define CHECK_ARGS(x) do { \
	if (!(x)) { \
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument); \
	} \
} while (false)

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateMaskNearEdgeDouble)(
		float fraction, double pixel_size, size_t num_data, double const x[],
		double const y[], double const *blc_x, double const *blc_y,
		double const *trc_x, double const *trc_y,
		bool mask[]) {
	// Argument check
	CHECK_ARGS(!std::isnan(fraction));
	CHECK_ARGS(!std::isinf(fraction));
	CHECK_ARGS(0.0 <= fraction && fraction <= 1.0);
	CHECK_ARGS(!std::isnan(pixel_size));
	CHECK_ARGS(!std::isinf(pixel_size));
	CHECK_ARGS(0.0 <= pixel_size);
	CHECK_ARGS(num_data != 1 || (num_data == 1 && pixel_size > 0.0));
	CHECK_ARGS(x != nullptr);
	CHECK_ARGS(y != nullptr);
	CHECK_ARGS(mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(x));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(y));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(mask));
	for (size_t i = 0; i < num_data; ++i) {
		CHECK_ARGS(std::isfinite(x[i]));
	}
	for (size_t i = 0; i < num_data; ++i) {
		CHECK_ARGS(std::isfinite(y[i]));
	}
	CHECK_ARGS((blc_x == nullptr || std::isfinite(*blc_x)));
	CHECK_ARGS((trc_x == nullptr || std::isfinite(*trc_x)));
	CHECK_ARGS((blc_y == nullptr || std::isfinite(*blc_y)));
	CHECK_ARGS((trc_y == nullptr || std::isfinite(*trc_y)));

	try {
		return CreateMaskNearEdge(fraction, pixel_size, num_data, x, y, blc_x,
				blc_y, trc_x, trc_y, mask);
	} catch (const std::bad_alloc &e) {
		// failed to allocate memory
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (...) {
		// any exception is thrown during interpolation
		assert(false);
		LOG4CXX_ERROR(logger, "Aborted due to unknown error");
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}

}
