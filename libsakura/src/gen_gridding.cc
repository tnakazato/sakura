#include <cassert>
#include <cstdlib>
#include <cmath>
#include <climits>

#define __STDC_LIMIT_MACROS 1
#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

namespace {
}

#define CHECK_ARGS(x) do { \
	if (!(x)) { \
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument); \
	} \
} while (false)

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GridConvolving)(
		size_t num_spectra, size_t start_spectrum, size_t end_spectrum,
		bool const spectrum_mask[/*num_spectra*/],
		double const x[/*num_spectra*/], double const y[/*num_spectra*/],
		size_t support, size_t sampling, size_t num_polarization,
		uint32_t const polarization_map[/*num_polarization*/],
		size_t num_channels, uint32_t const channel_map[/*num_channels*/],
		bool const mask/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const value/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const weight/*[num_spectra]*/[/*num_channels*/], bool weight_only,
		size_t num_convolution_table/*= ceil(sqrt(2.)*(support+1)*sampling)*/,
		float const convolution_table[/*num_convolution_table*/],
		size_t num_polarization_for_grid, size_t num_channels_for_grid,
		size_t width, size_t height,
		double weight_sum/*[num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float weight_of_grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/]) {
	CHECK_ARGS(spectrum_mask != nullptr);
	CHECK_ARGS(x != nullptr);
	CHECK_ARGS(y != nullptr);
	CHECK_ARGS(polarization_map != nullptr);
	CHECK_ARGS(channel_map != nullptr);
	CHECK_ARGS(mask != nullptr);
	CHECK_ARGS(weight_only || value != nullptr);
	CHECK_ARGS(weight != nullptr);
	CHECK_ARGS(convolution_table != nullptr);
	CHECK_ARGS(weight_sum != nullptr);
	CHECK_ARGS(weight_of_grid != nullptr);
	CHECK_ARGS(grid != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(spectrum_mask));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(x));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(y));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(polarization_map));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(channel_map));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(mask));
	CHECK_ARGS(weight_only || LIBSAKURA_SYMBOL(IsAligned)(value));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(weight));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(convolution_table));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(weight_sum));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(grid));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(weight_of_grid));
	CHECK_ARGS(
			0 <= start_spectrum && start_spectrum <= end_spectrum
					&& end_spectrum <= num_spectra);
	CHECK_ARGS(0 < support && support <= (INT32_MAX - 1) / 2);
	CHECK_ARGS(0 < sampling && sampling <= INT32_MAX);
	CHECK_ARGS(static_cast<size_t>(support) * sampling <= INT32_MAX / 32);
	CHECK_ARGS(0 < num_polarization && num_polarization <= INT32_MAX);
	CHECK_ARGS(0 < num_channels && num_channels <= INT32_MAX);
	CHECK_ARGS(
			static_cast<int>(ceil(sqrt(2.) * (support + 1) * sampling)) <= num_convolution_table && num_convolution_table <= INT32_MAX / 32);
	CHECK_ARGS(
			0 < num_polarization_for_grid && num_polarization_for_grid <= INT32_MAX);
	CHECK_ARGS(0 < num_channels_for_grid && num_channels_for_grid <= INT32_MAX);
	CHECK_ARGS(0 < width && width <= INT32_MAX);
	CHECK_ARGS(0 < height && height <= INT32_MAX);
	CHECK_ARGS(static_cast<size_t>(width) * height <= INT32_MAX);
#if !defined(NDEBUG)
	for (size_t i = 0; i < num_polarization; ++i) {
		CHECK_ARGS(
				0 <= polarization_map[i]
				&& polarization_map[i] < num_polarization_for_grid);
	}
	for (size_t i = 0; i < num_channels; ++i) {
		CHECK_ARGS(
				0 <= channel_map[i] && channel_map[i] < num_channels_for_grid);
	}
#endif /* !defined(NDEBUG) */

	try {
		auto gridImpl =
				::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetGriddingImpl();
		gridImpl->GridConvolving(num_spectra, start_spectrum, end_spectrum,
				spectrum_mask, x, y, support, sampling, num_polarization,
				polarization_map, num_channels, channel_map, mask, value,
				weight, weight_only, num_convolution_table, convolution_table,
				num_polarization_for_grid, num_channels_for_grid, width, height,
				weight_sum, weight_of_grid, grid);
	} catch (...) {
		assert(false); // no exception should not be raised for the current implementation.
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}
