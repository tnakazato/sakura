#include <iostream>
#include <cassert>
#include <cstdlib>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/localdef.h"
#include "libsakura/optimized_implementation_factory_impl.h"

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Interpolate1dFloat)(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		int polynomial_order, size_t num_base, double const x_base[],
		float const y_base[], size_t num_interpolated,
		double const x_interpolated[], float y_interpolated[]) {

	// get object optimized to run-time environment
	auto interpolator =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetInterpolationImpl();

	return interpolator->Interpolate1d(interpolation_method, polynomial_order,
			num_base, x_base, y_base, num_interpolated, x_interpolated,
			y_interpolated);
}

namespace {

template<typename DataType>
int LocateData(int start_position, int end_position, size_t num_base,
		DataType const x_base[], DataType x_located) {
	assert(end_position < static_cast<int>(num_base));
	assert(x_base != nullptr);

	// If length of the array is just 1, return 0
	if (num_base == 1)
		return 0;

	assert(LIBSAKURA_SYMBOL(IsAligned)(x_base));

	int left_index = start_position;
	int right_index = static_cast<int>(num_base);
	if (x_base[0] < x_base[num_base - 1]) {
		// ascending order
		if (x_located <= x_base[0]) {
			// out of range
			return 0;
		} else if (x_located > x_base[num_base - 1]) {
			// out of range
			return num_base;
		} else if (x_located < x_base[start_position]) {
			// x_located is not in the range (start_position, end_position)
			// call this function to search other location
			return LocateData(0, start_position, num_base, x_base, x_located);
		} else if (x_located > x_base[end_position]) {
			// x_located is not in the range (start_position, end_position)
			// call this function to search other location
			return LocateData(end_position, static_cast<int>(num_base - 1),
					num_base, x_base, x_located);
		} else {
			// do bisection
			int left_index = start_position;
			int right_index = end_position;
			while (right_index - left_index > 1) {
				int middle_index = (right_index + left_index) / 2;
				if (x_located > x_base[middle_index]) {
					left_index = middle_index;
				} else {
					right_index = middle_index;
				}
			}
			return right_index;
		}
	} else {
		// descending order
		if (x_located >= x_base[0]) {
			// out of range
			return 0;
		} else if (x_located < x_base[num_base - 1]) {
			// out of range
			return num_base;
		} else if (x_located > x_base[start_position]) {
			// x_located is not in the range (start_position, end_position)
			// call this function to search other location
			return LocateData(0, start_position, num_base, x_base, x_located);
		} else if (x_located < x_base[end_position]) {
			// x_located is not in the range (start_position, end_position)
			// call this function to search other location
			return LocateData(end_position, static_cast<int>(num_base - 1),
					num_base, x_base, x_located);
		} else {
			// do bisection
			int left_index = start_position;
			int right_index = end_position;
			while (right_index - left_index > 1) {
				int middle_index = (right_index + left_index) / 2;
				if (x_located < x_base[middle_index]) {
					left_index = middle_index;
				} else {
					right_index = middle_index;
				}
			}
			return right_index;
		}
	}
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {

// Returns right hand side index of the range that brackets x_located.
// For example, x_located locates between x_base[i] and x_base[i+1],
// Locate returns i+1. If x_located is out of range, Locate returns
// either 0 or num_base depending on which side x_located locates w.r.t.
// x_base array.
template<typename DataType>
int InterpolationImpl<DataType>::Locate(int start_position, int end_position,
		size_t num_base, double const x_base[/*num_base*/],
		double x_located) const {
	return ::LocateData<double>(start_position, end_position, num_base, x_base,
			x_located);
}

template class InterpolationImpl<float> ;
} /* LIBSAKURA_PREFIX */
