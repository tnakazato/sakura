#include <sstream>
#include <cassert>
#include <cstdlib>
#include <climits>
#include <memory>

#include "libsakura/sakura.h"
#include "libsakura/localdef.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/logger.h"

namespace {
// a logger for this module
auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("interpolation");

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
		if (LIBSAKURA_PREFIX::Logger::IsErrorEnabled(logger)) {
			std::ostringstream oss;
			oss << "Invalid interpolation method" << std::endl;
			LIBSAKURA_PREFIX::Logger::Error(logger, oss.str().c_str());
		}
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		process_data = false;
	}

	// num_base must be non-zero
	if (num_interpolation_axis == 0) {
		if (LIBSAKURA_PREFIX::Logger::IsErrorEnabled(logger)) {
			std::ostringstream oss;
			oss << "ERROR: num_base must be > 0" << std::endl;
			LIBSAKURA_PREFIX::Logger::Error(logger, oss.str().c_str());
		}
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		process_data = false;
	}

	// no interpolation will be done
	if (num_interpolated == 0 || num_array == 0) {
		// Nothing to do
		if (LIBSAKURA_PREFIX::Logger::IsInfoEnabled(logger)) {
			std::ostringstream oss;
			oss << "Nothing has been done since num_interpolated is 0"
					<< std::endl;
			LIBSAKURA_PREFIX::Logger::Info(logger, oss.str().c_str());
		}
		*status = LIBSAKURA_SYMBOL(Status_kOK);
		process_data = false;
	}

	// input arrays are not aligned
	if (!LIBSAKURA_SYMBOL(IsAligned)(base)
			|| !LIBSAKURA_SYMBOL(IsAligned)(data_base)
			|| !LIBSAKURA_SYMBOL(IsAligned)(interpolated)
			|| !LIBSAKURA_SYMBOL(IsAligned)(data_interpolated)) {
		if (LIBSAKURA_PREFIX::Logger::IsErrorEnabled(logger)) {
			std::ostringstream oss;
			oss << "ERROR: input arrays are not aligned" << std::endl;
			LIBSAKURA_PREFIX::Logger::Error(logger, oss.str().c_str());
		}
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		process_data = false;
	}

	// input arrays are null
	if (base == nullptr || data_base == nullptr || interpolated == nullptr
			|| data_interpolated == nullptr) {
		if (LIBSAKURA_PREFIX::Logger::IsErrorEnabled(logger)) {
			std::ostringstream oss;
			oss << "ERROR: input arrays are null" << std::endl;
			LIBSAKURA_PREFIX::Logger::Error(logger, oss.str().c_str());
		}
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		process_data = false;
	}
	return process_data;
}

}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Interpolate1DXFloat)(
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

	// get object optimized to run-time environment
	auto interpolator =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetInterpolationImpl();

	try {
		interpolator->Interpolate1DX(interpolation_method,
				polynomial_order, num_x_base, x_base, num_y, data_base,
				num_x_interpolated, x_interpolated, data_interpolated);
	} catch (...) {
		// any exception is thrown during interpolation
		if (LIBSAKURA_PREFIX::Logger::IsErrorEnabled(logger)) {
			std::ostringstream oss;
			oss << "ERROR: Aborted due to unknown error" << std::endl;
			LIBSAKURA_PREFIX::Logger::Error(logger, oss.str().c_str());
		}
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return status;
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Interpolate1DYFloat)(
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

	// get object optimized to run-time environment
	auto interpolator =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetInterpolationImpl();

	try {
		interpolator->Interpolate1DY(interpolation_method,
				polynomial_order, num_y_base, y_base, num_x, data_base,
				num_y_interpolated, y_interpolated, data_interpolated);
	} catch (...) {
		// any exception is thrown during interpolation
		if (LIBSAKURA_PREFIX::Logger::IsErrorEnabled(logger)) {
			std::ostringstream oss;
			oss << "ERROR: Aborted due to unknown error" << std::endl;
			LIBSAKURA_PREFIX::Logger::Error(logger, oss.str().c_str());
		}
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return status;
}

namespace {

template<typename DataType>
size_t LocateData(size_t start_position, size_t end_position, size_t num_base,
		DataType const x_base[], DataType x_located) {
	assert(num_base > 0);
	assert(start_position <= end_position);
	assert(end_position < num_base);
	assert(x_base != nullptr);

	// If length of the array is just 1, return 0
	if (num_base == 1)
		return 0;

	assert(LIBSAKURA_SYMBOL(IsAligned)(x_base));

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
			return LocateData(end_position, num_base - 1, num_base, x_base,
					x_located);
		} else {
			// do bisection
			size_t left_index = start_position;
			size_t right_index = end_position;
			while (right_index > left_index + 1) {
				size_t middle_index = (right_index + left_index) / 2;
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
			return LocateData(end_position, num_base - 1, num_base, x_base,
					x_located);
		} else {
			// do bisection
			size_t left_index = start_position;
			size_t right_index = end_position;
			while (right_index > left_index + 1) {
				size_t middle_index = (right_index + left_index) / 2;
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
template<class XDataType, class YDataType>
size_t InterpolationImpl<XDataType, YDataType>::Locate(size_t start_position,
		size_t end_position, size_t num_base,
		XDataType const x_base[/*num_base*/], XDataType x_located) const {
	return ::LocateData<XDataType>(start_position, end_position, num_base,
			x_base, x_located);
}

template class InterpolationImpl<double, float> ;
} /* LIBSAKURA_PREFIX */
