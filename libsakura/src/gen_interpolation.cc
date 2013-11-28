#include <sstream>
#include <cassert>
#include <cstdlib>
#include <climits>
#include <memory>
#include <new>

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
		LOG4CXX_ERROR(logger, "Invalid interpolation method\n");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		process_data = false;
	}

	// num_base must be non-zero
	if (num_interpolation_axis == 0) {
		LOG4CXX_ERROR(logger, "ERROR: num_base must be >0\n");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		process_data = false;
	}

	// no interpolation will be done
	if (num_interpolated == 0 || num_array == 0) {
		// Nothing to do
		LOG4CXX_INFO(logger, "Nothing has been done since num_interpolated is 0\n");
		*status = LIBSAKURA_SYMBOL(Status_kOK);
		process_data = false;
	}

	// input arrays are not aligned
	if (!LIBSAKURA_SYMBOL(IsAligned)(base)
			|| !LIBSAKURA_SYMBOL(IsAligned)(data_base)
			|| !LIBSAKURA_SYMBOL(IsAligned)(interpolated)
			|| !LIBSAKURA_SYMBOL(IsAligned)(data_interpolated)) {
		LOG4CXX_ERROR(logger, "ERROR: input arrays are not aligned\n");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		process_data = false;
	}

	// input arrays are null
	if (base == nullptr || data_base == nullptr || interpolated == nullptr
			|| data_interpolated == nullptr) {
		LOG4CXX_ERROR(logger, "ERROR: input arrays are null\n");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		process_data = false;
	}
	return process_data;
}

}

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

	// get object optimized to run-time environment
	auto interpolator =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetInterpolationImpl();

	try {
		interpolator->InterpolateXAxis(interpolation_method,
				polynomial_order, num_x_base, x_base, num_y, data_base,
				num_x_interpolated, x_interpolated, data_interpolated);
	} catch (const std::bad_alloc &e) {
		// failed to allocate memory
		LOG4CXX_ERROR(logger, "ERROR: Memory allocation failed.\n");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (...) {
		// any exception is thrown during interpolation
		assert(false);
		LOG4CXX_ERROR(logger, "ERROR: Aborted due to unknown error\n");
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

	// get object optimized to run-time environment
	auto interpolator =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetInterpolationImpl();

	try {
		interpolator->InterpolateYAxis(interpolation_method,
				polynomial_order, num_y_base, y_base, num_x, data_base,
				num_y_interpolated, y_interpolated, data_interpolated);
	} catch (const std::bad_alloc &e) {
		// failed to allocate memory
		LOG4CXX_ERROR(logger, "ERROR: Memory allocation failed.\n");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (...) {
		// any exception is thrown during interpolation
		assert(false);
		LOG4CXX_ERROR(logger, "ERROR: Aborted due to unknown error\n");
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return status;
}
