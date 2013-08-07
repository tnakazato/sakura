#include <iostream>
#include <cassert>
#include <cstdlib>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Interpolate1dFloat)(
		LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method, int polynomial_order,
		size_t num_base, double const x_base[], float const y_base[],
		size_t num_interpolated, double x_interpolated[], float y_interpolated[]) {

	// get object optimized to run-time environment
	auto interpolator =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetInterpolationImpl();

	switch (interpolation_method) {
	case LIBSAKURA_SYMBOL(InterpolationMethod_kNearest):
		interpolator->Interpolate1dFloatNearest(num_base, x_base, y_base,
				num_interpolated, x_interpolated, y_interpolated);
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kLinear):
		interpolator->Interpolate1dFloatLinear(num_base, x_base, y_base,
				num_interpolated, x_interpolated, y_interpolated);
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kPolynomial):
		interpolator->Interpolate1dFloatPolynomial(polynomial_order,
				num_base, x_base, y_base,
				num_interpolated, x_interpolated, y_interpolated);
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kSpline):
		interpolator->Interpolate1dFloatSpline(num_base, x_base, y_base,
				num_interpolated, x_interpolated, y_interpolated);
		break;
	default:
		std::cerr << "ERROR: Invalid interpolation method type" << std::endl;
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}
