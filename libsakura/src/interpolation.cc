#include <iostream>

#include <libsakura/sakura.h>
#include <libsakura/optimized_implementation_factory_impl.h>
#include <libsakura/localdef.h>

#if defined(__AVX__)

namespace {

} /* anonymous namespace */

#else

namespace {

} /* anonymous namespace */

#endif

namespace LIBSAKURA_PREFIX {

void ADDSUFFIX(Interpolation, ARCH_SUFFIX)::Interpolate1dFloatNearest(size_t num_base,
		double const x_base[/*num_base*/], float const y_base[/*num_base*/],
		size_t num_interpolated,
		double x_interpolated[/*num_interpolated*/],
		float y_interpolated[/*num_interpolated*/]) const {
	std::cout << "Interpolate1dFloatNearest is called" << std::endl;
}

void ADDSUFFIX(Interpolation, ARCH_SUFFIX)::Interpolate1dFloatLinear(size_t num_base,
		double const x_base[/*num_base*/], float const y_base[/*num_base*/],
		size_t num_interpolated,
		double x_interpolated[/*num_interpolated*/],
		float y_interpolated[/*num_interpolated*/]) const {
	std::cout << "Interpolate1dFloatLinear is called" << std::endl;
}

void ADDSUFFIX(Interpolation, ARCH_SUFFIX)::Interpolate1dFloatPolynomial(int polynomial_order,
		size_t num_base,
		double const x_base[/*num_base*/], float const y_base[/*num_base*/],
		size_t num_interpolated,
		double x_interpolated[/*num_interpolated*/],
		float y_interpolated[/*num_interpolated*/]) const {
	std::cout << "Interpolate1dFloatPolynomial is called" << std::endl;
}

void ADDSUFFIX(Interpolation, ARCH_SUFFIX)::Interpolate1dFloatSpline(size_t num_base,
		double const x_base[/*num_base*/], float const y_base[/*num_base*/],
		size_t num_interpolated,
		double x_interpolated[/*num_interpolated*/],
		float y_interpolated[/*num_interpolated*/]) const {
	std::cout << "Interpolate1dFloatSpline is called" << std::endl;
}

} /* namespace LIBSAKURA_PREFIX */
