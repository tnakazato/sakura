#include <cassert>
#include <iostream>

#include <libsakura/sakura.h>
#include <libsakura/optimized_implementation_factory_impl.h>
#include <libsakura/localdef.h>

#include <Eigen/Core>

using ::Eigen::Map;
using ::Eigen::Array;
using ::Eigen::Dynamic;
using ::Eigen::Aligned;

namespace LIBSAKURA_PREFIX {

template<typename DataType>
void ADDSUFFIX(Interpolation, ARCH_SUFFIX)<DataType>::Interpolate1dNearest(size_t num_base,
		double const x_base[/* num_base */], DataType const y_base[/* num_base */],
		size_t num_interpolated,
		double const x_interpolated[/* num_interpolated */],
		DataType y_interpolated[/*num_interpolated*/]) const
{
	assert(num_base > 0);
	assert(num_interpolated > 0);
	assert(LIBSAKURA_SYMBOL(IsAligned)(x_base));
	assert(LIBSAKURA_SYMBOL(IsAligned)(y_base));
	assert(LIBSAKURA_SYMBOL(IsAligned)(x_interpolated));
	assert(LIBSAKURA_SYMBOL(IsAligned)(y_interpolated));

	if (num_base == 1) {
		Map<Array<DataType, Dynamic, 1>, Aligned> y_interpolated_vector(y_interpolated, num_interpolated);
		y_interpolated_vector.setConstant(y_base[0]);
		return;
	}
	else {
		int location = 0;
		int end_position = static_cast<int>(num_base) - 1;
		for (size_t index = 0; index < num_interpolated; ++index) {
			location = this->Locate(location, end_position, num_base,
					x_base, x_interpolated[index]);
			if (location == 0) {
				y_interpolated[index] = y_base[0];
			}
			else if (location == static_cast<int>(num_base)) {
				y_interpolated[index] = y_base[location-1];
			}
			else {
				double dx_left = fabs(x_interpolated[index] - x_base[location-1]);
				double dx_right = fabs(x_interpolated[index] - x_base[location]);
				if (dx_left <= dx_right) {
					y_interpolated[index] = y_base[location-1];
				}
				else {
					y_interpolated[index] = y_base[location];
				}
			}
		}
	}
}

template<typename DataType>
void ADDSUFFIX(Interpolation, ARCH_SUFFIX)<DataType>::Interpolate1dLinear(size_t num_base,
		double const x_base[/*num_base*/], DataType const y_base[/*num_base*/],
		size_t num_interpolated,
		double const x_interpolated[/*num_interpolated*/],
		DataType y_interpolated[/*num_interpolated*/]) const {
	assert(num_base > 0);
	assert(num_interpolated > 0);
	assert(LIBSAKURA_SYMBOL(IsAligned)(x_base));
	assert(LIBSAKURA_SYMBOL(IsAligned)(y_base));
	assert(LIBSAKURA_SYMBOL(IsAligned)(x_interpolated));
	assert(LIBSAKURA_SYMBOL(IsAligned)(y_interpolated));

	if (num_base == 1) {
		Map<Array<DataType, Dynamic, 1>, Aligned> y_interpolated_vector(y_interpolated, num_interpolated);
		y_interpolated_vector.setConstant(y_base[0]);
		return;
	}
	else {
		int location = 0;
		int end_position = static_cast<int>(num_base) - 1;
		for (size_t index = 0; index < num_interpolated; ++index) {
			location = this->Locate(location, end_position, num_base,
					x_base, x_interpolated[index]);
			if (location == 0) {
				y_interpolated[index] = y_base[0];
			}
			else if (location == static_cast<int>(num_base)) {
				y_interpolated[index] = y_base[location-1];
			}
			else {
				y_interpolated[index] = y_base[location-1] + (y_base[location] - y_base[location-1])
						* (x_interpolated[index] - x_base[location-1])
						/ (x_base[location] - x_base[location-1]);
			}
		}
	}
}

template<typename DataType>
void ADDSUFFIX(Interpolation, ARCH_SUFFIX)<DataType>::Interpolate1dPolynomial(int polynomial_order,
		size_t num_base,
		double const x_base[/*num_base*/], DataType const y_base[/*num_base*/],
		size_t num_interpolated,
		double const x_interpolated[/*num_interpolated*/],
		DataType y_interpolated[/*num_interpolated*/]) const {
	assert(num_base > 0);
	assert(num_interpolated > 0);
	assert(LIBSAKURA_SYMBOL(IsAligned)(x_base));
	assert(LIBSAKURA_SYMBOL(IsAligned)(y_base));
	assert(LIBSAKURA_SYMBOL(IsAligned)(x_interpolated));
	assert(LIBSAKURA_SYMBOL(IsAligned)(y_interpolated));

	if (num_base == 1) {
		Map<Array<DataType, Dynamic, 1>, Aligned> y_interpolated_vector(y_interpolated, num_interpolated);
		y_interpolated_vector.setConstant(y_base[0]);
		return;
	}
	else if (polynomial_order + 1 >= static_cast<int>(num_base)){
		// use full region for interpolation
		// call polynomial interpolation
	}
	else {
		// use sub-region around the most nearest points
		int location = 0;
		int end_position = static_cast<int>(num_base) - 1;
		for (size_t index = 0; index < num_interpolated; ++index) {
			location = this->Locate(location, end_position, num_base,
					x_base, x_interpolated[index]);
			if (location == 0) {
				y_interpolated[index] = y_base[0];
			}
			else if (location == static_cast<int>(num_base)) {
				y_interpolated[index] = y_base[location-1];
			}
			else {
				int j = static_cast<int>(index) - 1 - polynomial_order / 2;
				unsigned int m = static_cast<unsigned int>(num_base) - 1
						- static_cast<unsigned int>(polynomial_order);
				unsigned int k = static_cast<unsigned int>((j > 0) ? j : 0);
				// call polynomial interpolation
			}
		}
	}
}

template<typename DataType>
void ADDSUFFIX(Interpolation, ARCH_SUFFIX)<DataType>::Interpolate1dSpline(size_t num_base,
		double const x_base[/*num_base*/], DataType const y_base[/*num_base*/],
		size_t num_interpolated,
		double const x_interpolated[/*num_interpolated*/],
		DataType y_interpolated[/*num_interpolated*/]) const {
}

template class ADDSUFFIX(Interpolation, ARCH_SUFFIX)<float>;
} /* namespace LIBSAKURA_PREFIX */
