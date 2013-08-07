#include <cassert>
#include <iostream>

#include <libsakura/sakura.h>
#include <libsakura/optimized_implementation_factory_impl.h>
#include <libsakura/localdef.h>

namespace {
} /* anonymous namespace */

#if defined(__AVX__)

namespace {

template<typename DataType>
int LocateDataSimd(int start_position, int end_position, size_t num_base,
		DataType const x_base[], DataType x_located)
{
	std::cout << "Sorry simd version is not implemented yet." << std::endl;
	return -1;
}

template<typename DataType>
void Interpolate1dNearestSimd(size_t num_base, double const x_base[], DataType const y_base[],
		size_t num_interpolated, double const x_interpolated[], DataType y_interpolated[])
{
	std::cout << "Sorry simd version is not implemented yet." << std::endl;
}

} /* anonymous namespace */

#else

#include <Eigen/Core>

using ::Eigen::Map;
using ::Eigen::Array;
using ::Eigen::Dynamic;
using ::Eigen::Aligned;

namespace {

template<typename DataType>
int LocateDataEigen(int start_position, int end_position, size_t num_base,
		DataType const x_base[], DataType x_located)
{
	assert(end_position < static_cast<int>(num_base));

	// If length of the array is just 1, return 0
	if (num_base == 1)
		return 0;

	assert(LIBSAKURA_SYMBOL(IsAligned)(x_base));

	Map<Array<DataType, Dynamic, 1>, Aligned> x_vector(const_cast<DataType *>(x_base), num_base);

	int left_index = start_position;
	int right_index = static_cast<int>(num_base);
	if (x_vector[0] < x_vector[num_base-1])  {
		// ascending order
		if (x_located <= x_vector[0]) {
			// out of range
			return 0;
		}
		else if (x_located > x_vector[num_base-1]) {
			// out of range
			return num_base;
		}
		else if (x_located < x_vector[start_position]) {
			// x_located is not in the range (start_position, end_position)
			// call this function to search other location
			return LocateDataEigen(0, start_position, num_base, x_base, x_located);
		}
		else if (x_located > x_vector[end_position]) {
			// x_located is not in the range (start_position, end_position)
			// call this function to search other location
			return LocateDataEigen(end_position, static_cast<int>(num_base-1), num_base, x_base, x_located);
		}
		else {
			// do bisection
			int left_index = start_position;
			int right_index = end_position;
			while (right_index - left_index > 1) {
				int middle_index = (right_index + left_index) / 2;
				if (x_located > x_vector[middle_index]) {
					left_index = middle_index;
				}
				else {
					right_index = middle_index;
				}
			}
			return right_index;
		}
	}
	else {
		// descending order
		if (x_located >= x_vector[0]) {
			// out of range
			return 0;
		}
		else if (x_located < x_vector[num_base-1]) {
			// out of range
			return num_base;
		}
		else if (x_located > x_vector[start_position]) {
			// x_located is not in the range (start_position, end_position)
			// call this function to search other location
			return LocateDataEigen(0, start_position, num_base, x_base, x_located);
		}
		else if (x_located < x_vector[end_position]) {
			// x_located is not in the range (start_position, end_position)
			// call this function to search other location
			return LocateDataEigen(end_position, static_cast<int>(num_base-1), num_base, x_base, x_located);
		}
		else {
			// do bisection
			int left_index = start_position;
			int right_index = end_position;
			while (right_index - left_index > 1) {
				int middle_index = (right_index + left_index) / 2;
				if (x_located < x_vector[middle_index]) {
					left_index = middle_index;
				}
				else {
					right_index = middle_index;
				}
			}
			return right_index;
		}
	}
}

template<typename DataType>
void Interpolate1dNearestEigen(size_t num_base, double const x_base[], DataType const y_base[],
		size_t num_interpolated, double const x_interpolated[], DataType y_interpolated[])
{
	assert(num_base > 0);
	assert(num_interpolated > 0);
	assert(LIBSAKURA_SYMBOL(IsAligned)(x_base));
	assert(LIBSAKURA_SYMBOL(IsAligned)(y_base));
	assert(LIBSAKURA_SYMBOL(IsAligned)(x_interpolated));
	assert(LIBSAKURA_SYMBOL(IsAligned)(y_interpolated));

	if (num_base == 1) {
		Map<Array<DataType, Dynamic, 1>, Aligned> y_interpolated_vector(y_interpolated, num_base);
		for (size_t index = 0; index < num_interpolated; ++index) {
			y_interpolated_vector[index] = y_base[0];
		}
		return;
	}
	else {
		int location = 0;
		int end_position = static_cast<int>(num_base) - 1;
		Map<Array<double, Dynamic, 1>, Aligned> x_interpolated_vector(const_cast<double *>(x_interpolated), num_interpolated);
		Map<Array<double, Dynamic, 1>, Aligned> x_base_vector(const_cast<double *>(x_base), num_base);
		Map<Array<DataType, Dynamic, 1>, Aligned> y_interpolated_vector(y_interpolated, num_interpolated);
		Map<Array<DataType, Dynamic, 1>, Aligned> y_base_vector(const_cast<DataType *>(y_base), num_base);
		for (size_t index = 0; index < num_interpolated; ++index) {
			location = LocateDataEigen<double>(location, end_position, num_base,
					x_base, x_interpolated[index]/*x_interpolated_vector[index]*/);
			if (location == 0) {
				y_interpolated_vector[index] = y_base_vector[0];
			}
			else if (location == static_cast<int>(num_base)) {
				y_interpolated_vector[index] = y_base_vector[location-1];
			}
			else {
				double dx_left = fabs(x_interpolated[index] - x_base[location-1]);
				double dx_right = fabs(x_interpolated[index] - x_base[location]);
				if (dx_left < dx_right) {
					y_interpolated_vector[index] = y_base_vector[location-1];
				}
				else {
					y_interpolated_vector[index] = y_base_vector[location];
				}
			}
		}
	}
}


} /* anonymous namespace */

#endif

namespace LIBSAKURA_PREFIX {

void ADDSUFFIX(Interpolation, ARCH_SUFFIX)::Interpolate1dFloatNearest(size_t num_base,
		double const x_base[/*num_base*/], float const y_base[/*num_base*/],
		size_t num_interpolated,
		double const x_interpolated[/*num_interpolated*/],
		float y_interpolated[/*num_interpolated*/]) const {
	std::cout << "Interpolate1dFloatNearest is called" << std::endl;
#if defined(__AVX__)
	Interpolate1dNearestSimd<float>(num_base, x_base, y_base,
			num_interpolated, x_interpolated, y_interpolated);
#else
	Interpolate1dNearestEigen<float>(num_base, x_base, y_base,
			num_interpolated, x_interpolated, y_interpolated);
#endif
}

void ADDSUFFIX(Interpolation, ARCH_SUFFIX)::Interpolate1dFloatLinear(size_t num_base,
		double const x_base[/*num_base*/], float const y_base[/*num_base*/],
		size_t num_interpolated,
		double const x_interpolated[/*num_interpolated*/],
		float y_interpolated[/*num_interpolated*/]) const {
	std::cout << "Interpolate1dFloatLinear is called" << std::endl;
}

void ADDSUFFIX(Interpolation, ARCH_SUFFIX)::Interpolate1dFloatPolynomial(int polynomial_order,
		size_t num_base,
		double const x_base[/*num_base*/], float const y_base[/*num_base*/],
		size_t num_interpolated,
		double const x_interpolated[/*num_interpolated*/],
		float y_interpolated[/*num_interpolated*/]) const {
	std::cout << "Interpolate1dFloatPolynomial is called" << std::endl;
}

void ADDSUFFIX(Interpolation, ARCH_SUFFIX)::Interpolate1dFloatSpline(size_t num_base,
		double const x_base[/*num_base*/], float const y_base[/*num_base*/],
		size_t num_interpolated,
		double const x_interpolated[/*num_interpolated*/],
		float y_interpolated[/*num_interpolated*/]) const {
	std::cout << "Interpolate1dFloatSpline is called" << std::endl;
}

int ADDSUFFIX(Interpolation, ARCH_SUFFIX)::Locate(int start_position, int end_position,
		size_t num_base, double const x_base[], double x_located) const {
#if defined(__AVX__)
	return LocateDataSimd<double>(start_position, end_position, num_base, x_base, x_located);
#else
	return LocateDataEigen<double>(start_position, end_position, num_base, x_base, x_located);
#endif
}

} /* namespace LIBSAKURA_PREFIX */
