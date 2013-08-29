#include <cassert>
#include <iostream>
#include <memory>
#include <cstdalign>

#include <libsakura/sakura.h>
#include <libsakura/optimized_implementation_factory_impl.h>
#include <libsakura/localdef.h>

//#include <Eigen/Core>
//
//using ::Eigen::Map;
//using ::Eigen::Array;
//using ::Eigen::Dynamic;
//using ::Eigen::Aligned;

namespace {

using namespace LIBSAKURA_PREFIX;

template<typename DataType>
int DoNevillePolynomial(double const x_base[], DataType const y_base[],
		unsigned int left_index, int num_elements, double x_interpolated,
		DataType &y_interpolated) {

	// working pointers
	double const *x_ptr = &x_base[left_index];
	DataType const *y_ptr = &y_base[left_index];

	// storage for C and D in Neville's algorithm
	size_t sakura_alignment = LIBSAKURA_SYMBOL(GetAlignment)();
	DataType *storage_for_c = new DataType[num_elements + sakura_alignment - 1];
	DataType *storage_for_d = new DataType[num_elements + sakura_alignment - 1];
	DataType *c =
			const_cast<DataType*>(reinterpret_cast<const DataType*>(LIBSAKURA_SYMBOL(AlignAny)(
					sizeof(DataType), storage_for_c, num_elements)));
	DataType *d =
			const_cast<DataType*>(reinterpret_cast<const DataType*>(LIBSAKURA_SYMBOL(AlignAny)(
					sizeof(DataType), storage_for_d, num_elements)));

	for (int i = 0; i < num_elements; ++i) {
		c[i] = y_ptr[i];
		d[i] = y_ptr[i];
	}

	// Neville's algorithm
	y_interpolated = c[0];
	for (int m = 1; m < num_elements; ++m) {
		// Evaluate Cm1, Cm2, Cm3, ... Cm[n-m] and Dm1, Dm2, Dm3, ... Dm[n-m].
		// Those are stored to c[0], c[1], ..., c[n-m-1] and d[0], d[1], ...,
		// d[n-m-1].
		for (int i = 0; i < num_elements - m; ++i) {
			DataType cd = c[i + 1] - d[i];
			double dx = x_ptr[i] - x_ptr[i + m];
			try {
				cd /= static_cast<DataType>(dx);
			} catch (...) {
				delete[] storage_for_c;
				delete[] storage_for_d;
				std::cerr << "x_base has duplicate elements" << std::endl;
				return 1;
			}
			c[i] = (x_ptr[i] - x_interpolated) * cd;
			d[i] = (x_ptr[i + m] - x_interpolated) * cd;
		}

		// In each step, c[0] holds Cm1 which is a correction between
		// P12...m and P12...[m+1]. Thus, the following repeated update
		// corresponds to the route P1 -> P12 -> P123 ->...-> P123...n.
		y_interpolated += c[0];
	}

	delete[] storage_for_c;
	delete[] storage_for_d;

	return 0;
}
}

namespace LIBSAKURA_PREFIX {

template<typename DataType>
void ADDSUFFIX(Interpolation, ARCH_SUFFIX)<DataType>::Interpolate1dNearest(
		size_t num_base, double const x_base[/* num_base */],
		DataType const y_base[/* num_base */], size_t num_interpolated,
		double const x_interpolated[/* num_interpolated */],
		DataType y_interpolated[/*num_interpolated*/]) const {
	assert(num_base > 0);
	assert(num_interpolated > 0);
	assert(LIBSAKURA_SYMBOL(IsAligned)(x_base));
	assert(LIBSAKURA_SYMBOL(IsAligned)(y_base));
	assert(LIBSAKURA_SYMBOL(IsAligned)(x_interpolated));
	assert(LIBSAKURA_SYMBOL(IsAligned)(y_interpolated));

	if (num_base == 1) {
//		Map<Array<DataType, Dynamic, 1>, Aligned> y_interpolated_vector(
//				y_interpolated, num_interpolated);
//		y_interpolated_vector.setConstant(y_base[0]);
		for (size_t index = 0; index < num_interpolated; ++index) {
			y_interpolated[index] = y_base[0];
		}
	} else {
		int location = 0;
		int end_position = static_cast<int>(num_base) - 1;
		for (size_t index = 0; index < num_interpolated; ++index) {
			location = this->Locate(location, end_position, num_base, x_base,
					x_interpolated[index]);
			if (location == 0) {
				y_interpolated[index] = y_base[0];
			} else if (location == static_cast<int>(num_base)) {
				y_interpolated[index] = y_base[location - 1];
			} else {
				double dx_left = fabs(
						x_interpolated[index] - x_base[location - 1]);
				double dx_right = fabs(
						x_interpolated[index] - x_base[location]);
				if (dx_left <= dx_right) {
					y_interpolated[index] = y_base[location - 1];
				} else {
					y_interpolated[index] = y_base[location];
				}
			}
		}
	}
}

template<typename DataType>
void ADDSUFFIX(Interpolation, ARCH_SUFFIX)<DataType>::Interpolate1dLinear(
		size_t num_base, double const x_base[/*num_base*/],
		DataType const y_base[/*num_base*/], size_t num_interpolated,
		double const x_interpolated[/*num_interpolated*/],
		DataType y_interpolated[/*num_interpolated*/]) const {
	assert(num_base > 0);
	assert(num_interpolated > 0);
	assert(LIBSAKURA_SYMBOL(IsAligned)(x_base));
	assert(LIBSAKURA_SYMBOL(IsAligned)(y_base));
	assert(LIBSAKURA_SYMBOL(IsAligned)(x_interpolated));
	assert(LIBSAKURA_SYMBOL(IsAligned)(y_interpolated));

	if (num_base == 1) {
//		Map<Array<DataType, Dynamic, 1>, Aligned> y_interpolated_vector(
//				y_interpolated, num_interpolated);
//		y_interpolated_vector.setConstant(y_base[0]);
		for (size_t index = 0; index < num_interpolated; ++index) {
			y_interpolated[index] = y_base[0];
		}
	} else {
		int location = 0;
		int end_position = static_cast<int>(num_base) - 1;
		for (size_t index = 0; index < num_interpolated; ++index) {
			location = this->Locate(location, end_position, num_base, x_base,
					x_interpolated[index]);
			if (location == 0) {
				y_interpolated[index] = y_base[0];
			} else if (location == static_cast<int>(num_base)) {
				y_interpolated[index] = y_base[location - 1];
			} else {
				y_interpolated[index] = y_base[location - 1]
						+ (y_base[location] - y_base[location - 1])
								* (x_interpolated[index] - x_base[location - 1])
								/ (x_base[location] - x_base[location - 1]);
			}
		}
	}
}

template<typename DataType>
void ADDSUFFIX(Interpolation, ARCH_SUFFIX)<DataType>::Interpolate1dPolynomial(
		int polynomial_order, size_t num_base,
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
//		Map<Array<DataType, Dynamic, 1>, Aligned> y_interpolated_vector(
//				y_interpolated, num_interpolated);
//		y_interpolated_vector.setConstant(y_base[0]);
		for (size_t index = 0; index < num_interpolated; ++index) {
			y_interpolated[index] = y_base[0];
		}
	} else if (polynomial_order == 0) {
		// This is special case: 0-th polynomial interpolation acts like nearest interpolation
		Interpolate1dNearest(num_base, x_base, y_base, num_interpolated,
				x_interpolated, y_interpolated);
	} else if (polynomial_order + 1 >= static_cast<int>(num_base)) {
		// use full region for interpolation
		// call polynomial interpolation
		int location = 0;
		int end_position = static_cast<int>(num_base) - 1;
		for (size_t index = 0; index < num_interpolated; ++index) {
			location = this->Locate(location, end_position, num_base, x_base,
					x_interpolated[index]);
			if (location == 0) {
				y_interpolated[index] = y_base[0];
			} else if (location == static_cast<int>(num_base)) {
				y_interpolated[index] = y_base[location - 1];
			} else {
				// call polynomial interpolation
				int status = DoNevillePolynomial<DataType>(x_base, y_base, 0,
						num_base, x_interpolated[index], y_interpolated[index]);
			}
		}
	} else {
		// use sub-region around the nearest points
		int location = 0;
		int end_position = static_cast<int>(num_base) - 1;
		for (size_t index = 0; index < num_interpolated; ++index) {
			location = this->Locate(location, end_position, num_base, x_base,
					x_interpolated[index]);
			if (location == 0) {
				y_interpolated[index] = y_base[0];
			} else if (location == static_cast<int>(num_base)) {
				y_interpolated[index] = y_base[location - 1];
			} else {
				int j = location - 1 - polynomial_order / 2;
				unsigned int m = static_cast<unsigned int>(num_base) - 1
						- static_cast<unsigned int>(polynomial_order);
				unsigned int k = static_cast<unsigned int>((j > 0) ? j : 0);
				// call polynomial interpolation
				int status = DoNevillePolynomial<DataType>(x_base, y_base, k,
						polynomial_order + 1, x_interpolated[index],
						y_interpolated[index]);
			}
		}
	}
}

template<typename DataType>
void ADDSUFFIX(Interpolation, ARCH_SUFFIX)<DataType>::Interpolate1dSpline(
		size_t num_base, double const x_base[/*num_base*/],
		DataType const y_base[/*num_base*/], size_t num_interpolated,
		double const x_interpolated[/*num_interpolated*/],
		DataType y_interpolated[/*num_interpolated*/]) const {
}

template class ADDSUFFIX(Interpolation, ARCH_SUFFIX)<float> ;
} /* namespace LIBSAKURA_PREFIX */
