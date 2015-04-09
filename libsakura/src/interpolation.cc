/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2014
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
#include <cassert>
#include <sstream>
#include <memory>
#include <utility>
#include <vector>

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
#include <libsakura/logger.h>
#include <libsakura/memory_manager.h>

namespace {

// a logger for this module
auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("interpolation");

template<class DataType>
struct StorageAndAlignedPointer {
	StorageAndAlignedPointer() :
			storage(nullptr), pointer(nullptr) {
	}
	virtual ~StorageAndAlignedPointer() {
		if (storage != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(storage);
			storage = nullptr;
		}
	}
	void *storage;
	DataType *pointer;
};

template<class DataType>
inline void AllocateAndAlign(size_t num_elements,
		StorageAndAlignedPointer<DataType> *holder) {
	holder->storage = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
			num_elements * sizeof(DataType), &(holder->pointer));
}

// Locator
template<class XDataType>
size_t LocateData(size_t start_position, size_t end_position, size_t num_base,
		XDataType const base_position[], XDataType located_position) {
	assert(num_base > 0);
	assert(start_position <= end_position);
	assert(end_position < num_base);
	assert(base_position != nullptr);

	// If length of the array is just 1, return 0
	if (num_base == 1)
		return 0;

	// base_position must be sorted in ascending order
	if (located_position <= base_position[0]) {
		// out of range
		return 0;
	} else if (located_position > base_position[num_base - 1]) {
		// out of range
		return num_base;
	} else if (located_position < base_position[start_position]) {
		// x_located is not in the range (start_position, end_position)
		// call this function to search other location
		return LocateData(0, start_position, num_base, base_position,
				located_position);
	} else if (located_position > base_position[end_position]) {
		// located_position is not in the range (start_position, end_position)
		// call this function to search other location
		return LocateData(end_position, num_base - 1, num_base, base_position,
				located_position);
	} else {
		// do bisection
		size_t left_index = start_position;
		size_t right_index = end_position;
		while (right_index > left_index + 1) {
			size_t const middle_index = (right_index + left_index) / 2;
			if (located_position > base_position[middle_index]) {
				left_index = middle_index;
			} else {
				right_index = middle_index;
			}
		}
		return right_index;
	}
}

template<class XDataType>
size_t Locate(size_t num_base, size_t num_located,
		XDataType const base_position[], XDataType const located_position[],
		size_t location_list[]) {
	// input arrays must be sorted in ascending order
	size_t num_location_list = 0;
	if (num_base == 1) {
		if (located_position[num_located - 1] <= base_position[0]) {
			num_location_list = 1;
			location_list[0] = 0;
		} else if (located_position[0] >= base_position[0]) {
			num_location_list = 1;
			location_list[0] = 1;
		} else {
			num_location_list = 2;
			location_list[0] = 0;
			location_list[1] = 1;
		}
	} else if (located_position[num_located - 1] <= base_position[0]) {
		// all located_position's are on the left (lower) side of base_position
		num_location_list = 1;
		location_list[0] = 0;
	} else if (located_position[0] >= base_position[num_base - 1]) {
		// all located_position's are on the right (upper) side of base_position
		num_location_list = 1;
		location_list[0] = num_base;
	} else {
		size_t start_position = 0;
		size_t end_position = num_base - 1;
		size_t previous_location = num_base + 1;
		for (size_t i = 0; i < num_located && start_position <= end_position;
				++i) {
			size_t const location = LocateData(start_position, end_position,
					num_base, base_position, located_position[i]);
			if (location != previous_location) {
				location_list[num_location_list] = location;
				num_location_list += 1;
				start_position = location;
				previous_location = location;
			}
		}
	}
	return num_location_list;
}

// Initializer
template<class XDataType, class YDataType>
struct NullWorkingData {
	static NullWorkingData<XDataType, YDataType> * Initialize(
			uint8_t polynomial_order, size_t num_base, size_t num_array,
			XDataType const base_position[], YDataType const base_data[]) {
		return nullptr;
	}
};

template<class DataType>
struct PolynomialWorkingData {
	void InitializePolynomial(uint8_t order, size_t num_base, size_t array_size,
			size_t num_array) {
		polynomial_order =
				(order + 1u >= num_base) ?
						static_cast<uint8_t>(num_base - 1) : order;
		num_elements = polynomial_order + 1;
		xholder.resize(array_size);
		for (size_t i = 0; i < array_size; ++i) {
			AllocateAndAlign<DataType>(num_elements * num_array, &(xholder[i]));
		}
	}
	uint8_t polynomial_order;
	size_t num_elements;
	std::vector<StorageAndAlignedPointer<DataType> > xholder;
};

template<class XDataType, class YDataType>
struct PolynomialXWorkingData: public PolynomialWorkingData<XDataType> {
	static PolynomialXWorkingData<XDataType, YDataType> * Initialize(
			uint8_t polynomial_order, size_t num_base, size_t num_array,
			XDataType const base_position[], YDataType const base_data[]) {
		PolynomialXWorkingData<XDataType, YDataType> *work_data =
				new PolynomialXWorkingData<XDataType, YDataType>();
		work_data->InitializePolynomial(polynomial_order, num_base, 2, 1);
		return work_data;
	}
};

template<class XDataType, class YDataType>
struct PolynomialYWorkingData: public PolynomialWorkingData<XDataType> {
	static PolynomialYWorkingData<XDataType, YDataType> * Initialize(
			uint8_t polynomial_order, size_t num_base, size_t num_array,
			XDataType const base_position[], YDataType const base_data[]) {
		PolynomialYWorkingData<XDataType, YDataType> *work_data =
				new PolynomialYWorkingData<XDataType, YDataType>();
		work_data->InitializePolynomial(polynomial_order, num_base, 3,
				num_array);
		return work_data;
	}
};

template<class XDataType, class YDataType>
inline void DeriveSplineCorrectionTermImpl(size_t num_base,
		XDataType const base_position[], size_t num_array,
		YDataType const base_data[], YDataType d2ydx2[],
		YDataType upper_triangular[]) {
	// This is a condition of natural cubic spline
	for (size_t i = 0; i < num_array; ++i) {
		d2ydx2[i] = 0.0;
		d2ydx2[(num_base - 1) * num_array + i] = 0.0;
		upper_triangular[i] = 0.0;
	}

	// Solve tridiagonal system.
	// Here tridiagonal matrix is decomposed to upper triangular matrix.
	// upper_tridiangular stores upper triangular elements, while
	// d2ydx2 stores right-hand-side vector. The diagonal
	// elements are normalized to 1.

	// x_base is ascending order
	XDataType a1 = base_position[1] - base_position[0];
	for (size_t i = 2; i < num_base; ++i) {
		XDataType const a2 = base_position[i] - base_position[i - 1];
		XDataType const b1 = 1.0 / (base_position[i] - base_position[i - 2]);
		for (size_t j = 0; j < num_array; ++j) {
			size_t i0 = num_array * i + j;
			size_t i1 = i0 - num_array;
			size_t i2 = i1 - num_array;
			d2ydx2[i1] = 3.0 * b1
					* ((base_data[i0] - base_data[i1]) / a2
							- (base_data[i1] - base_data[i2]) / a1
							- d2ydx2[i2] * 0.5 * a1);
			XDataType a3 = 1.0 / (1.0 - upper_triangular[i2] * 0.5 * a1 * b1);
			d2ydx2[i1] *= a3;
			upper_triangular[i1] = 0.5 * a2 * b1 * a3;
		}
		a1 = a2;
	}

	// Solve the system by backsubstitution and store solution to d2ydx2
	for (size_t k = num_base; k >= 3; --k) {
		for (size_t j = 0; j < num_array; ++j) {
			size_t const index = (k - 2) * num_array + j;
			d2ydx2[index] -= upper_triangular[index]
					* d2ydx2[index + num_array];
		}
	}
}

template<class XDataType, class YDataType>
struct SplineXWorkingData: public StorageAndAlignedPointer<YDataType> {
	static SplineXWorkingData<XDataType, YDataType> * Initialize(
			uint8_t polynomial_order, size_t num_base, size_t num_array,
			XDataType const base_position[], YDataType const base_data[]) {
		SplineXWorkingData<XDataType, YDataType> *work_data =
				new SplineXWorkingData<XDataType, YDataType>();
		AllocateAndAlign<YDataType>(num_base * num_array, work_data);
		work_data->DeriveSplineCorrectionTerm(num_base, base_position,
				num_array, base_data);
		return work_data;
	}
	void DeriveSplineCorrectionTerm(size_t num_base,
			XDataType const base_position[], size_t num_array,
			YDataType const base_data[]) {
		StorageAndAlignedPointer<YDataType> holder_for_u;
		AllocateAndAlign<YDataType>(num_base * num_array, &holder_for_u);
		YDataType *upper_triangular = holder_for_u.pointer;
		for (size_t i = 0; i < num_array; ++i) {
			size_t const index = i * num_base;
			YDataType const *base_data_work = &(base_data[index]);
			YDataType *d2ydx2 = &(this->pointer[index]);
			DeriveSplineCorrectionTermImpl(num_base, base_position, 1,
					base_data_work, d2ydx2, upper_triangular);
		}
	}
};

template<class XDataType, class YDataType>
struct SplineYWorkingData: public StorageAndAlignedPointer<YDataType> {
	static SplineYWorkingData<XDataType, YDataType> * Initialize(
			uint8_t polynomial_order, size_t num_base, size_t num_array,
			XDataType const base_position[], YDataType const base_data[]) {
		SplineYWorkingData<XDataType, YDataType> *work_data =
				new SplineYWorkingData<XDataType, YDataType>();
		AllocateAndAlign<YDataType>(num_base * num_array, work_data);
		work_data->DeriveSplineCorrectionTerm(num_base, base_position,
				num_array, base_data);
		return work_data;
	}
	void DeriveSplineCorrectionTerm(size_t num_base,
			XDataType const base_position[], size_t num_array,
			YDataType const base_data[]) {
		StorageAndAlignedPointer<YDataType> holder_for_u;
		AllocateAndAlign<YDataType>(num_base * num_array, &holder_for_u);
		YDataType *upper_triangular = holder_for_u.pointer;
		DeriveSplineCorrectionTermImpl<XDataType, YDataType>(num_base,
				base_position, num_array, base_data, this->pointer,
				upper_triangular);
	}
};

// Interface class for Interpolator
template<class InterpolatorImpl, class WorkData, class XDataType,
		class YDataType>
struct InterpolatorInterface {
	typedef WorkData WorkingData;
	static void Interpolate1D(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t num_location,
			size_t const location[], size_t offset,
			WorkingData const * const work_data) {
		for (size_t k = 1; k < num_location; ++k) {
			size_t const left_index = offset + k - 1;
			InterpolatorImpl::DoInterpolate(num_base, base_position, num_array,
					base_data, num_interpolated, interpolated_position,
					interpolated_data, location, k, left_index, work_data);
		}
	}
};

// TODO: documentation needs to be improved.
/**
 * Nearest interpolation engine classes
 *
 * nearest condition
 * - interpolated_position[i] == midpoint
 *     ---> nearest is left_value (left side)
 * - interpolated_position[i] < midpoint and ascending order
 *     ---> nearest is left_value (left side)
 * - interpolated_position[i] > midpoint and ascending order
 *     ---> nearest is right_value (right side)
 * - interpolated_position[i] < midpoint and descending order
 *     ---> nearest is right_value (right side)
 * - interpolated_position[i] > midpoint and descending order
 *     ---> nearest is left_value (left side)
 */
template<class XDataType, class YDataType>
struct NearestXInterpolatorImpl: public InterpolatorInterface<
		NearestXInterpolatorImpl<XDataType, YDataType>,
		NullWorkingData<XDataType, YDataType>, XDataType, YDataType> {
	typedef typename InterpolatorInterface<
			NearestXInterpolatorImpl<XDataType, YDataType>,
			NullWorkingData<XDataType, YDataType>, XDataType, YDataType>::WorkingData WorkingData;
	static void DoInterpolate(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t const location[], size_t k,
			size_t left_index, WorkingData const * const work_data) {
		XDataType const midpoint = 0.5
				* (base_position[left_index + 1] + base_position[left_index]);
		for (size_t j = 0; j < num_array; ++j) {
			YDataType const left_value = base_data[j * num_base + left_index];
			YDataType const right_value = base_data[j * num_base + left_index
					+ 1];
			YDataType *work = &interpolated_data[j * num_interpolated];
			for (size_t i = location[k - 1]; i < location[k]; ++i) {
				if ((interpolated_position[i] != midpoint)
						& ((interpolated_position[i] > midpoint))) {
					work[i] = right_value;
				} else {
					work[i] = left_value;
				}
			}
		}
	}
};

template<class XDataType, class YDataType>
struct NearestYInterpolatorImpl: public InterpolatorInterface<
		NearestYInterpolatorImpl<XDataType, YDataType>,
		NullWorkingData<XDataType, YDataType>, XDataType, YDataType> {
	typedef typename InterpolatorInterface<
			NearestYInterpolatorImpl<XDataType, YDataType>,
			NullWorkingData<XDataType, YDataType>, XDataType, YDataType>::WorkingData WorkingData;
	static void DoInterpolate(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t const location[], size_t k,
			size_t left_index, WorkingData const * const work_data) {
		XDataType const midpoint = 0.5
				* (base_position[left_index + 1] + base_position[left_index]);
		for (size_t i = location[k - 1]; i < location[k]; ++i) {
			size_t offset_index = 0;
			if ((interpolated_position[i] != midpoint)
					& ((interpolated_position[i] > midpoint))) {
				offset_index = 1;
			}
			YDataType *work = &interpolated_data[num_array * i];
			YDataType const * const nearest = &base_data[num_array
					* (left_index + offset_index)];
			for (size_t j = 0; j < num_array; ++j) {
				work[j] = nearest[j];
			}
		}
	}
};

// TODO: documentation must be added
template<class XDataType, class YDataType>
struct LinearXInterpolatorImpl: public InterpolatorInterface<
		LinearXInterpolatorImpl<XDataType, YDataType>,
		NullWorkingData<XDataType, YDataType>, XDataType, YDataType> {
	typedef typename InterpolatorInterface<
			LinearXInterpolatorImpl<XDataType, YDataType>,
			NullWorkingData<XDataType, YDataType>, XDataType, YDataType>::WorkingData WorkingData;
	static void DoInterpolate(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t const location[], size_t k,
			size_t left_index, WorkingData const * const work_data) {
		for (size_t j = 0; j < num_array; ++j) {
			size_t const offset_index_left = j * num_base + left_index;
			XDataType const dydx =
					static_cast<XDataType>(base_data[offset_index_left + 1]
							- base_data[offset_index_left])
							/ (base_position[left_index + 1]
									- base_position[left_index]);
			YDataType *y_work = &interpolated_data[j * num_interpolated];
			for (size_t i = location[k - 1]; i < location[k]; ++i) {
				y_work[i] = base_data[offset_index_left]
						+ static_cast<YDataType>(dydx
								* (interpolated_position[i]
										- base_position[left_index]));
			}
		}
	}
};

template<class XDataType, class YDataType>
struct LinearYInterpolatorImpl: public InterpolatorInterface<
		LinearYInterpolatorImpl<XDataType, YDataType>,
		NullWorkingData<XDataType, YDataType>, XDataType, YDataType> {
	typedef typename InterpolatorInterface<
			LinearYInterpolatorImpl<XDataType, YDataType>,
			NullWorkingData<XDataType, YDataType>, XDataType, YDataType>::WorkingData WorkingData;
	static void DoInterpolate(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t const location[], size_t k,
			size_t left_index, WorkingData const * const work_data) {
		for (size_t i = location[k - 1]; i < location[k]; ++i) {
			YDataType const fraction =
					static_cast<YDataType>((interpolated_position[i]
							- base_position[left_index])
							/ (base_position[left_index + 1]
									- base_position[left_index]));
			YDataType *work = &interpolated_data[num_array * i];
			YDataType const *left_value = &base_data[num_array * left_index];
			YDataType const *right_value = &base_data[num_array
					* (left_index + 1)];
			for (size_t j = 0; j < num_array; ++j) {
				work[j] = left_value[j]
						+ fraction * (right_value[j] - left_value[j]);
			}
		}
	}
};

// TODO: documentation must be added
template<class XDataType, class YDataType>
struct PolynomialXInterpolatorImpl: public InterpolatorInterface<
		PolynomialXInterpolatorImpl<XDataType, YDataType>,
		PolynomialXWorkingData<XDataType, YDataType>, XDataType, YDataType> {
	typedef typename InterpolatorInterface<
			PolynomialXInterpolatorImpl<XDataType, YDataType>,
			PolynomialXWorkingData<XDataType, YDataType>, XDataType, YDataType>::WorkingData WorkingData;
	static void DoInterpolate(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t const location[], size_t k,
			size_t left_index, WorkingData const * const work_data) {
		int const left_edge1 = left_index - work_data->polynomial_order / 2;
		size_t const left_edge2 = num_base - work_data->num_elements;
		size_t left_edge =
				static_cast<size_t>((left_edge1 > 0) ? left_edge1 : 0);
		left_edge = (left_edge > left_edge2) ? left_edge2 : left_edge;
		for (size_t j = 0; j < num_array; ++j) {
			for (size_t i = location[k - 1]; i < location[k]; ++i) {
				PerformNevilleAlgorithm(num_base, base_position, base_data,
						left_edge, j, interpolated_position[i],
						&interpolated_data[j * num_interpolated + i],
						work_data);
			}
		}
	}
private:
	static void PerformNevilleAlgorithm(size_t num_base,
			XDataType const base_position[], YDataType const base_data[],
			size_t left_index, size_t array_index,
			XDataType interpolated_position, YDataType *interpolated_data,
			WorkingData const * const work_data) {

		// working pointers
		XDataType const *x_ptr = &(base_position[left_index]);
		YDataType const *y_ptr =
				&(base_data[num_base * array_index + left_index]);

		XDataType *c = work_data->xholder[0].pointer;
		XDataType *d = work_data->xholder[1].pointer;

		for (size_t i = 0; i < work_data->num_elements; ++i) {
			c[i] = static_cast<XDataType>(y_ptr[i]);
			d[i] = c[i];
		}

		// Neville's algorithm
		XDataType work = c[0];
		for (size_t m = 1; m < work_data->num_elements; ++m) {
			// Evaluate Cm1, Cm2, Cm3, ... Cm[n-m] and Dm1, Dm2, Dm3, ... Dm[n-m].
			// Those are stored to c[0], c[1], ..., c[n-m-1] and d[0], d[1], ...,
			// d[n-m-1].
			for (size_t i = m; i < work_data->num_elements; ++i) {
				XDataType cd = c[i + 1 - m] - d[i - m];
				XDataType const dx = x_ptr[i - m] - x_ptr[i];
				assert(dx != 0);
				cd /= dx;
				c[i - m] = (x_ptr[i - m] - interpolated_position) * cd;
				d[i - m] = (x_ptr[i] - interpolated_position) * cd;
			}

			// In each step, c[0] holds Cm1 which is a correction between
			// P12...m and P12...[m+1]. Thus, the following repeated update
			// corresponds to the route P1 -> P12 -> P123 ->...-> P123...n.
			work += c[0];
		}

		*interpolated_data = static_cast<YDataType>(work);
	}
};

template<class XDataType, class YDataType>
struct PolynomialYInterpolatorImpl: public InterpolatorInterface<
		PolynomialYInterpolatorImpl<XDataType, YDataType>,
		PolynomialYWorkingData<XDataType, YDataType>, XDataType, YDataType> {
	typedef typename InterpolatorInterface<
			PolynomialYInterpolatorImpl<XDataType, YDataType>,
			PolynomialYWorkingData<XDataType, YDataType>, XDataType, YDataType>::WorkingData WorkingData;
	static void DoInterpolate(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t const location[], size_t k,
			size_t left_index, WorkingData const * const work_data) {
		int left_edge1 = left_index - work_data->polynomial_order / 2;
		size_t left_edge2 = num_base - work_data->num_elements;
		size_t left_edge =
				static_cast<size_t>((left_edge1 > 0) ? left_edge1 : 0);
		left_edge = (left_edge > left_edge2) ? left_edge2 : left_edge;
		for (size_t i = location[k - 1]; i < location[k]; ++i) {
			PerformNevilleAlgorithm(num_array, base_position, base_data,
					left_edge, interpolated_position[i],
					&interpolated_data[num_array * i], work_data);
		}
	}
private:
	static void PerformNevilleAlgorithm(size_t num_array,
			XDataType const base_position[], YDataType const base_data[],
			size_t left_index, XDataType interpolated_position,
			YDataType interpolated_data[],
			WorkingData const * const work_data) {
		// working pointers
		XDataType const *x_ptr = &(base_position[left_index]);
		XDataType *c = work_data->xholder[0].pointer;
		XDataType *d = work_data->xholder[1].pointer;
		XDataType *work = work_data->xholder[2].pointer;

		for (size_t i = 0; i < work_data->num_elements * num_array; ++i) {
			c[i] =
					static_cast<XDataType>(base_data[left_index * num_array + i]);
			d[i] = c[i];
		}

		// Neville's algorithm
		for (size_t i = 0; i < num_array; ++i) {
			work[i] = c[i];
		}
		for (size_t m = 1; m < work_data->num_elements; ++m) {
			// Evaluate Cm1, Cm2, Cm3, ... Cm[n-m] and Dm1, Dm2, Dm3, ... Dm[n-m].
			// Those are stored to c[0], c[1], ..., c[n-m-1] and d[0], d[1], ...,
			// d[n-m-1].
			for (size_t i = m; i < work_data->num_elements; ++i) {
				XDataType const dx = x_ptr[i - m] - x_ptr[i];
				assert(dx != 0);
				size_t const offset = (i - m) * num_array;
				for (size_t j = 0; j < num_array; ++j) {
					XDataType cd = (c[offset + num_array + j] - d[offset + j])
							/ dx;
					c[offset + j] = (x_ptr[i - m] - interpolated_position) * cd;
					d[offset + j] = (x_ptr[i] - interpolated_position) * cd;
				}
			}

			// In each step, c[0] holds Cm1 which is a correction between
			// P12...m and P12...[m+1]. Thus, the following repeated update
			// corresponds to the route P1 -> P12 -> P123 ->...-> P123...n.
			for (size_t i = 0; i < num_array; ++i) {
				work[i] += c[i];
			}
		}

		for (size_t i = 0; i < num_array; ++i) {
			interpolated_data[i] = static_cast<YDataType>(work[i]);
		}
	}
};

// TODO: documentation must be added
template<class XDataType, class YDataType>
struct SplineXInterpolatorImpl: public InterpolatorInterface<
		SplineXInterpolatorImpl<XDataType, YDataType>,
		SplineXWorkingData<XDataType, YDataType>, XDataType, YDataType> {
	typedef typename InterpolatorInterface<
			SplineXInterpolatorImpl<XDataType, YDataType>,
			SplineXWorkingData<XDataType, YDataType>, XDataType, YDataType>::WorkingData WorkingData;
	static void DoInterpolate(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t const location[], size_t k,
			size_t left_index, WorkingData const * const work_data) {
		YDataType const * const d2ydx2 = work_data->pointer;
		XDataType const dx = base_position[left_index + 1]
				- base_position[left_index];
		XDataType const dx_factor = dx * dx / 6.0;
		for (size_t j = 0; j < num_array; ++j) {
			size_t const offset_index_left = j * num_base + left_index;
			YDataType *work = &interpolated_data[j * num_interpolated];
			for (size_t i = location[k - 1]; i < location[k]; ++i) {
				XDataType const a = (base_position[left_index + 1]
						- interpolated_position[i]) / dx;
				XDataType const b = 1.0 - a;
				work[i] = static_cast<YDataType>(a
						* base_data[offset_index_left]
						+ b * base_data[offset_index_left + 1]
						+ ((a * a * a - a) * d2ydx2[offset_index_left]
								+ (b * b * b - b)
										* d2ydx2[offset_index_left + 1])
								* dx_factor);
			}
		}
	}
};

template<class XDataType, class YDataType>
struct SplineYInterpolatorImpl: public InterpolatorInterface<
		SplineYInterpolatorImpl<XDataType, YDataType>,
		SplineYWorkingData<XDataType, YDataType>, XDataType, YDataType> {
	typedef typename InterpolatorInterface<
			SplineYInterpolatorImpl<XDataType, YDataType>,
			SplineYWorkingData<XDataType, YDataType>, XDataType, YDataType>::WorkingData WorkingData;
	static void DoInterpolate(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t const location[], size_t k,
			size_t left_index, WorkingData const * const work_data) {
		YDataType const * const d2ydx2 = work_data->pointer;
		XDataType const dx = base_position[left_index + 1]
				- base_position[left_index];
		YDataType const * const left_value = &base_data[left_index * num_array];
		YDataType const * const right_value = &base_data[(left_index + 1)
				* num_array];
		YDataType const * const d2ydx2_left = &d2ydx2[left_index * num_array];
		YDataType const * const d2ydx2_right = &d2ydx2[(left_index + 1)
				* num_array];
		for (size_t i = location[k - 1]; i < location[k]; ++i) {
			XDataType const a = (base_position[left_index + 1]
					- interpolated_position[i]) / dx;
			XDataType const b = 1.0 - a;
			YDataType *work = &interpolated_data[i * num_array];
			XDataType const aaa = (a * a * a - a) * dx * dx / 6.0;
			XDataType const bbb = (b * b * b - b) * dx * dx / 6.0;
			for (size_t j = 0; j < num_array; ++j) {
				work[j] = static_cast<YDataType>(a * left_value[j]
						+ b * right_value[j] + aaa * d2ydx2_left[j]
						+ bbb * d2ydx2_right[j]);
			}
		}
	}
};

template<class XDataType, class YDataType>
struct XInterpolatorHelper {
//	typedef NearestXInterpolatorImpl<XDataType, YDataType> NearestInterpolator;
//	typedef LinearXInterpolatorImpl<XDataType, YDataType> LinearInterpolator;
//	typedef PolynomialXInterpolatorImpl<XDataType, YDataType> PolynomialInterpolator;
//	typedef SplineXInterpolatorImpl<XDataType, YDataType> SplineInterpolator;

	static size_t FillDataAsAscending(size_t num_base, size_t num_array,
			size_t iarray, XDataType const position[], YDataType const data[],
			bool const mask[], bool is_ascending, XDataType x[],
			YDataType y[]) {
		assert(iarray < num_array);
		size_t n = 0;
		for (size_t i = 0; i < num_base; ++i) {
			size_t index = i + iarray * num_base;
			XDataType p;
			YDataType d;
			bool m;
			if (is_ascending) {
				p = position[i];
				d = data[index];
				m = mask[index];
			} else {
				p = position[num_base - 1 - i];
				d = data[num_base - 1 - i + iarray * num_base];
				m = mask[num_base - 1 - i + iarray * num_base];
			}
			if (m) {
				x[n] = p;
				y[n] = d;
				++n;
			}
		}
		return n;
	}

	static void MaskAll(size_t num_interpolated, size_t num_array,
			size_t iarray, bool mask[]) {
		assert(iarray < num_array);
		for (size_t i = 0; i < num_interpolated; ++i) {
			mask[num_interpolated * iarray + i] = false;
		}
	}

	static void SubstituteSingleBaseDataPerArray(size_t num_interpolated,
			size_t num_array, size_t iarray, YDataType value,
			YDataType data[]) {
		assert(iarray < num_array);
		for (size_t i = 0; i < num_interpolated; ++i) {
			data[num_interpolated * iarray + i] = value;
		}
	}

	static void FillResult(size_t num_interpolated, size_t num_array,
			size_t iarray, bool is_ascending, YDataType y1[],
			YDataType data[]) {
		assert(iarray < num_array);
		if (is_ascending) {
			for (size_t i = 0; i < num_interpolated; ++i) {
				data[i + iarray * num_interpolated] = y1[i];
			}
		} else {
			for (size_t i = 0; i < num_interpolated; ++i) {
				data[num_interpolated - 1 - i + num_interpolated * iarray] =
						y1[i];
			}
		}

	}

	template<class DataType>
	static inline void Reorder(size_t num_x, size_t num_y,
			DataType const input_data[], DataType output_data[]) {
		for (size_t i = 0; i < num_y; ++i) {
			size_t const start_position = num_x * i;
			size_t const end_position = start_position + num_x;
			for (size_t j = start_position; j < end_position; ++j) {
				output_data[j] = input_data[end_position - (j - start_position)
						- 1];
			}
		}
	}
	static inline void SubstituteSingleBaseData(size_t num_base,
			size_t num_array, size_t num_interpolated,
			YDataType const base_data[], YDataType interpolated_data[]) {
		Substitute(0, num_interpolated, num_array, 0, num_interpolated,
				num_base, base_data, interpolated_data);
	}
	static inline void SubstituteLeftMostData(size_t location, size_t num_base,
			size_t num_array, size_t num_interpolated,
			YDataType const base_data[], YDataType interpolated_data[]) {
		Substitute(0, location, num_array, 0, num_interpolated, num_base,
				base_data, interpolated_data);
	}
	static inline void SubstituteRightMostData(size_t location, size_t num_base,
			size_t num_array, size_t num_interpolated,
			YDataType const base_data[], YDataType interpolated_data[]) {
		Substitute(location, num_interpolated, num_array, num_base - 1,
				num_interpolated, num_base, base_data, interpolated_data);
	}
	static inline void SwapResult(size_t num_array, size_t num_interpolated,
			YDataType interpolated_data[]) {
		size_t const midpoint = num_interpolated / 2;
		size_t const right_edge = num_interpolated - 1;
		for (size_t j = 0; j < num_array; ++j) {
			YDataType *work = &interpolated_data[j * num_interpolated];
			for (size_t i = 0; i < midpoint; ++i) {
				std::swap<YDataType>(work[i], work[right_edge - i]);
			}
		}
	}
private:
	static inline void Substitute(size_t start, size_t end, size_t num,
			size_t offset, size_t step_in, size_t step_out,
			YDataType const in[], YDataType out[]) {
		for (size_t j = 0; j < num; ++j) {
			for (size_t i = start; i < end; ++i) {
				out[j * step_in + i] = in[j * step_out + offset];
			}
		}
	}
};

template<class XDataType, class YDataType>
struct YInterpolatorHelper {
//	typedef NearestYInterpolatorImpl<XDataType, YDataType> NearestInterpolator;
//	typedef LinearYInterpolatorImpl<XDataType, YDataType> LinearInterpolator;
//	typedef PolynomialYInterpolatorImpl<XDataType, YDataType> PolynomialInterpolator;
//	typedef SplineYInterpolatorImpl<XDataType, YDataType> SplineInterpolator;

	static size_t FillDataAsAscending(size_t num_base, size_t num_array,
			size_t iarray, XDataType const position[], YDataType const data[],
			bool const mask[], bool is_ascending, XDataType x[],
			YDataType y[]) {
		assert(iarray < num_array);
		size_t n = 0;
		for (size_t i = 0; i < num_base; ++i) {
			size_t const index = num_array * i + iarray;
			XDataType p;
			YDataType d;
			bool m;
			if (is_ascending) {
				p = position[i];
				d = data[index];
				m = mask[index];
			} else {
				p = position[num_base - 1 - i];
				d = data[num_array * (num_base - 1 - i) + iarray];
				m = mask[num_array * (num_base - 1 - i) + iarray];
			}
			if (m) {
				x[n] = p;
				y[n] = d;
				++n;
			}
		}
		return n;
	}

	static void MaskAll(size_t num_interpolated, size_t num_array,
			size_t iarray, bool mask[]) {
		assert(iarray < num_array);
		for (size_t i = 0; i < num_interpolated; ++i) {
			mask[iarray + num_array * i] = false;
		}
	}

	static void SubstituteSingleBaseDataPerArray(size_t num_interpolated,
			size_t num_array, size_t iarray, YDataType value,
			YDataType data[]) {
		assert(iarray < num_array);
		for (size_t i = 0; i < num_interpolated; ++i) {
			data[iarray + num_array * i] = value;
		}
	}

	static void FillResult(size_t num_interpolated, size_t num_array,
			size_t iarray, bool is_ascending, YDataType y1[],
			YDataType data[]) {
		assert(iarray < num_array);
		if (is_ascending) {
			for (size_t i = 0; i < num_interpolated; ++i) {
				data[iarray + i * num_array] = y1[i];
			}
		} else {
			for (size_t i = 0; i < num_interpolated; ++i) {
				data[num_array * (num_interpolated - 1 - i) + iarray] = y1[i];
			}
		}
	}

	template<class DataType>
	static inline void Reorder(size_t num_x, size_t num_y,
			DataType const input_data[], DataType output_data[]) {
		for (size_t i = 0; i < num_x; ++i) {
			DataType *out_storage = &output_data[i * num_y];
			DataType const *in_storage = &input_data[(num_x - 1 - i) * num_y];
			for (size_t j = 0; j < num_y; ++j) {
				out_storage[j] = in_storage[j];
			}
		}
	}
	static inline void SubstituteSingleBaseData(size_t num_base,
			size_t num_array, size_t num_interpolated,
			YDataType const base_data[], YDataType interpolated_data[]) {
		Substitute(0, num_interpolated, num_array, 0, base_data,
				interpolated_data);
	}
	static inline void SubstituteLeftMostData(size_t location, size_t num_base,
			size_t num_array, size_t num_interpolated,
			YDataType const base_data[], YDataType interpolated_data[]) {
		Substitute(0, location, num_array, 0, base_data, interpolated_data);
	}
	static inline void SubstituteRightMostData(size_t location, size_t num_base,
			size_t num_array, size_t num_interpolated,
			YDataType const base_data[], YDataType interpolated_data[]) {
		Substitute(location, num_interpolated, num_array, num_base - 1,
				base_data, interpolated_data);
	}
	static inline void SwapResult(size_t num_array, size_t num_interpolated,
			YDataType interpolated_data[]) {
		size_t midpoint = num_interpolated / 2;
		for (size_t i = 0; i < midpoint; ++i) {
			YDataType *a = &interpolated_data[i * num_array];
			YDataType *b = &interpolated_data[(num_interpolated - 1 - i)
					* num_array];
			for (size_t j = 0; j < num_array; ++j) {
				std::swap<YDataType>(a[j], b[j]);
			}
		}
	}
private:
	static inline void Substitute(size_t start, size_t end, size_t num,
			size_t offset, YDataType const in[], YDataType out[]) {
		for (size_t i = start; i < end; ++i) {
			for (size_t j = 0; j < num; ++j) {
				out[i * num + j] = in[offset * num + j];
			}
		}
	}
};

template<class Interpolator, class Helper, class XDataType, class YDataType>
void Interpolate1D(uint8_t polynomial_order, size_t num_base,
		XDataType const base_position[/*num_base*/], size_t num_array,
		YDataType const base_data[/*num_base*num_array*/],
		bool const base_mask[/*num_base*num_array*/], size_t num_interpolated,
		XDataType const interpolated_position[/*num_interpolated*/],
		YDataType interpolated_data[/*num_interpolated*num_array*/],
		bool interpolated_mask[/*num_interpolated*num_array*/]) {
	typedef typename Interpolator::WorkingData WorkingData;
	assert(num_base > 0);
	assert(num_array > 0);
	assert(num_interpolated > 0);
	assert(base_position != nullptr);
	assert(base_data != nullptr);
	assert(interpolated_position != nullptr);
	assert(interpolated_data != nullptr);
	assert(base_mask != nullptr);
	assert(interpolated_mask != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(base_position));
	assert(LIBSAKURA_SYMBOL(IsAligned)(base_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(interpolated_position));
	assert(LIBSAKURA_SYMBOL(IsAligned)(interpolated_data));

	// basic check
	bool const is_base_ascending = (base_position[num_base - 1]
			> base_position[0]);
	bool const is_interp_ascending =
			(interpolated_position[num_interpolated - 1]
					> interpolated_position[0]);

	// Initialize interpolated_mask to true
	for (size_t i = 0; i < num_interpolated * num_array; ++i) {
		interpolated_mask[i] = true;
	}

	// Working arrays
	StorageAndAlignedPointer<XDataType> x0_storage, x1_storage;
	StorageAndAlignedPointer<YDataType> y0_storage, y1_storage;
	AllocateAndAlign(num_base, &x0_storage);
	AllocateAndAlign(num_base, &y0_storage);
	XDataType *x0 = x0_storage.pointer;
	YDataType *y0 = y0_storage.pointer;
	XDataType const *x1 = interpolated_position;
	if (!is_interp_ascending) {
		AllocateAndAlign(num_interpolated, &x1_storage);
		x1 = x1_storage.pointer;
		XDataType *_x = const_cast<XDataType *>(x1);
		for (size_t i = 0; i < num_interpolated; ++i) {
			_x[i] = interpolated_position[num_interpolated - 1 - i];
		}
	}
	AllocateAndAlign(num_interpolated, &y1_storage);
	YDataType *y1 = y1_storage.pointer;

	for (size_t i = 0; i < num_array; ++i) {
		size_t const n = Helper::FillDataAsAscending(num_base, num_array, i,
				base_position, base_data, base_mask, is_base_ascending, x0, y0);
		if (n == 0) {
			// cannot perform interpolation, mask all data
			Helper::MaskAll(num_interpolated, num_array, i, interpolated_mask);
		} else if (n == 1) {
			// no need to interpolate, just substitute single base data
			// to all elements in interpolated data, keep input mask
			Helper::SubstituteSingleBaseDataPerArray(num_interpolated,
					num_array, i, y0[0], interpolated_data);
		} else {
			// perform interpolation, keep input mask

			// Perform 1-dimensional interpolation
			// Any preparation for interpolation should be done here
			std::unique_ptr<WorkingData> wdata_storage(
					WorkingData::Initialize(polynomial_order, n, 1, x0, y0));
			WorkingData const * const work_data = wdata_storage.get();

			// Locate each element in base_position against interpolated_position
			StorageAndAlignedPointer<size_t> size_t_holder;
			AllocateAndAlign<size_t>(n, &size_t_holder);
			size_t *location_base = size_t_holder.pointer;
			size_t const num_location_base = Locate<XDataType>(num_interpolated, n,
					x1, x0, location_base);

			// Outside of base_position[0]
			Helper::SubstituteLeftMostData(location_base[0], n, 1,
					num_interpolated, y0, y1);

			// Between base_position[0] and base_position[num_base-1]
			size_t offset = 0;
			if (x0[0] < x1[0]) {
				for (size_t i = 1; i < n; ++i) {
					if (x0[offset + 1] < x1[0]) {
						offset++;
					} else {
						break;
					}
				}
			}
			Interpolator::Interpolate1D(n, x0, 1, y0, num_interpolated, x1, y1,
					num_location_base, location_base, offset, work_data);

			// Outside of base_position[num_base-1]
			Helper::SubstituteRightMostData(
					location_base[num_location_base - 1], n, 1,
					num_interpolated, y0, y1);

			Helper::FillResult(num_interpolated, num_array, i,
					is_interp_ascending, y1, interpolated_data);
		}
	}
}

template<class InterpolatorHelper, class XDataType, class YDataType>
void ExecuteInterpolate1D(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_base,
		XDataType const base_position[], size_t num_array,
		YDataType const base_data[], bool const base_mask[],
		size_t num_interpolated, XDataType const interpolated_position[],
		YDataType interpolated_data[], bool interpolated_mask[]) {
	typedef NearestXInterpolatorImpl<XDataType, YDataType> NearestInterpolator;
	typedef LinearXInterpolatorImpl<XDataType, YDataType> LinearInterpolator;
	typedef PolynomialXInterpolatorImpl<XDataType, YDataType> PolynomialInterpolator;
	typedef SplineXInterpolatorImpl<XDataType, YDataType> SplineInterpolator;
	typedef void (*Interpolate1DFunc)(uint8_t, size_t, XDataType const *,
			size_t, YDataType const *, bool const *, size_t, XDataType const *,
			YDataType *, bool *);
	Interpolate1DFunc func = [](uint8_t, size_t, XDataType const *,
			size_t, YDataType const *, bool const *, size_t,
			XDataType const *, YDataType *, bool *) {
		throw LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	};
	switch (interpolation_method) {
	case LIBSAKURA_SYMBOL(InterpolationMethod_kNearest):
		func = Interpolate1D<NearestInterpolator, InterpolatorHelper, XDataType,
				YDataType>;
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kLinear):
		func = Interpolate1D<LinearInterpolator, InterpolatorHelper, XDataType,
				YDataType>;
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kPolynomial):
		if (polynomial_order == 0) {
			// This is special case: 0-th polynomial interpolation
			// acts like nearest interpolation
			func = Interpolate1D<NearestInterpolator, InterpolatorHelper,
					XDataType, YDataType>;
		} else {
			func = Interpolate1D<PolynomialInterpolator, InterpolatorHelper,
					XDataType, YDataType>;
		}
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kSpline):
		func = Interpolate1D<SplineInterpolator, InterpolatorHelper, XDataType,
				YDataType>;
		break;
	default:
		// invalid interpolation method type
		break;
	}
	(*func)(polynomial_order, num_base, base_position, num_array, base_data,
			base_mask, num_interpolated, interpolated_position,
			interpolated_data, interpolated_mask);
}

// basic check of arguments
bool IsValidArguments(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_interpolation_axis,
		double const base[], size_t num_array, float const data_base[],
		bool const mask_base[], size_t num_interpolated,
		double const interpolated[], float const data_interpolated[],
		bool const mask_interpolated[], LIBSAKURA_SYMBOL(
				Status)*status) {
	// check interpolation_method
	if ((interpolation_method < 0)
			&& (interpolation_method
					>= LIBSAKURA_SYMBOL(InterpolationMethod_kNumElements))) {
		LOG4CXX_ERROR(logger, "Invalid interpolation method");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		return false;
	}

	// num_base must be non-zero
	if (num_interpolation_axis == 0) {
		LOG4CXX_ERROR(logger, "num_base must be >0");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		return false;
	}

	// no interpolation will be done
	if (num_interpolated == 0 || num_array == 0) {
		// Nothing to do
		LOG4CXX_INFO(logger,
				"Nothing has been done since num_interpolated is 0");
		*status = LIBSAKURA_SYMBOL(Status_kOK);
		return false;
	}

	// input arrays are not aligned
	if (!LIBSAKURA_SYMBOL(IsAligned)(base)
			|| !LIBSAKURA_SYMBOL(IsAligned)(data_base) || !LIBSAKURA_SYMBOL(
					IsAligned)(interpolated)
			|| !LIBSAKURA_SYMBOL(IsAligned)(data_interpolated)) {
		LOG4CXX_ERROR(logger, "input arrays are not aligned");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		return false;
	}

	// input arrays are null
	if (base == nullptr || data_base == nullptr || interpolated == nullptr
			|| data_interpolated == nullptr || mask_base == nullptr
			|| mask_interpolated == nullptr) {
		LOG4CXX_ERROR(logger, "input arrays are null");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		return false;
	}
	return true;
}

template<typename Func>
LIBSAKURA_SYMBOL(Status) DoInterpolate(Func func,
LIBSAKURA_SYMBOL(
		InterpolationMethod) interpolation_method, uint8_t polynomial_order,
		size_t num_base, double const base[/*num_base*/], size_t num_array,
		float const data_base[/*num_base*num_array*/],
		bool const mask_base[/*num_base*num_array*/], size_t num_interpolated,
		double const interpolated[/*num_x_interpolated*/],
		float data_interpolated[/*num_x_interpolated*num_y*/],
		bool mask_interpolated[/*num_x_interpolated*num_y*/]) {
	// check arguments
	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Status_kOK);
	if (!IsValidArguments(interpolation_method, polynomial_order, num_base, base,
			num_array, data_base, mask_base, num_interpolated, interpolated,
			data_interpolated, mask_interpolated, &status)) {
		return status;
	}

	try {
		func(interpolation_method, polynomial_order, num_base, base, num_array,
				data_base, mask_base, num_interpolated, interpolated,
				data_interpolated, mask_interpolated);
	} catch (const std::bad_alloc &e) {
		// failed to allocate memory
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const LIBSAKURA_SYMBOL(Status)&LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		// failed to allocate memory
		LOG4CXX_ERROR(logger, "Invalid interpolation type.");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	} catch (...) {
		// any exception is thrown during interpolation
		assert(false);
		LOG4CXX_ERROR(logger, "Aborted due to unknown error");
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return status;
}

} /* anonymous namespace */

/**
 * InterpolateXAxisFloat performs 1D interpolation along column based on x_base
 * and data_base. data_base is a serial array of column-major matrix data.
 * Its memory layout is assumed to be:
 *     data_base[0]     = data[0][0]
 *     data_base[1]     = data[1][0]
 *     data_base[2]     = data[2][0]
 *     ...
 *     data_base[N]     = data[N][0]
 *     data_base[N+1]   = data[1][0]
 *     ...
 *     data_base[N*M-1] = data[N][M]
 * where N and M correspond to num_x_base and num_y respectively.
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InterpolateXAxisFloat)(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_x_base,
		double const x_base[/*num_x_base*/], size_t num_y,
		float const data_base[/*num_x_base*num_y*/],
		bool const mask_base[/*num_x_base*num_y*/], size_t num_x_interpolated,
		double const x_interpolated[/*num_x_interpolated*/],
		float data_interpolated[/*num_x_interpolated*num_y*/],
		bool mask_interpolated[/*num_x_interpolated*num_y*/]) noexcept {

	return DoInterpolate(
			ExecuteInterpolate1D<XInterpolatorHelper<double, float>, double,
					float>, interpolation_method, polynomial_order, num_x_base,
			x_base, num_y, data_base, mask_base, num_x_interpolated,
			x_interpolated, data_interpolated, mask_interpolated);

}

/**
 * Interpolate1DYAxisFloat performs 1D interpolation along column based on y_base
 * and data_base. data_base is a serial array of row-major matrix data.
 * Its memory layout is assumed to be:
 *     data_base[0]     = data[0][0]
 *     data_base[1]     = data[0][1]
 *     data_base[2]     = data[0][2]
 *     ...
 *     data_base[M]     = data[0][M]
 *     data_base[M+1]   = data[1][0]
 *     ...
 *     data_base[M*N-1] = data[N][M]
 * where N and M correspond to num_y_base and num_x respectively.
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InterpolateYAxisFloat)(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_y_base,
		double const y_base[/*num_y_base*/], size_t num_x,
		float const data_base[/*num_y_base*num_x*/],
		bool const mask_base[/*num_y_base*num_x*/], size_t num_y_interpolated,
		double const y_interpolated[/*num_y_interpolated*/],
		float data_interpolated[/*num_y_interpolated*num_x*/],
		bool mask_interpolated[/*num_y_interpolated*num_x*/]) noexcept {

	return DoInterpolate(
			ExecuteInterpolate1D<YInterpolatorHelper<double, float>, double,
					float>, interpolation_method, polynomial_order, num_y_base,
			y_base, num_x, data_base, mask_base, num_y_interpolated,
			y_interpolated, data_interpolated, mask_interpolated);

}
