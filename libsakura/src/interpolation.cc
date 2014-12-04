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

#include "locator.tcc"

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
inline void AllocateAndAlign(size_t num_array,
		StorageAndAlignedPointer<DataType> *holder) {
	holder->storage = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
			num_array * sizeof(DataType), &(holder->pointer));
}

template<class Reorderer, class XDataType, class YDataType>
inline void GetAscendingArray(size_t num_base,
		XDataType const base_array[/*num_base*/], size_t num_array,
		YDataType const unordered_array[/*num_base*num_array*/],
		StorageAndAlignedPointer<YDataType> *storage) {
	if (base_array[0] < base_array[num_base - 1]) {
		storage->pointer = const_cast<YDataType *>(unordered_array);
	} else {
		AllocateAndAlign(num_base * num_array, storage);
		Reorderer::Reorder(num_base, num_array, unordered_array,
				storage->pointer);
	}
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
		polynomial_order_ =
				(order + 1u >= num_base) ?
						static_cast<uint8_t>(num_base - 1) : order;
		num_elements_ = polynomial_order_ + 1;
		xholder_.resize(array_size);
		for (size_t i = 0; i < array_size; ++i) {
			AllocateAndAlign<DataType>(num_elements_ * num_array,
					&(xholder_[i]));
		}
	}
	uint8_t polynomial_order_;
	size_t num_elements_;
	std::vector<StorageAndAlignedPointer<DataType> > xholder_;
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
		XDataType a2 = base_position[i] - base_position[i - 1];
		XDataType b1 = 1.0 / (base_position[i] - base_position[i - 2]);
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

	// Solve the system by backsubstitution and store solution to d2ydx2_
	for (size_t k = num_base; k >= 3; --k) {
		for (size_t j = 0; j < num_array; ++j) {
			size_t index = (k - 2) * num_array + j;
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
			size_t index = i * num_base;
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
	static WorkingData * PrepareForInterpolation(uint8_t polynomial_order,
			size_t num_base, size_t num_array, XDataType const base_position[],
			YDataType const base_data[]) {
		return WorkingData::Initialize(polynomial_order, num_base, num_array,
				base_position, base_data);
	}
	static void Interpolate1D(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t num_location,
			size_t const location[], size_t offset, WorkingData *work_data) {
		for (size_t k = 1; k < num_location; ++k) {
			size_t left_index = offset + k - 1;
			XDataType middle_point =
					0.5
							* (base_position[left_index + 1]
									+ base_position[left_index]);
			InterpolatorImpl::DoInterpolate(num_base, base_position, num_array,
					base_data, num_interpolated, interpolated_position,
					interpolated_data, num_location, location, offset, k,
					left_index, middle_point, work_data);
		}
	}
};

// TODO: documentation needs to be improved.
/**
 * Nearest interpolation engine classes
 *
 * nearest condition
 * - interpolated_position[i] == middle_point
 *     ---> nearest is y_left (left side)
 * - interpolated_position[i] < middle_point and ascending order
 *     ---> nearest is y_left (left side)
 * - interpolated_position[i] > middle_point and ascending order
 *     ---> nearest is y_right (right side)
 * - interpolated_position[i] < middle_point and descending order
 *     ---> nearest is y_right (right side)
 * - interpolated_position[i] > middle_point and descending order
 *     ---> nearest is y_left (left side)
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
			YDataType interpolated_data[], size_t num_location,
			size_t const location[], size_t offset, size_t k, size_t left_index,
			XDataType middle_point, WorkingData *work_data) {
		for (size_t j = 0; j < num_array; ++j) {
			YDataType left_value = base_data[j * num_base + left_index];
			YDataType right_value = base_data[j * num_base + left_index + 1];
			YDataType *work = &interpolated_data[j * num_interpolated];
			for (size_t i = location[k - 1]; i < location[k]; ++i) {
				if ((interpolated_position[i] != middle_point)
						& ((interpolated_position[i] > middle_point))) {
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
			YDataType interpolated_data[], size_t num_location,
			size_t const location[], size_t offset, size_t k, size_t left_index,
			XDataType middle_point, WorkingData *work_data) {
		for (size_t i = location[k - 1]; i < location[k]; ++i) {
			size_t offset_index = 0;
			if ((interpolated_position[i] != middle_point)
					& ((interpolated_position[i] > middle_point))) {
				offset_index = 1;
			}
			YDataType *work = &interpolated_data[num_array * i];
			YDataType const *nearest = &base_data[num_array
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
			YDataType interpolated_data[], size_t num_location,
			size_t const location[], size_t offset, size_t k, size_t left_index,
			XDataType middle_point, WorkingData *work_data) {
		for (size_t j = 0; j < num_array; ++j) {
			size_t offset_index_left = j * num_base + left_index;
			XDataType dydx =
					static_cast<XDataType>(base_data[offset_index_left + 1]
							- base_data[offset_index_left])
							/ (base_position[left_index + 1]
									- base_position[left_index]);
			YDataType base_term = base_data[offset_index_left];
			YDataType *y_work = &interpolated_data[j * num_interpolated];
			for (size_t i = location[k - 1]; i < location[k]; ++i) {
				y_work[i] = base_term
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
			YDataType interpolated_data[], size_t num_location,
			size_t const location[], size_t offset, size_t k, size_t left_index,
			XDataType middle_point, WorkingData *work_data) {
		for (size_t i = location[k - 1]; i < location[k]; ++i) {
			YDataType fraction =
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
			YDataType interpolated_data[], size_t num_location,
			size_t const location[], size_t offset, size_t k, size_t left_index,
			XDataType middle_point, WorkingData *work_data) {
		int left_edge1 = left_index - work_data->polynomial_order_ / 2;
		size_t left_edge2 = num_base - work_data->num_elements_;
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
			WorkingData *work_data) {

		// working pointers
		XDataType const *x_ptr = &(base_position[left_index]);
		YDataType const *y_ptr =
				&(base_data[num_base * array_index + left_index]);

		XDataType *c = work_data->xholder_[0].pointer;
		XDataType *d = work_data->xholder_[1].pointer;

		for (size_t i = 0; i < work_data->num_elements_; ++i) {
			c[i] = static_cast<XDataType>(y_ptr[i]);
			d[i] = c[i];
		}

		// Neville's algorithm
		XDataType work = c[0];
		for (size_t m = 1; m < work_data->num_elements_; ++m) {
			// Evaluate Cm1, Cm2, Cm3, ... Cm[n-m] and Dm1, Dm2, Dm3, ... Dm[n-m].
			// Those are stored to c[0], c[1], ..., c[n-m-1] and d[0], d[1], ...,
			// d[n-m-1].
			for (size_t i = m; i < work_data->num_elements_; ++i) {
				XDataType cd = c[i + 1 - m] - d[i - m];
				XDataType dx = x_ptr[i - m] - x_ptr[i];
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
			YDataType interpolated_data[], size_t num_location,
			size_t const location[], size_t offset, size_t k, size_t left_index,
			XDataType middle_point, WorkingData *work_data) {
		int left_edge1 = left_index - work_data->polynomial_order_ / 2;
		size_t left_edge2 = num_base - work_data->num_elements_;
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
			YDataType interpolated_data[], WorkingData *work_data) {
		// working pointers
		XDataType const *x_ptr = &(base_position[left_index]);
		XDataType *c = work_data->xholder_[0].pointer;
		XDataType *d = work_data->xholder_[1].pointer;
		XDataType *work = work_data->xholder_[2].pointer;

		size_t start = left_index * num_array;
		size_t num_elements = work_data->num_elements_ * num_array;
		for (size_t i = 0; i < num_elements; ++i) {
			c[i] = static_cast<XDataType>(base_data[start + i]);
			d[i] = c[i];
		}

		// Neville's algorithm
		for (size_t i = 0; i < num_array; ++i) {
			work[i] = c[i];
		}
		for (size_t m = 1; m < work_data->num_elements_; ++m) {
			// Evaluate Cm1, Cm2, Cm3, ... Cm[n-m] and Dm1, Dm2, Dm3, ... Dm[n-m].
			// Those are stored to c[0], c[1], ..., c[n-m-1] and d[0], d[1], ...,
			// d[n-m-1].
			for (size_t i = m; i < work_data->num_elements_; ++i) {
				XDataType dx = x_ptr[i - m] - x_ptr[i];
				assert(dx != 0);
				size_t offset = (i - m) * num_array;
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
			YDataType interpolated_data[], size_t num_location,
			size_t const location[], size_t offset, size_t k, size_t left_index,
			XDataType middle_point, WorkingData *work_data) {
		YDataType *d2ydx2 = work_data->pointer;
		XDataType dx = base_position[left_index + 1]
				- base_position[left_index];
		XDataType dx_factor = dx * dx / 6.0;
		for (size_t j = 0; j < num_array; ++j) {
			size_t offset_index_left = j * num_base + left_index;
			YDataType *work = &interpolated_data[j * num_interpolated];
			for (size_t i = location[k - 1]; i < location[k]; ++i) {
				XDataType a = (base_position[left_index + 1]
						- interpolated_position[i]) / dx;
				XDataType b = 1.0 - a;
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
			YDataType interpolated_data[], size_t num_location,
			size_t const location[], size_t offset, size_t k, size_t left_index,
			XDataType middle_point, WorkingData *work_data) {
		YDataType *d2ydx2 = work_data->pointer;
		XDataType dx = base_position[left_index + 1]
				- base_position[left_index];
		YDataType const *left_value = &base_data[left_index * num_array];
		YDataType const *right_value = &base_data[(left_index + 1) * num_array];
		YDataType const *d2ydx2_left = &d2ydx2[left_index * num_array];
		YDataType const *d2ydx2_right = &d2ydx2[(left_index + 1) * num_array];
		for (size_t i = location[k - 1]; i < location[k]; ++i) {
			XDataType a = (base_position[left_index + 1]
					- interpolated_position[i]) / dx;
			XDataType b = 1.0 - a;
			YDataType *work = &interpolated_data[i * num_array];
			XDataType aaa = (a * a * a - a) * dx * dx / 6.0;
			XDataType bbb = (b * b * b - b) * dx * dx / 6.0;
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
	typedef NearestXInterpolatorImpl<XDataType, YDataType> NearestInterpolator;
	typedef LinearXInterpolatorImpl<XDataType, YDataType> LinearInterpolator;
	typedef PolynomialXInterpolatorImpl<XDataType, YDataType> PolynomialInterpolator;
	typedef SplineXInterpolatorImpl<XDataType, YDataType> SplineInterpolator;

	template<class DataType>
	static inline void Reorder(size_t num_x, size_t num_y,
			DataType const input_data[], DataType output_data[]) {
		for (size_t i = 0; i < num_y; ++i) {
			size_t start_position = num_x * i;
			size_t end_position = start_position + num_x;
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
		size_t middle_point = num_interpolated / 2;
		size_t right_edge = num_interpolated - 1;
		for (size_t j = 0; j < num_array; ++j) {
			YDataType *work = &interpolated_data[j * num_interpolated];
			for (size_t i = 0; i < middle_point; ++i) {
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
	typedef NearestYInterpolatorImpl<XDataType, YDataType> NearestInterpolator;
	typedef LinearYInterpolatorImpl<XDataType, YDataType> LinearInterpolator;
	typedef PolynomialYInterpolatorImpl<XDataType, YDataType> PolynomialInterpolator;
	typedef SplineYInterpolatorImpl<XDataType, YDataType> SplineInterpolator;

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
		size_t middle_point = num_interpolated / 2;
		for (size_t i = 0; i < middle_point; ++i) {
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
		XDataType const base_position[/*num_x_base*/], size_t num_array,
		YDataType const base_data[/*num_x_base*num_y*/],
		size_t num_interpolated,
		XDataType const interpolated_position[/*num_x_interpolated*/],
		YDataType interpolated_data[/*num_x_interpolated*num_y*/]) {
	assert(num_base > 0);
	assert(num_array > 0);
	assert(num_interpolated > 0);
	assert(base_position != nullptr);
	assert(base_data != nullptr);
	assert(interpolated_position != nullptr);
	assert(interpolated_data != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(base_position));
	assert(LIBSAKURA_SYMBOL(IsAligned)(base_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(interpolated_position));
	assert(LIBSAKURA_SYMBOL(IsAligned)(interpolated_data));

	if (num_base == 1) {
		// No need to interpolate, just substitute base_data
		// to all elements in y_interpolated
		Helper::SubstituteSingleBaseData(num_base, num_array, num_interpolated,
				base_data, interpolated_data);
		return;
	}

	std::vector<StorageAndAlignedPointer<XDataType> > xdatatype_holder(2);
	StorageAndAlignedPointer<YDataType> ydatatype_holder;
	GetAscendingArray<Helper, XDataType, XDataType>(num_base, base_position, 1,
			base_position, &xdatatype_holder[0]);
	GetAscendingArray<Helper, XDataType, YDataType>(num_base, base_position,
			num_array, base_data, &ydatatype_holder);
	GetAscendingArray<Helper, XDataType, XDataType>(num_interpolated,
			interpolated_position, 1, interpolated_position,
			&xdatatype_holder[1]);
	XDataType const *base_position_work = xdatatype_holder[0].pointer;
	YDataType const *base_data_work = ydatatype_holder.pointer;
	XDataType const *interpolated_position_work = xdatatype_holder[1].pointer;

	// Perform 1-dimensional interpolation
	// Any preparation for interpolation should be done here
	std::unique_ptr<typename Interpolator::WorkingData> wdata_storage(
			Interpolator::PrepareForInterpolation(polynomial_order, num_base,
					num_array, base_position_work, base_data_work));
	typename Interpolator::WorkingData *work_data = wdata_storage.get();

	// Locate each element in x_base against x_interpolated
	StorageAndAlignedPointer<size_t> size_t_holder;
	AllocateAndAlign<size_t>(num_base, &size_t_holder);
	size_t *location_base = size_t_holder.pointer;
	size_t num_location_base = Locate<XDataType>(num_interpolated, num_base,
			interpolated_position_work, base_position_work, location_base);

	// Outside of x_base[0]
	Helper::SubstituteLeftMostData(location_base[0], num_base, num_array,
			num_interpolated, base_data_work, interpolated_data);

	// Between x_base[0] and x_base[num_x_base-1]
	size_t offset = 0;
	if (base_position_work[0] < interpolated_position_work[0]) {
		for (size_t i = 1; i < num_base; ++i) {
			if (base_position_work[offset + 1]
					< interpolated_position_work[0]) {
				offset++;
			} else {
				break;
			}
		}
	}
	Interpolator::Interpolate1D(num_base, base_position_work, num_array,
			base_data_work, num_interpolated, interpolated_position_work,
			interpolated_data, num_location_base, location_base, offset,
			work_data);

	// Outside of x_base[num_x_base-1]
	Helper::SubstituteRightMostData(location_base[num_location_base - 1],
			num_base, num_array, num_interpolated, base_data_work,
			interpolated_data);

	// swap output array
	if (interpolated_position[0]
			> interpolated_position[num_interpolated - 1]) {
		Helper::SwapResult(num_array, num_interpolated, interpolated_data);
	}
}

template<class InterpolatorHelper, class XDataType, class YDataType>
void ExecuteInterpolate1D(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_base,
		XDataType const base_position[], size_t num_array,
		YDataType const base_data[], size_t num_interpolated,
		XDataType const interpolated_position[],
		YDataType interpolated_data[]) {
	typedef void (*Interpolate1DFunc)(uint8_t, size_t, XDataType const *,
			size_t, YDataType const *, size_t, XDataType const *, YDataType *);
	Interpolate1DFunc func = [](uint8_t, size_t, XDataType const *,
			size_t, YDataType const *, size_t, XDataType const *, YDataType *) {
		throw LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	};
	switch (interpolation_method) {
	case LIBSAKURA_SYMBOL(InterpolationMethod_kNearest):
		func = Interpolate1D<typename InterpolatorHelper::NearestInterpolator,
				InterpolatorHelper, XDataType, YDataType>;
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kLinear):
		func = Interpolate1D<typename InterpolatorHelper::LinearInterpolator,
				InterpolatorHelper, XDataType, YDataType>;
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kPolynomial):
		if (polynomial_order == 0) {
			// This is special case: 0-th polynomial interpolation
			// acts like nearest interpolation
			func = Interpolate1D<
					typename InterpolatorHelper::NearestInterpolator,
					InterpolatorHelper, XDataType, YDataType>;
		} else {
			func = Interpolate1D<
					typename InterpolatorHelper::PolynomialInterpolator,
					InterpolatorHelper, XDataType, YDataType>;
		}
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kSpline):
		func = Interpolate1D<typename InterpolatorHelper::SplineInterpolator,
				InterpolatorHelper, XDataType, YDataType>;
		break;
	default:
		// invalid interpolation method type
		break;
	}
	(*func)(polynomial_order, num_base, base_position, num_array, base_data,
			num_interpolated, interpolated_position, interpolated_data);
}

// basic check of arguments
bool CheckArguments(LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_interpolation_axis,
		double const base[], size_t num_array, float const data_base[],
		size_t num_interpolated, double const interpolated[],
		float const data_interpolated[], LIBSAKURA_SYMBOL(
				Status)*status) {

	bool process_data = true;

	// check interpolation_method
	if (interpolation_method != LIBSAKURA_SYMBOL(InterpolationMethod_kNearest)
			&& interpolation_method
					!= LIBSAKURA_SYMBOL(InterpolationMethod_kLinear)
			&& interpolation_method
					!= LIBSAKURA_SYMBOL(InterpolationMethod_kPolynomial)
			&& interpolation_method
					!= LIBSAKURA_SYMBOL(InterpolationMethod_kSpline)) {
		LOG4CXX_ERROR(logger, "Invalid interpolation method");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		process_data = false;
	}

	// num_base must be non-zero
	if (num_interpolation_axis == 0) {
		LOG4CXX_ERROR(logger, "num_base must be >0");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		process_data = false;
	}

	// no interpolation will be done
	if (num_interpolated == 0 || num_array == 0) {
		// Nothing to do
		LOG4CXX_INFO(logger,
				"Nothing has been done since num_interpolated is 0");
		*status = LIBSAKURA_SYMBOL(Status_kOK);
		process_data = false;
	}

	// input arrays are not aligned
	if (!LIBSAKURA_SYMBOL(IsAligned)(base)
			|| !LIBSAKURA_SYMBOL(IsAligned)(data_base) || !LIBSAKURA_SYMBOL(
					IsAligned)(interpolated)
			|| !LIBSAKURA_SYMBOL(IsAligned)(data_interpolated)) {
		LOG4CXX_ERROR(logger, "input arrays are not aligned");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		process_data = false;
	}

	// input arrays are null
	if (base == nullptr || data_base == nullptr || interpolated == nullptr
			|| data_interpolated == nullptr) {
		LOG4CXX_ERROR(logger, "input arrays are null");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		process_data = false;
	}
	return process_data;
}

template<typename Func>
LIBSAKURA_SYMBOL(Status) DoInterpolate(Func func,
LIBSAKURA_SYMBOL(
		InterpolationMethod) interpolation_method, uint8_t polynomial_order,
		size_t num_base, double const base[/*num_x_base*/], size_t num_array,
		float const data_base[/*num_x_base*num_y*/], size_t num_interpolated,
		double const interpolated[/*num_x_interpolated*/],
		float data_interpolated[/*num_x_interpolated*num_y*/]) {
	// check arguments
	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Status_kOK);
	if (!CheckArguments(interpolation_method, polynomial_order, num_base, base,
			num_array, data_base, num_interpolated, interpolated,
			data_interpolated, &status)) {
		return status;
	}

	try {
		func(interpolation_method, polynomial_order, num_base, base, num_array,
				data_base, num_interpolated, interpolated, data_interpolated);
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
		float const data_base[/*num_x_base*num_y*/], size_t num_x_interpolated,
		double const x_interpolated[/*num_x_interpolated*/],
		float data_interpolated[/*num_x_interpolated*num_y*/]) {

	return DoInterpolate(
			ExecuteInterpolate1D<XInterpolatorHelper<double, float>, double,
					float>, interpolation_method, polynomial_order, num_x_base,
			x_base, num_y, data_base, num_x_interpolated, x_interpolated,
			data_interpolated);

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
		float const data_base[/*num_y_base*num_x*/], size_t num_y_interpolated,
		double const y_interpolated[/*num_y_interpolated*/],
		float data_interpolated[/*num_y_interpolated*num_x*/]) {

	return DoInterpolate(
			ExecuteInterpolate1D<YInterpolatorHelper<double, float>, double,
					float>, interpolation_method, polynomial_order, num_y_base,
			y_base, num_x, data_base, num_y_interpolated, y_interpolated,
			data_interpolated);

}
