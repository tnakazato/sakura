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
#include <libsakura/sakura.h>

#include <cassert>
#include <climits>
#include <memory>
#include <utility>
#include <vector>

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
template<class DataType>
size_t FindIndexOfClosestElementFromSortedArray(size_t start_index,
		size_t end_index, size_t num_reference, DataType const reference[],
		DataType data) {
	assert(num_reference > 0);
	assert(start_index <= end_index);
	assert(end_index < num_reference);
	assert(reference != nullptr);

	// If length of the array is just 1, return 0
	if (num_reference == 1)
		return 0;

	// base_position must be sorted in ascending order
	if (data <= reference[0]) {
		// out of range
		return 0;
	} else if (data > reference[num_reference - 1]) {
		// out of range
		return num_reference;
	} else if (data < reference[start_index]) {
		// x_located is not in the range (start_position, end_position)
		// call this function to search other location
		return FindIndexOfClosestElementFromSortedArray(0, start_index,
				num_reference, reference, data);
	} else if (data > reference[end_index]) {
		// located_position is not in the range (start_position, end_position)
		// call this function to search other location
		return FindIndexOfClosestElementFromSortedArray(end_index,
				num_reference - 1, num_reference, reference, data);
	} else {
		// do bisection
		size_t left_index = start_index;
		size_t right_index = end_index;
		while (left_index + 1 < right_index) {
			size_t const middle_index = (right_index + left_index) / 2;
			if (data > reference[middle_index]) {
				left_index = middle_index;
			} else {
				right_index = middle_index;
			}
		}
		return right_index;
	}
}

template<class DataType>
size_t Locate(size_t num_reference, size_t num_data, DataType const reference[],
		DataType const data[], size_t location_list[]) {
	// input arrays must be sorted in ascending order
	size_t num_location_list = 0;
	if (data[num_data - 1] <= reference[0]) {
		// all located_position's are on the left (lower) side of base_position
		num_location_list = 1;
		location_list[0] = 0;
	} else if (data[0] >= reference[num_reference - 1]) {
		// all located_position's are on the right (upper) side of base_position
		num_location_list = 1;
		location_list[0] = num_reference;
	} else if (num_reference == 1) {
		num_location_list = 2;
		location_list[0] = 0;
		location_list[1] = 1;
	} else {
		size_t start_position = 0;
		size_t end_position = num_reference - 1;
		size_t previous_location = num_reference + 1;
		for (size_t i = 0; i < num_data && start_position <= end_position;
				++i) {
			size_t const location = FindIndexOfClosestElementFromSortedArray(
					start_position, end_position, num_reference, reference,
					data[i]);
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
	static NullWorkingData<XDataType, YDataType> * Allocate(
			uint8_t polynomial_order, size_t num_base, size_t num_array) {
		return nullptr;
	}
	static void Initialize(size_t num_base, size_t num_array,
			XDataType const base_position[], YDataType const base_data[],
			NullWorkingData<XDataType, YDataType> *work_data) {
	}
};

template<class DataType>
struct PolynomialWorkingData {
	PolynomialWorkingData(uint8_t order, size_t num_base, size_t array_size,
			size_t num_array) :
			xholder(array_size) {
		assert(num_base > 0);
		polynomial_order = static_cast<uint8_t>(std::min(num_base - 1,
				static_cast<size_t>(order)));
		num_elements = polynomial_order + 1;
		for (size_t i = 0; i < array_size; ++i) {
			AllocateAndAlign<DataType>(num_elements * num_array, &(xholder[i]));
		}
	}
	std::vector<StorageAndAlignedPointer<DataType> > xholder;
	size_t num_elements;
	uint8_t polynomial_order;
};

template<class XDataType, class YDataType>
struct PolynomialXWorkingData: public PolynomialWorkingData<XDataType> {
	static PolynomialXWorkingData<XDataType, YDataType> * Allocate(
			uint8_t polynomial_order, size_t num_base, size_t num_array) {
		return new PolynomialXWorkingData<XDataType, YDataType>(
				polynomial_order, num_base, 2, 1);
	}
	static void Initialize(

	size_t num_base, size_t num_array, XDataType const base_position[],
			YDataType const base_data[],
			PolynomialXWorkingData<XDataType, YDataType> *work_datas) {
	}
	PolynomialXWorkingData(uint8_t order, size_t num_base, size_t array_size,
			size_t num_array) :
			PolynomialWorkingData<XDataType>(order, num_base, array_size,
					num_array) {
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
	static SplineXWorkingData<XDataType, YDataType> * Allocate(
			uint8_t polynomial_order, size_t num_base, size_t num_array) {
		SplineXWorkingData<XDataType, YDataType> *work_data =
				new SplineXWorkingData<XDataType, YDataType>();
		AllocateAndAlign<YDataType>(num_base * num_array, work_data);
		return work_data;
	}
	static void Initialize(size_t num_base, size_t num_array,
			XDataType const base_position[], YDataType const base_data[],
			SplineXWorkingData<XDataType, YDataType> *work_data) {
		work_data->DeriveSplineCorrectionTerm(num_base, base_position,
				num_array, base_data);
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

// Interface class for Interpolator
template<class InterpolatorImpl, class WorkData, class XDataType,
		class YDataType>
struct InterpolatorInterface {
	typedef WorkData WorkingData;
	static void Interpolate1D(size_t num_base, XDataType const base_position[],
			YDataType const base_data[], size_t num_interpolated,
			XDataType const interpolated_position[], size_t num_location,
			size_t const location[], size_t offset,
			WorkingData const * const work_data,
			YDataType interpolated_data[]) {
		for (size_t k = 1; k < num_location; ++k) {
			size_t const left_index = offset + k - 1;
			InterpolatorImpl::DoInterpolate(num_base, base_position, base_data,
					num_interpolated, interpolated_position, location, k,
					left_index, work_data, interpolated_data);
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
			YDataType const base_data[], size_t num_interpolated,
			XDataType const interpolated_position[], size_t const location[],
			size_t k, size_t left_index, WorkingData const * const work_data,
			YDataType interpolated_data[]) {
		XDataType const midpoint = 0.5
				* (base_position[left_index + 1] + base_position[left_index]);
		YDataType const left_value = base_data[left_index];
		YDataType const right_value = base_data[left_index + 1];
		//YDataType *work = &interpolated_data[0];
		for (size_t i = location[k - 1]; i < location[k]; ++i) {
			interpolated_data[i] =
					(interpolated_position[i] > midpoint) ?
							right_value : left_value;
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
			YDataType const base_data[], size_t num_interpolated,
			XDataType const interpolated_position[], size_t const location[],
			size_t k, size_t left_index, WorkingData const * const work_data,
			YDataType interpolated_data[]) {
		size_t const offset_index_left = left_index;
		XDataType const dydx =
				static_cast<XDataType>(base_data[offset_index_left + 1]
						- base_data[offset_index_left])
						/ (base_position[left_index + 1]
								- base_position[left_index]);
		YDataType *y_work = &interpolated_data[0];
		for (size_t i = location[k - 1]; i < location[k]; ++i) {
			y_work[i] = base_data[offset_index_left]
					+ static_cast<YDataType>(dydx
							* (interpolated_position[i]
									- base_position[left_index]));
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
			YDataType const base_data[], size_t num_interpolated,
			XDataType const interpolated_position[], size_t const location[],
			size_t k, size_t left_index, WorkingData const * const work_data,
			YDataType interpolated_data[]) {
		assert(num_base >= work_data->num_elements);
		size_t const left_edge2 = num_base - work_data->num_elements;
		size_t left_edge =
				(left_index >= work_data->polynomial_order / 2u) ?
						left_index - work_data->polynomial_order / 2u : 0;
		left_edge = (left_edge > left_edge2) ? left_edge2 : left_edge;
		for (size_t i = location[k - 1]; i < location[k]; ++i) {
			PerformNevilleAlgorithm(num_base, base_position, base_data,
					left_edge, interpolated_position[i], work_data,
					&interpolated_data[i]);
		}
	}
private:
	static void PerformNevilleAlgorithm(size_t num_base,
			XDataType const base_position[], YDataType const base_data[],
			size_t left_index, XDataType interpolated_position,
			WorkingData const * const work_data, YDataType *interpolated_data) {

		// working pointers
		XDataType const *x_ptr = &(base_position[left_index]);
		YDataType const *y_ptr = &(base_data[left_index]);

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

// TODO: documentation must be added
template<class XDataType, class YDataType>
struct SplineXInterpolatorImpl: public InterpolatorInterface<
		SplineXInterpolatorImpl<XDataType, YDataType>,
		SplineXWorkingData<XDataType, YDataType>, XDataType, YDataType> {
	typedef typename InterpolatorInterface<
			SplineXInterpolatorImpl<XDataType, YDataType>,
			SplineXWorkingData<XDataType, YDataType>, XDataType, YDataType>::WorkingData WorkingData;
	static void DoInterpolate(size_t num_base, XDataType const base_position[],
			YDataType const base_data[], size_t num_interpolated,
			XDataType const interpolated_position[], size_t const location[],
			size_t k, size_t left_index, WorkingData const * const work_data,
			YDataType interpolated_data[]) {
		YDataType const * const d2ydx2 = work_data->pointer;
		XDataType const dx = base_position[left_index + 1]
				- base_position[left_index];
		XDataType const dx_factor = dx * dx / 6.0;
		size_t const offset_index_left = left_index;
		for (size_t i = location[k - 1]; i < location[k]; ++i) {
			XDataType const a = (base_position[left_index + 1]
					- interpolated_position[i]) / dx;
			XDataType const b = 1.0 - a;
			interpolated_data[i] = static_cast<YDataType>(a
					* base_data[offset_index_left]
					+ b * base_data[offset_index_left + 1]
					+ ((a * a * a - a) * d2ydx2[offset_index_left]
							+ (b * b * b - b) * d2ydx2[offset_index_left + 1])
							* dx_factor);
		}
	}
};

template<class XDataType, class YDataType>
inline size_t FillDataAsAscendingImpl(size_t position_index_start,
		long long position_index_increment, size_t data_index_start,
		long long data_index_increment, size_t num_base,
		XDataType const position[], YDataType const data[],
		bool const mask[], XDataType x[], YDataType y[]) {
	size_t n = 0;
	for (size_t i = 0; i < num_base; ++i) {
		size_t data_index = data_index_start + data_index_increment * i;
		size_t position_index = position_index_start
				+ position_index_increment * i;
		if (mask[data_index]) {
			x[n] = position[position_index];
			y[n] = data[data_index];
			++n;
		}
	}
	return n;
}

template<class YDataType>
inline void FillResultImpl(size_t index_start, long long index_increment,
		size_t num_interpolated, YDataType y1[], YDataType data[]) {
	for (size_t i = 0; i < num_interpolated; ++i) {
		data[index_start + i * index_increment] = y1[i];
	}
}

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
		assert(num_base > 0);
		assert(iarray < num_array);
		size_t position_index_start = is_ascending ? 0 : num_base - 1;
		long long index_increment = is_ascending ? 1LL : -1LL;
		size_t data_index_start =
				is_ascending ? iarray * num_base : iarray * num_base + num_base - 1;
		assert(data_index_start < num_base * num_array);
		return FillDataAsAscendingImpl(position_index_start, index_increment,
				data_index_start, index_increment, num_base, position, data,
				mask, x, y);
	}

	template<class T>
	static void FillOneRowWithValue(size_t num_column, size_t num_row, T value,
			size_t the_row, T data[/*num_row][num_column*/]) {
		assert(the_row < num_row);
		for (size_t i = 0; i < num_column; ++i) {
			data[num_column * the_row + i] = value;
		}
	}

	static void FillResult(size_t num_interpolated, size_t num_array,
			size_t iarray, bool is_ascending, YDataType y1[],
			YDataType data[]) {
		assert(num_interpolated > 0);
		assert(iarray < num_array);
		size_t start =
				is_ascending ?
						iarray * num_interpolated :
						iarray * num_interpolated + num_interpolated - 1;
		assert(start < num_interpolated * num_array);
		long long increment = is_ascending ? 1LL : -1LL;
		FillResultImpl(start, increment, num_interpolated, y1, data);
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
		Substitute(num_array, 0, num_base, base_data, 0, num_interpolated,
				interpolated_data);
	}
	static inline void SubstituteLeftMostData(size_t location, size_t num_base,
			size_t num_array, size_t num_interpolated,
			YDataType const base_data[], YDataType interpolated_data[]) {
		Substitute(num_array, 0, num_base, base_data, 0, location,
				num_interpolated, interpolated_data);
	}
	static inline void SubstituteRightMostData(size_t location, size_t num_base,
			size_t num_array, size_t num_interpolated,
			YDataType const base_data[], YDataType interpolated_data[]) {
		Substitute(num_array, num_base - 1, num_base, base_data, location,
				num_interpolated, num_interpolated, interpolated_data);
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
	static inline void Substitute(size_t num_column, size_t source_row,
			size_t in_increment, YDataType const in[], size_t start_row,
			size_t end_row, size_t out_increment, YDataType out[]) {
		for (size_t j = 0; j < num_column; ++j) {
			for (size_t i = start_row; i < end_row; ++i) {
				out[j * out_increment + i] = in[j * in_increment + source_row];
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
		assert(num_base > 0);
		size_t position_index_start = is_ascending ? 0 : num_base - 1;
		long long position_index_increment = is_ascending ? 1LL : -1LL;
		size_t data_index_start =
				is_ascending ? iarray : num_array * (num_base - 1) + iarray;
		assert(num_array <= LLONG_MAX);
		long long data_index_increment = is_ascending ? num_array : -num_array;
		return FillDataAsAscendingImpl(position_index_start,
				position_index_increment, data_index_start,
				data_index_increment, num_base, position, data, mask, x, y);
	}

	template<class T>
	static void FillOneRowWithValue(size_t num_column, size_t num_row, T value,
			size_t the_row, T data[/*num_column][num_row*/]) {
		assert(the_row < num_row);
		for (size_t i = 0; i < num_column; ++i) {
			data[the_row + num_row * i] = value;
		}
	}

	static void FillResult(size_t num_interpolated, size_t num_array,
			size_t iarray, bool is_ascending, YDataType y1[],
			YDataType data[]) {
		assert(iarray < num_array);
		assert(num_interpolated > 0);
		size_t start =
				is_ascending ?
						iarray : num_array * (num_interpolated - 1) + iarray;
		assert(num_array < LLONG_MAX);
		long long increment = is_ascending ? num_array : -num_array;
		FillResultImpl(start, increment, num_interpolated, y1, data);
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
		Substitute(num_array, 0, base_data, 0, num_interpolated,
				interpolated_data);
	}
	static inline void SubstituteLeftMostData(size_t location, size_t num_base,
			size_t num_array, size_t num_interpolated,
			YDataType const base_data[], YDataType interpolated_data[]) {
		Substitute(num_array, 0, base_data, 0, location, interpolated_data);
	}
	static inline void SubstituteRightMostData(size_t location, size_t num_base,
			size_t num_array, size_t num_interpolated,
			YDataType const base_data[], YDataType interpolated_data[]) {
		Substitute(num_array, num_base - 1, base_data, location,
				num_interpolated, interpolated_data);
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
	static inline void Substitute(size_t num_column, size_t source_row,
			YDataType const in[], size_t start_row, size_t end_row,
			YDataType out[]) {
		for (size_t i = start_row; i < end_row; ++i) {
			for (size_t j = 0; j < num_column; ++j) {
				out[i * num_column + j] = in[source_row * num_column + j];
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
	assert(LIBSAKURA_SYMBOL(IsAligned)(base_mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(interpolated_mask));

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
		XDataType *x_writable = x1_storage.pointer;
		for (size_t i = 0; i < num_interpolated; ++i) {
			x_writable[i] = interpolated_position[num_interpolated - 1 - i];
		}
		x1 = x_writable;
	}
	AllocateAndAlign(num_interpolated, &y1_storage);
	YDataType *y1 = y1_storage.pointer;

	std::unique_ptr<WorkingData> wdata_storage(
			WorkingData::Allocate(polynomial_order, num_base, 1));
	WorkingData *work_data = wdata_storage.get();

	for (size_t iarray = 0; iarray < num_array; ++iarray) {
		size_t const n = Helper::FillDataAsAscending(num_base, num_array,
				iarray, base_position, base_data, base_mask, is_base_ascending,
				x0, y0);
		if (n == 0) {
			// cannot perform interpolation, mask all data
			Helper::FillOneRowWithValue(num_interpolated, num_array, false,
					iarray, interpolated_mask);
		} else if (n == 1) {
			// no need to interpolate, just substitute single base data
			// to all elements in interpolated data, keep input mask
			Helper::FillOneRowWithValue(num_interpolated, num_array, y0[0],
					iarray, interpolated_data);
		} else {
			// perform interpolation, keep input mask

			// Perform 1-dimensional interpolation
			// Any preparation for interpolation should be done here
			WorkingData::Initialize(n, 1, x0, y0, work_data);

			// Locate each element in base_position against interpolated_position
			StorageAndAlignedPointer<size_t> size_t_holder;
			AllocateAndAlign<size_t>(n, &size_t_holder);
			size_t *location_base = size_t_holder.pointer;
			size_t const num_location_base = Locate<XDataType>(num_interpolated,
					n, x1, x0, location_base);

			// Outside of base_position[0]
			Helper::SubstituteLeftMostData(location_base[0], n, 1,
					num_interpolated, y0, y1);

			// Between base_position[0] and base_position[num_base-1]
			size_t offset = 0;
			while (offset + 1 < n && x0[offset + 1] < x1[0]) {
				offset++;
			}
			Interpolator::Interpolate1D(n, x0, y0, num_interpolated, x1,
					num_location_base, location_base, offset, work_data, y1);

			// Outside of base_position[num_base-1]
			Helper::SubstituteRightMostData(
					location_base[num_location_base - 1], n, 1,
					num_interpolated, y0, y1);

			Helper::FillResult(num_interpolated, num_array, iarray,
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
		uint8_t polynomial_order, size_t num_base, double const base_position[],
		size_t num_array, float const base_data[],
		bool const base_mask[], size_t num_interpolated,
		double const interpolated_position[], float const interpolated_data[],
		bool const interpolated_mask[], LIBSAKURA_SYMBOL(
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
	if (num_base == 0) {
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
	if (!LIBSAKURA_SYMBOL(IsAligned)(base_position)
			|| !LIBSAKURA_SYMBOL(IsAligned)(base_data) || !LIBSAKURA_SYMBOL(
					IsAligned)(interpolated_position)
			|| !LIBSAKURA_SYMBOL(IsAligned)(interpolated_data)
			|| !LIBSAKURA_SYMBOL(IsAligned)(base_mask)
			|| !LIBSAKURA_SYMBOL(IsAligned)(interpolated_mask)) {
		LOG4CXX_ERROR(logger, "input arrays are not aligned");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		return false;
	}

	// input arrays are null
	if (base_position == nullptr || base_data == nullptr
			|| interpolated_position == nullptr || interpolated_data == nullptr
			|| base_mask == nullptr || interpolated_mask == nullptr) {
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
		size_t num_base, double const base_position[/*num_base*/],
		size_t num_array, float const base_data[/*num_base*num_array*/],
		bool const base_mask[/*num_base*num_array*/], size_t num_interpolated,
		double const interpolated_position[/*num_interpolated*/],
		float interpolated_data[/*num_interpolated*num_array*/],
		bool interpolated_mask[/*num_interpolated*num_array*/]) {
	// check arguments
	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Status_kOK);
	if (!IsValidArguments(interpolation_method, polynomial_order, num_base,
			base_position, num_array, base_data, base_mask, num_interpolated,
			interpolated_position, interpolated_data, interpolated_mask,
			&status)) {
		return status;
	}

	try {
		func(interpolation_method, polynomial_order, num_base, base_position,
				num_array, base_data, base_mask, num_interpolated,
				interpolated_position, interpolated_data, interpolated_mask);
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

}
/* anonymous namespace */

/**
 * InterpolateXAxisFloat performs 1D interpolation along column based on
 * base_position and base_data. base_data is a serial array of column-major
 * matrix data. Its memory layout is assumed to be:
 *     base_data[0]     = data[0][0]
 *     base_data[1]     = data[1][0]
 *     base_data[2]     = data[2][0]
 *     ...
 *     base_data[N]     = data[N][0]
 *     base_data[N+1]   = data[1][0]
 *     ...
 *     base_data[N*M-1] = data[N][M]
 * where N and M correspond to num_base and num_array respectively.
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InterpolateXAxisFloat)(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_base,
		double const base_position[/*num_base*/], size_t num_array,
		float const base_data[/*num_base*num_array*/],
		bool const base_mask[/*num_base*num_array*/], size_t num_interpolated,
		double const interpolated_position[/*num_interpolated*/],
		float interpolated_data[/*num_interpolated*num_array*/],
		bool interpolated_mask[/*num_interpolated*num_array*/]) noexcept {

	return DoInterpolate(
			ExecuteInterpolate1D<XInterpolatorHelper<double, float>, double,
					float>, interpolation_method, polynomial_order, num_base,
			base_position, num_array, base_data, base_mask, num_interpolated,
			interpolated_position, interpolated_data, interpolated_mask);

}

/**
 * Interpolate1DYAxisFloat performs 1D interpolation along column based on
 * base_position and base_data. base_data is a serial array of row-major
 * matrix data. Its memory layout is assumed to be:
 *     base_data[0]     = data[0][0]
 *     base_data[1]     = data[0][1]
 *     base_data[2]     = data[0][2]
 *     ...
 *     base_data[M]     = data[0][M]
 *     base_data[M+1]   = data[1][0]
 *     ...
 *     base_data[M*N-1] = data[N][M]
 * where N and M correspond to num_y_base and num_x respectively.
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InterpolateYAxisFloat)(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_base,
		double const base_position[/*num_base*/], size_t num_array,
		float const base_data[/*num_base*num_array*/],
		bool const base_mask[/*num_base*num_array*/], size_t num_interpolated,
		double const interpolated_position[/*num_interpolated*/],
		float interpolated_data[/*num_interpolated*num_array*/],
		bool interpolated_mask[/*num_interpolated*num_array*/]) noexcept {

	return DoInterpolate(
			ExecuteInterpolate1D<YInterpolatorHelper<double, float>, double,
					float>, interpolation_method, polynomial_order, num_base,
			base_position, num_array, base_data, base_mask, num_interpolated,
			interpolated_position, interpolated_data, interpolated_mask);

}
