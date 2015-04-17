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

#include <algorithm>
#include <cassert>
#include <climits>
#include <memory>
#include <utility>

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
				assert(num_location_list < num_data);
				location_list[num_location_list] = location;
				num_location_list += 1;
				start_position = location;
				previous_location = location;
			}
		}
	}
	return num_location_list;
}

// Working data
// NullWorkingData is a null implementation of working data,
// It has nothing and it does nothing
template<class XDataType, class YDataType>
struct NullWorkingData {
	static NullWorkingData<XDataType, YDataType> * Allocate(
			uint8_t polynomial_order, size_t num_base) {
		return nullptr;
	}
	static void Initialize(size_t num_base, XDataType const base_position[],
			YDataType const base_data[],
			NullWorkingData<XDataType, YDataType> *work_data) {
	}
};

struct CustomizedMemoryManagement {
	static void *operator new(size_t size) {
		//std::cout << "This is overloaded new" << std::endl;
		void *p = LIBSAKURA_PREFIX::Memory::Allocate(size);
		if (p == nullptr) {
			throw std::bad_alloc();
		}
		return p;
	}
	static void operator delete(void *p) {
		//std::cout << "This is overloaded delete" << std::endl;
		LIBSAKURA_PREFIX::Memory::Free(p);
	}
};

// Working data for polynomial interpolation
// It stores working array for Neville's algorithm
template<class XDataType, class YDataType>
struct PolynomialXWorkingData: CustomizedMemoryManagement {
	static PolynomialXWorkingData<XDataType, YDataType> * Allocate(
			uint8_t polynomial_order, size_t num_base) {
		return new PolynomialXWorkingData<XDataType, YDataType>(
				polynomial_order, num_base);
	}
	static void Initialize(size_t num_base, XDataType const base_position[],
			YDataType const base_data[],
			PolynomialXWorkingData<XDataType, YDataType> *work_datas) {
	}
	PolynomialXWorkingData(uint8_t order, size_t num_base) {
		assert(num_base > 0);
		polynomial_order = static_cast<uint8_t>(std::min(num_base - 1,
				static_cast<size_t>(order)));
		num_elements = polynomial_order + 1;
		AllocateAndAlign<XDataType>(num_elements, &(xholder[0]));
		AllocateAndAlign<XDataType>(num_elements, &(xholder[1]));
	}
	StorageAndAlignedPointer<XDataType> xholder[2];
	size_t num_elements;
	uint8_t polynomial_order;
};

// Working data for spline interpolation
// It stores working array for spline correction term
template<class XDataType, class YDataType>
struct SplineXWorkingData: public StorageAndAlignedPointer<YDataType>,
		public CustomizedMemoryManagement {
	static SplineXWorkingData<XDataType, YDataType> * Allocate(
			uint8_t polynomial_order, size_t num_base) {
		SplineXWorkingData<XDataType, YDataType> *work_data =
				new SplineXWorkingData<XDataType, YDataType>();
		try {
			AllocateAndAlign<YDataType>(num_base, work_data);
		} catch (...) {
			delete work_data;
			throw;
		}
		return work_data;
	}
	static void Initialize(size_t num_base, XDataType const base_position[],
			YDataType const base_data[],
			SplineXWorkingData<XDataType, YDataType> *work_data) {
		StorageAndAlignedPointer<YDataType> holder_for_u;
		AllocateAndAlign<YDataType>(num_base, &holder_for_u);
		YDataType *upper_triangular = holder_for_u.pointer;
		YDataType *d2ydx2 = work_data->pointer;

		// This is a condition of natural cubic spline
		d2ydx2[0] = YDataType(0.0);
		d2ydx2[num_base - 1] = YDataType(0.0);
		upper_triangular[0] = YDataType(0.0);

		// Solve tridiagonal system.
		// Here tridiagonal matrix is decomposed to upper triangular matrix.
		// upper_tridiangular stores upper triangular elements, while
		// d2ydx2 stores right-hand-side vector. The diagonal
		// elements are normalized to 1.

		// x_base is ascending order
		auto a1 = base_position[1] - base_position[0];
		assert(a1 != decltype(a1)(0));
		for (size_t i = 2; i < num_base; ++i) {
			auto const a2 = base_position[i] - base_position[i - 1];
			assert(a2 != decltype(a2)(0));
			XDataType const b1 = XDataType(1.0)
					/ (base_position[i] - base_position[i - 2]);
			size_t i0 = i;
			size_t i1 = i0 - 1;
			size_t i2 = i1 - 1;
			d2ydx2[i1] = YDataType(3.0) * b1
					* ((base_data[i0] - base_data[i1]) / a2
							- (base_data[i1] - base_data[i2]) / a1
							- d2ydx2[i2] * a1 / YDataType(2.0));
			XDataType a3 = XDataType(1.0)
					/ (XDataType(1.0)
							- upper_triangular[i2] * a1 * b1 / XDataType(2.0));
			d2ydx2[i1] *= a3;
			upper_triangular[i1] = a2 * b1 * a3 / YDataType(2.0);
			a1 = a2;
		}

		// Solve the system by backsubstitution and store solution to d2ydx2
		for (size_t k = num_base; k >= 3; --k) {
			d2ydx2[k - 2] -= upper_triangular[k - 2] * d2ydx2[k - 1];
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
struct NearestInterpolatorImpl: public InterpolatorInterface<
		NearestInterpolatorImpl<XDataType, YDataType>,
		NullWorkingData<XDataType, YDataType>, XDataType, YDataType> {
	typedef typename InterpolatorInterface<
			NearestInterpolatorImpl<XDataType, YDataType>,
			NullWorkingData<XDataType, YDataType>, XDataType, YDataType>::WorkingData WorkingData;
	static void DoInterpolate(size_t num_base, XDataType const base_position[],
			YDataType const base_data[], size_t num_interpolated,
			XDataType const interpolated_position[], size_t const location[],
			size_t k, size_t left_index, WorkingData const * const work_data,
			YDataType interpolated_data[]) {
		auto const midpoint = (base_position[left_index + 1]
				+ base_position[left_index])
				/ decltype(base_position[left_index])(2.0);
		auto const left_value = base_data[left_index];
		auto const right_value = base_data[left_index + 1];
		for (size_t i = location[k - 1]; i < location[k]; ++i) {
			interpolated_data[i] =
					(interpolated_position[i] > midpoint) ?
							right_value : left_value;
		}
	}
};

// TODO: documentation must be added
template<class XDataType, class YDataType>
struct LinearInterpolatorImpl: public InterpolatorInterface<
		LinearInterpolatorImpl<XDataType, YDataType>,
		NullWorkingData<XDataType, YDataType>, XDataType, YDataType> {
	typedef typename InterpolatorInterface<
			LinearInterpolatorImpl<XDataType, YDataType>,
			NullWorkingData<XDataType, YDataType>, XDataType, YDataType>::WorkingData WorkingData;
	static void DoInterpolate(size_t num_base, XDataType const base_position[],
			YDataType const base_data[], size_t num_interpolated,
			XDataType const interpolated_position[], size_t const location[],
			size_t k, size_t left_index, WorkingData const * const work_data,
			YDataType interpolated_data[]) {
		XDataType const dydx = static_cast<XDataType>(base_data[left_index + 1]
				- base_data[left_index])
				/ (base_position[left_index + 1] - base_position[left_index]);
		for (size_t i = location[k - 1]; i < location[k]; ++i) {
			interpolated_data[i] = base_data[left_index]
					+ static_cast<YDataType>(dydx
							* (interpolated_position[i]
									- base_position[left_index]));
		}
	}
};

// TODO: documentation must be added
template<class XDataType, class YDataType>
struct PolynomialInterpolatorImpl: public InterpolatorInterface<
		PolynomialInterpolatorImpl<XDataType, YDataType>,
		PolynomialXWorkingData<XDataType, YDataType>, XDataType, YDataType> {
	typedef typename InterpolatorInterface<
			PolynomialInterpolatorImpl<XDataType, YDataType>,
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
			interpolated_data[i] = PerformNevilleAlgorithm(num_base,
					base_position, base_data, left_edge,
					interpolated_position[i], work_data);
		}
	}
private:
	static YDataType PerformNevilleAlgorithm(size_t num_base,
			XDataType const base_position[], YDataType const base_data[],
			size_t left_index, XDataType interpolated_position,
			WorkingData const * const work_data) {

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
			for (size_t i = 0; i < work_data->num_elements - m; ++i) {
				XDataType cd = c[i + 1] - d[i];
				XDataType const dx = x_ptr[i] - x_ptr[i + m];
				assert(dx != 0);
				cd /= dx;
				c[i] = (x_ptr[i] - interpolated_position) * cd;
				d[i] = (x_ptr[i + m] - interpolated_position) * cd;
			}

			// In each step, c[0] holds Cm1 which is a correction between
			// P12...m and P12...[m+1]. Thus, the following repeated update
			// corresponds to the route P1 -> P12 -> P123 ->...-> P123...n.
			work += c[0];
		}

		return static_cast<YDataType>(work);
	}
};

// TODO: documentation must be added
template<class XDataType, class YDataType>
struct SplineInterpolatorImpl: public InterpolatorInterface<
		SplineInterpolatorImpl<XDataType, YDataType>,
		SplineXWorkingData<XDataType, YDataType>, XDataType, YDataType> {
	typedef typename InterpolatorInterface<
			SplineInterpolatorImpl<XDataType, YDataType>,
			SplineXWorkingData<XDataType, YDataType>, XDataType, YDataType>::WorkingData WorkingData;
	static void DoInterpolate(size_t num_base, XDataType const base_position[],
			YDataType const base_data[], size_t num_interpolated,
			XDataType const interpolated_position[], size_t const location[],
			size_t k, size_t left_index, WorkingData const * const work_data,
			YDataType interpolated_data[]) {
		YDataType const * const d2ydx2 = work_data->pointer;
		auto const dx = base_position[left_index + 1]
				- base_position[left_index];
		auto const dx_factor = dx * dx / decltype(dx)(6.0);
		for (size_t i = location[k - 1]; i < location[k]; ++i) {
			auto const a = (base_position[left_index + 1]
					- interpolated_position[i]) / dx;
			auto const b = decltype(a)(1.0) - a;
			interpolated_data[i] = static_cast<YDataType>(a
					* base_data[left_index] + b * base_data[left_index + 1]
					+ ((a * a * a - a) * d2ydx2[left_index]
							+ (b * b * b - b) * d2ydx2[left_index + 1])
							* dx_factor);
		}
	}
};

template<class Indexer, class XDataType, class YDataType>
inline size_t FillDataAsAscending(size_t num_column, XDataType const position[],
		size_t num_row, size_t the_row, YDataType const data[],
		bool const mask[], XDataType out_position[], YDataType out_data[]) {
	size_t n = 0;
	for (size_t i = 0; i < num_column; ++i) {
		size_t data_index = Indexer::GetIndex(num_column, num_row, the_row, i);
		size_t position_index = Indexer::GetIndex(num_column, 1, 0, i);
		if (mask[data_index]) {
			out_position[n] = position[position_index];
			out_data[n] = data[data_index];
			++n;
		}
	}
	return n;
}

template<class Indexer, class DataType>
inline void CopyResult(size_t num_column, size_t num_row, size_t the_row,
		DataType in_data[], DataType out_data[]) {
	for (size_t i = 0; i < num_column; ++i) {
		size_t index = Indexer::GetIndex(num_column, num_row, the_row, i);
		assert(index < num_row * num_column);
		out_data[index] = in_data[i];
	}
}

template<class DataType>
inline void FillOutOfRangeAreaWithEdgeValue(size_t edge_index,
		DataType const in[], size_t start, size_t end, DataType out[]) {
	for (size_t i = start; i < end; ++i) {
		out[i] = in[edge_index];
	}
}

template<class Indexer, class DataType>
static void FillOneRowWithValue(size_t num_column, size_t num_row,
		DataType value, size_t the_row,
		DataType data[/*num_row][num_column*/]) {
	assert(the_row < num_row);
	for (size_t i = 0; i < num_column; ++i) {
		size_t index = Indexer::GetIndex(num_column, num_row, the_row, i);
		assert(index < num_row * num_column);
		data[index] = value;
	}
}

struct XAscendingIndexer {
	static size_t GetIndex(size_t num_column, size_t num_row, size_t the_row,
			size_t index) {
		return num_column * the_row + index;
	}
};

struct XDescendingIndexer {
	static size_t GetIndex(size_t num_column, size_t num_row, size_t the_row,
			size_t index) {
		return num_column * the_row + num_column - 1 - index;
	}
};

struct YAscendingIndexer {
	static size_t GetIndex(size_t num_column, size_t num_row, size_t the_row,
			size_t index) {
		return the_row + num_row * index;
	}
};

struct YDescendingIndexer {
	static size_t GetIndex(size_t num_column, size_t num_row, size_t the_row,
			size_t index) {
		return the_row + num_row * (num_column - 1 - index);
	}
};

struct XInterpolatorHelper {
	typedef XAscendingIndexer AscendingIndexer;
	typedef XDescendingIndexer DescendingIndexer;
};

struct YInterpolatorHelper {
	typedef YAscendingIndexer AscendingIndexer;
	typedef YDescendingIndexer DescendingIndexer;
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
	typedef typename Helper::AscendingIndexer AscendingIndexer;
	typedef typename Helper::DescendingIndexer DescendingIndexer;
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

	// Define helper function here
	auto data_filler =
			(is_base_ascending) ?
					FillDataAsAscending<AscendingIndexer, XDataType, YDataType> :
					FillDataAsAscending<DescendingIndexer, XDataType, YDataType>;
	auto data_copier =
			(is_interp_ascending) ?
					CopyResult<AscendingIndexer, YDataType> :
					CopyResult<DescendingIndexer, YDataType>;

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
			WorkingData::Allocate(polynomial_order, num_base));
	WorkingData *work_data = wdata_storage.get();

	StorageAndAlignedPointer<size_t> size_t_holder;

	for (size_t iarray = 0; iarray < num_array; ++iarray) {
		// Pick up valid data and associating position from base_position and base_data
		// by referring base_mask. The data and position are reversed if sort order is
		// descending.
		size_t const n = data_filler(num_base, base_position, num_array, iarray,
				base_data, base_mask, x0, y0);
		if (n == 0) {
			// cannot perform interpolation, mask all data
			FillOneRowWithValue<AscendingIndexer, bool>(num_interpolated,
					num_array, false, iarray, interpolated_mask);
		} else if (n == 1) {
			// no need to interpolate, just substitute single base data
			// to all elements in interpolated data, keep input mask
			FillOneRowWithValue<AscendingIndexer, YDataType>(num_interpolated,
					num_array, y0[0], iarray, interpolated_data);
		} else {
			// perform interpolation, keep input mask

			// Locate each element in base_position against interpolated_position
			StorageAndAlignedPointer<size_t> size_t_holder;
			AllocateAndAlign<size_t>(n, &size_t_holder);
			size_t *location_base = size_t_holder.pointer;
			size_t const num_location_base = Locate<XDataType>(num_interpolated,
					n, x1, x0, location_base);

			// Outside of base_position[0]
			FillOutOfRangeAreaWithEdgeValue(0, y0, 0, location_base[0], y1);

			// Perform 1-dimensional interpolation
			// between base_position[0] and base_position[num_base-1]
			size_t offset = 0;
			while (offset + 1 < n && x0[offset + 1] < x1[0]) {
				offset++;
			}
			WorkingData::Initialize(n, x0, y0, work_data);
			Interpolator::Interpolate1D(n, x0, y0, num_interpolated, x1,
					num_location_base, location_base, offset, work_data, y1);

			// Outside of base_position[num_base-1]
			FillOutOfRangeAreaWithEdgeValue(n - 1, y0,
					location_base[num_location_base - 1], num_interpolated, y1);

			// Copy results in working array, y1, into result array
			data_copier(num_interpolated, num_array, iarray, y1,
					interpolated_data);
		}
	}
}

#ifndef NDEBUG
template<class DataType, class Checker>
bool IsDuplicateOrNotSorted(Checker checker, size_t num_elements,
		DataType const data[]) {
	for (size_t i = 0; i < num_elements - 1; ++i) {
		if (checker(data[i], data[i + 1])) {
			return true;
		}
	}
	return false;
}

template<class DataType>
bool IsDuplicateOrNotSorted(size_t num_elements, DataType const data[]) {
	if (num_elements <= 1) {
		return false;
	} else if (data[0] < data[num_elements - 1]) {
		// ascending check
		return IsDuplicateOrNotSorted(
				[](DataType x, DataType y) {return x >= y;}, num_elements,
				data);
	} else if (data[0] > data[num_elements - 1]) {
		// descending check
		return IsDuplicateOrNotSorted(
				[](DataType x, DataType y) {return x <= y;}, num_elements,
				data);
	}
	return true;
}
#endif

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

#ifndef NDEBUG
	// base_position is not sorted or base_position has duplicated element
	if (IsDuplicateOrNotSorted(num_base, base_position)) {
		LOG4CXX_ERROR(logger,
				"base_position is not sorted or has duplicated elements");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		return false;
	}

	// interpolated_position is not sorted or base_position has duplicated element
	if (IsDuplicateOrNotSorted(num_interpolated, interpolated_position)) {
		LOG4CXX_ERROR(logger,
				"interpolated_position is not sorted or has duplicated elements");
		*status = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		return false;
	}
#endif

	return true;
}

template<typename XDataType, typename YDataType, typename InterpolatorHelper>
LIBSAKURA_SYMBOL(Status) DoInterpolate(
LIBSAKURA_SYMBOL(
		InterpolationMethod) interpolation_method, uint8_t polynomial_order,
		size_t num_base, XDataType const base_position[/*num_base*/],
		size_t num_array, YDataType const base_data[/*num_base*num_array*/],
		bool const base_mask[/*num_base*num_array*/], size_t num_interpolated,
		XDataType const interpolated_position[/*num_interpolated*/],
		YDataType interpolated_data[/*num_interpolated*num_array*/],
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
		typedef NearestInterpolatorImpl<XDataType, YDataType> NearestInterpolator;
		typedef LinearInterpolatorImpl<XDataType, YDataType> LinearInterpolator;
		typedef PolynomialInterpolatorImpl<XDataType, YDataType> PolynomialInterpolator;
		typedef SplineInterpolatorImpl<XDataType, YDataType> SplineInterpolator;
		typedef void (*Interpolate1DFunc)(uint8_t, size_t, XDataType const *,
				size_t, YDataType const *, bool const *, size_t,
				XDataType const *, YDataType *, bool *);
		Interpolate1DFunc func = [](uint8_t, size_t, XDataType const *,
				size_t, YDataType const *, bool const *, size_t,
				XDataType const *, YDataType *, bool *) {
			throw LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		};
		switch (interpolation_method) {
		case LIBSAKURA_SYMBOL(InterpolationMethod_kNearest):
			func = Interpolate1D<NearestInterpolator, InterpolatorHelper,
					XDataType, YDataType>;
			break;
		case LIBSAKURA_SYMBOL(InterpolationMethod_kLinear):
			func = Interpolate1D<LinearInterpolator, InterpolatorHelper,
					XDataType, YDataType>;
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
			func = Interpolate1D<SplineInterpolator, InterpolatorHelper,
					XDataType, YDataType>;
			break;
		default:
			// invalid interpolation method type
			break;
		}
		(*func)(polynomial_order, num_base, base_position, num_array, base_data,
				base_mask, num_interpolated, interpolated_position,
				interpolated_data, interpolated_mask);
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

	return DoInterpolate<double, float, XInterpolatorHelper>(
			interpolation_method, polynomial_order, num_base, base_position,
			num_array, base_data, base_mask, num_interpolated,
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

	return DoInterpolate<double, float, YInterpolatorHelper>(
			interpolation_method, polynomial_order, num_base, base_position,
			num_array, base_data, base_mask, num_interpolated,
			interpolated_position, interpolated_data, interpolated_mask);

}
