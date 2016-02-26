/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2016
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
#include <algorithm>
#include <cassert>
#include <climits>
#include <memory>
#include <type_traits>
#include <utility>

#include <libsakura/localdef.h>
#include <libsakura/logger.h>
#include <libsakura/memory_manager.h>
#include <libsakura/sakura.h>

namespace {

// a logger for this module
auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("interpolation");

/**
 * @struct StorageAndAlignedPointer
 * @brief Storage for aligned array
 * @details
 * @a StorageAndAlignedPointer holds two pointers, @a storage and @a pointer.
 * @a storage points an initial address of allocated memory area while
 * @a pointer keeps a minimum aligned address of allocated memory area.
 * Its destructor will automatically release allocated area that @a storage points to.
 *
 * @tparam DataType data type
 */
template<class DataType>
struct StorageAndAlignedPointer {
	/**
	 * Constructor. It initializes @a storage and @a pointer by nullptr.
	 */
	StorageAndAlignedPointer() :
			storage(nullptr), pointer(nullptr) {
	}
	/**
	 * Destructur. It releases allocated area that @a storage points to.
	 */
	~StorageAndAlignedPointer() {
		if (storage != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(storage);
			storage = nullptr;
		}
	}
	/**
	 * An initial address of allocated memory area.
	 */
	void *storage;
	/**
	 * A minimum aligned address of allocated memory area.
	 */
	DataType *pointer;
};

/**
 * AllocateAndAlign allocates required memory area that is necessary to get aligned array
 * with number of elements specified by @a num_elements. Required size with mergin for
 * alignment will autocatically be calculated. Throws the std::bad_alloc if allocation
 * failed or size of allocated memory area is not enough.
 * Initial address and aligned address will be held by @a holder.
 *
 * @tparam DataType data type
 *
 * @param[in] num_elements number of elements of aligned array
 * @param[in] holder[out] struct holding both initial address and aligned address of
 * allocated memory area
 * @exception std::bad_alloc failed to allocate memory
 */
template<class DataType>
inline void AllocateAndAlign(size_t num_elements,
		StorageAndAlignedPointer<DataType> *holder) {
	holder->storage = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
			num_elements * sizeof(DataType), &(holder->pointer));
}

/**
 * For given @a data, it returns the array index, @a i,  for @a reference whose
 * element @a reference[i-1] and @a reference[i] brackets @a data. Reference
 * array @a reference must be sorted in an ascending order. If @a data is less
 * than @a reference[0], 0 will be returned. Also, if @a data is greater than
 * @a reference[num_reference-1], @a num_reference will be returned.
 *
 * @tparam DataType data type for @a reference and @a data
 *
 * @param[in] start_index start index of @a reference for the search
 * @param[in] end_index end index of @a reference for the search
 * @param[in] num_reference number of elements for @a reference
 * @param[in] reference reference array. Its length must be @a num_reference
 * @param[in] data data to be located.
 * @return array index for @a reference that brackets @a data
 */
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
	assert(reference[0] < reference[num_reference - 1]);

	if (data <= reference[0]) {
		// out of range
		return 0;
	} else if (data > reference[num_reference - 1]) {
		// out of range
		return num_reference;
	} else if (data < reference[start_index]) {
		// data is not in the range between start_index and end_index
		// call this function to search other location (< start_index)
		return FindIndexOfClosestElementFromSortedArray(0, start_index,
				num_reference, reference, data);
	} else if (data > reference[end_index]) {
		// data is not in the range between start_index and end_index
		// call this function to search other location (> end_index)
		return FindIndexOfClosestElementFromSortedArray(end_index,
				num_reference - 1, num_reference, reference, data);
	} else {
		// data will be located between start_index and end_index
		// perform bisection search
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

/**
 * Compare elements in @a data with the ones in @a reference and determine
 * where each element in @a data is located in between the elements in
 * @a reference.
 * The result is stored in @a location_list as a upper index of @a reference
 * that brackets the element in @a data. The index 0 or @a num_reference in
 * @a location_list indicate that the data is out of range.
 *
 * For example, supporse that @a data and @a reference are as follows:
 * @verbatim
 * float data[] = {-0.5, 1.8, 1.9, 2.2};
 * float reference[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
 * @endverbatim
 * The @a location_list will store the indices 0, 4, and 5. And the
 * @a lower_index_list will be 0, 2, and 3. These results are visually
 * shown below.
 *
 * @verbartim
 *
 * data:         0                       1    2         3
 * reference:        0    1    2    3              4         5    6
 *               |   x    x    x    x    |    |    x    |    x    x
 * location:         ^                             ^         ^
 * lower_index:  ^                            ^         ^
 *
 * @endverbatim
 *
 * @pre @a data and @a reference must be sorted in an ascending order
 *
 * @tparam DataType data type for @a reference and @a data
 *
 * @param[in] num_reference number of reference
 * @param[in] reference reference. Number of elements must be @a num_reference.
 * must-be-aligned
 * @param[in] num_data number of data
 * @param[in] data data. Number of elements must be @a num_data.
 * must-be-aligned
 * @param[out] location_list location list. Number of elements must be @a num_data.
 * must-be-aligned
 * @param[out] lower_index_list index list for lower side data with respect to
 * reference locations that are specified by @a location_list. must-be-aligned
 * @return number of locations
 */
template<class DataType>
size_t Locate(size_t num_reference, DataType const reference[], size_t num_data,
		DataType const data[], size_t location_list[],
		size_t lower_index_list[]) {
	// input arrays must be sorted in ascending order and no duplicates
	assert(num_reference == 1 || reference[0] < reference[num_reference - 1]);
	assert(num_data == 1 || data[0] < data[num_data - 1]);

	size_t num_location_list = 0;
	if (data[num_data - 1] <= reference[0]) {
		// all data are on the left (lower) side of reference
		num_location_list = 1;
		location_list[0] = 0;
		lower_index_list[0] = 0;
	} else if (data[0] >= reference[num_reference - 1]) {
		// all data are on the right (upper) side of reference
		num_location_list = 1;
		location_list[0] = num_reference;
		lower_index_list[0] = 0;
	} else if (num_reference == 1) {
		// data are divided at a certain point by reference
		num_location_list = 2;
		location_list[0] = 0;
		location_list[1] = 1;
		size_t const rloc = FindIndexOfClosestElementFromSortedArray(0,
				num_data - 1, num_data, data, reference[0]);
		lower_index_list[0] = rloc - ((data[rloc] == reference[0]) ? 0 : 1);
		lower_index_list[1] = num_data - 1;
	} else {
		// Evaluate location by FindIndexOfClosestElementFromSortedArray, which
		// basically performs bisection search to locate data
		// start_position and end_position are the cache for searched range
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
				if (location < num_reference) {
					size_t const rloc =
							FindIndexOfClosestElementFromSortedArray(0,
									num_data - 1, num_data, data,
									reference[location]);
					lower_index_list[num_location_list] = rloc
							- ((data[rloc] == reference[location]) ? 0 : 1);
				} else {
					lower_index_list[num_location_list] = num_data - 1;
				}
				num_location_list += 1;
				start_position = location;
				previous_location = location;
			}
		}
	}
	return num_location_list;
}

// Working data
// Working data should have static methods, Allocate and Initialize.
// It keeps necessary working data inside.

/**
 * @a NullWorkingData is a null implementation of working data.
 * It has nothing and it does nothing.
 *
 * @tparam XDataType data type for position
 * @tparam YDataType data type for data
 */
template<class XDataType, class YDataType>
struct NullWorkingData {
	/**
	 * The @a Allocate method does nothing.
	 * @param[in] polynomial_order no effect
	 * @param[in] num_base no effect
	 * @return nullptr
	 */
	static NullWorkingData<XDataType, YDataType> * Allocate(
			uint8_t polynomial_order, size_t num_base) {
		return nullptr;
	}
	/**
	 * The @a Initialize method does nothing.
	 * @param[in] num_base no effect
	 * @param[in] base_position no effect
	 * @param[in] base_data no effect
	 * @param[out] work_data no effect
	 */
	static void Initialize(size_t num_base, XDataType const base_position[],
			YDataType const base_data[],
			NullWorkingData<XDataType, YDataType> *work_data) {
	}
};

/**
 * @a CustomizedMemoryManagement defines overloaded operator new and delete
 * that depends on user-specified allocation and deallocation functions.
 * Working data can use overloaded new/delete by inheriting this class.
 */
struct CustomizedMemoryManagement {
	/**
	 * Overloaded operator new. It uses
	 * @link sakura::Memory::Allocate Memory::Allocate @endlink
	 * to allocate memory.
	 * @param[in] size size of allocated memory
	 */
	static void *operator new(size_t size) {
		//std::cout << "This is overloaded new" << std::endl;
		void *p = LIBSAKURA_PREFIX::Memory::Allocate(size);
		if (p == nullptr) {
			throw std::bad_alloc();
		}
		return p;
	}
	/**
	 * Overloaded operator delete. It uses
	 * @link sakura::Memory::Free Memory::Free @endlink
	 * to release memory.
	 * @param[in] p initial address of the allocated memory
	 */
	static void operator delete(void *p) {
		//std::cout << "This is overloaded delete" << std::endl;
		LIBSAKURA_PREFIX::Memory::Free(p);
	}
};

/**
 * @a PolynomialWorkingData is a working data for polynomial interpolation.
 * It stores working array for Neville's algorithm.
 *
 * @tparam XDataType data type for position
 * @tparam YDataType data type for data
 */
template<class XDataType, class YDataType>
struct PolynomialWorkingData: CustomizedMemoryManagement {
	typedef typename std::common_type<XDataType, YDataType>::type WDataType;

	/**
	 * The @a Allocate method creates @a PolynomialWorkingData instance and return it.
	 * @param[in] polynomial_order polynomial order for polynomial interpolation
	 * @param[in] num_base number of base data
	 * @return working data
	 */
	static PolynomialWorkingData<XDataType, YDataType> * Allocate(
			uint8_t polynomial_order, size_t num_base) {
		return new PolynomialWorkingData<XDataType, YDataType>(polynomial_order,
				num_base);
	}
	/**
	 * The @a Initialize method does nothing.
	 * @param[in] num_base no effect
	 * @param[in] base_position no effect
	 * @param[in] base_data no effect
	 * @param[out] work_datas no effect
	 */
	static void Initialize(size_t num_base, XDataType const base_position[],
			YDataType const base_data[],
			PolynomialWorkingData<XDataType, YDataType> *work_datas) {
	}
	/**
	 * Constructor allocates working array for polynomial interpolation.
	 * @param[in] order polynomial order for polynomial interpolation
	 * @param[in] num_base number of base data
	 */
	PolynomialWorkingData(uint8_t order, size_t num_base) {
		assert(num_base > 0);
		polynomial_order = static_cast<uint8_t>(std::min(num_base - 1,
				static_cast<size_t>(order)));
		num_elements = polynomial_order + 1;
		AllocateAndAlign<WDataType>(num_elements, &(xholder[0]));
		AllocateAndAlign<WDataType>(num_elements, &(xholder[1]));
	}
	/**
	 * Working arrays for polynomial interpolation. The arrays are kept as
	 * @a StorageAndAlignedPointer instance since working arrays should be
	 * aligned.
	 */
	StorageAndAlignedPointer<WDataType> xholder[2];
	/**
	 * Number of elements of each working array.
	 */
	size_t num_elements;
	/**
	 * Specified polynomial order for the polynomial interpolation.
	 */
	uint8_t polynomial_order;
};

/**
 * @a SplineWorkingData is a working data for spline interpolation.
 * It stores working array for spline correction term.
 *
 * @tparam XDataType data type for position
 * @tparam YDataType data type for data
 */
template<class XDataType, class YDataType>
struct SplineWorkingData: public CustomizedMemoryManagement {
	typedef typename std::common_type<XDataType, YDataType>::type WDataType;

	/**
	 * The @a Allocate method creates @a SplineWorkingData instance and return it.
	 * @param[in] polynomial_order no effect
	 * @param[in] num_base number of base data
	 * @return working data
	 */
	static SplineWorkingData<XDataType, YDataType> * Allocate(
			uint8_t polynomial_order, size_t num_base) {
		SplineWorkingData<XDataType, YDataType> *work_data =
				new SplineWorkingData<XDataType, YDataType>(num_base);
		return work_data;
	}
	/**
	 * The @a Initialize method computes spline correction term.
	 *
	 * Condition of the spline correction term is that the first derivative is smooth
	 * and the second derivative is contiguous within the interpolation range.
	 * Spline curve of each segment can be constructed from linearly interpolated
	 * curve, which is a straight line that connects two points, by adding spline
	 * correction term.
	 *
	 * Spline correction term is a cubic polynomial that satisfies the following
	 * conditions:
	 *
	 *    -# its value must be zero at each data point
	 *    -# its second derivative linearly changes between two base points
	 *
	 * These conditions can be formalized as a tri-diagonal system with a certain
	 * boundary condition, which defines the second derivatives of edge points.
	 * The solution of the system is the second derivatives of each base points.
	 * The method solves tri-diagonal system with the natural cubic spline condition,
	 * which defines the second derivatives of edge points as zero.
	 *
	 * @param[in] num_base number of base data
	 * @param[in] base_position positions of base data. Number of elements must be
	 * @a num_base. must-be-aligned
	 * @param[in] base_data base data. Number of elements must be @a num_base.
	 * must-be-aligned
	 * @param[out] work_datas working data
	 */
	static void Initialize(size_t num_base, XDataType const base_position[],
			YDataType const base_data[],
			SplineWorkingData<XDataType, YDataType> *work_data) {
		auto upper_triangular = work_data->upper_triangular.pointer;
		auto d2ydx2 = work_data->second_derivative.pointer;

		// This is a condition of natural cubic spline: second derivative of
		// interpolated spline curve should be zero.
		d2ydx2[0] = WDataType(0);
		d2ydx2[num_base - 1] = WDataType(0);
		upper_triangular[0] = WDataType(0);

		// Solve tridiagonal system.
		// Here tridiagonal matrix is decomposed to upper triangular matrix.
		// upper_tridiangular stores upper triangular elements, while
		// d2ydx2 stores right-hand-side vector. The diagonal
		// elements are normalized to 1.

		// base_position is ascending order
		constexpr auto kOne = WDataType(1);
		constexpr auto kTwo = WDataType(2);
		constexpr auto kThree = WDataType(3);
		auto a1 = WDataType(base_position[1] - base_position[0]);
		assert(a1 != decltype(a1)(0));
		auto d1 = WDataType(base_data[1] - base_data[0]);
		for (size_t i = 2; i < num_base; ++i) {
			auto const a2 = WDataType(base_position[i] - base_position[i - 1]);
			assert(a2 != decltype(a2)(0));
			auto const b1 = WDataType(base_position[i] - base_position[i - 2]);
			assert(b1 != decltype(b1)(0));
			auto const d2 = WDataType(base_data[i] - base_data[i - 1]);
			auto const u = a2 / (kTwo * b1);
			auto const l = a1 / (kTwo * b1);
			auto const r = (d2 / a2 - d1 / a1) * kThree / b1;
			auto const denominator = kOne - upper_triangular[i - 2] * l;
			upper_triangular[i - 1] = u / denominator;
			d2ydx2[i - 1] = (r - l * d2ydx2[i - 2]) / denominator;
			a1 = a2;
			d1 = d2;
		}

		// Solve the system by backsubstitution and store solution to d2ydx2
		for (size_t k = num_base; k >= 3; --k) {
			d2ydx2[k - 2] -= upper_triangular[k - 2] * d2ydx2[k - 1];
		}
	}
	/**
	 * Constructor allocates memory for spline correction term.
	 * @param[in] num_base number of elements of spline correction term
	 */
	SplineWorkingData(size_t num_base) {
		AllocateAndAlign<WDataType>(num_base, &second_derivative);
		AllocateAndAlign<WDataType>(num_base, &upper_triangular);
	}
	/**
	 * Working array for spline interpolation. The array is kept as
	 * @a StorageAndAlignedPointer instance since working arrays should be
	 * aligned. It will store a second derivative of spline curve for each
	 * base position.
	 */
	StorageAndAlignedPointer<WDataType> second_derivative;
	/**
	 * Working array for spline interpolation. The array is kept as
	 * @a StorageAndAlignedPointer instance since working arrays should be
	 * aligned. It will store elements of upper triangular matrix that
	 * is necessary to calculate spline correction term.
	 */
	StorageAndAlignedPointer<WDataType> upper_triangular;
};

/**
 * Interface class for Interpolator. It defines an template method, @a Interpolate1D,
 * that will be called from @link ::Interpolate1D Interpolate1D @endlink.
 * Implementation classes should inherit it according to CRTP and implement
 * @a DoInterpolate method.
 *
 * @tparam InterolatorImpl implementation class of interpolator
 * @tparam WorkData working data type that corresponds to @a InterpolatorImpl
 * @tparam XDataType data type for position
 * @tparam YDataType data type for data
 */
template<class InterpolatorImpl, class WorkData, class XDataType,
		class YDataType>
struct InterpolatorInterface {
	/**
	 * Type of working data
	 */
	typedef WorkData WorkingData;
	/**
	 * Template method for interpolation. It calls @a DoInterpolate method of its
	 * implementation class for each interpolation position.
	 *
	 * @param[in] num_base number of base data
	 * @param[in] base_position positions of base data. Number of elemments must be
	 * @a num_data. must-be-aligned
	 * @param[in] base_data base data. Number of elements must be @a num_data.
	 * must-be-aligned
	 * @param[in] num_interpolated number of interpolation data
	 * @param[in] interpolated_position positions of interpolation data. Number of elements
	 * must be @a num_interpolated. must-be-aligned
	 * @param[in] num_location Length of @a location.
	 * @param[in] location list of indices that indicates the location of each position
	 * in @a base_position with respect to @a interpolated_position. Its length must be
	 * @a num_location. must-be-aligned
	 * @param[in] lower_index index array of lower side of base data for interpolated data
	 * segment. lower_index[k] is an index for segment between location[k] and location[k+1].
	 * must-be-aligned
	 * @param[in] work_data working data for interpolation
	 * @param[out] interpolated_data interpolated value for each element whose location
	 * is stored in @a interpolated_position
	 */
	template<class Indexer>
	static void Interpolate1D(size_t num_base, XDataType const base_position[],
			YDataType const base_data[], size_t num_interpolated,
			XDataType const interpolated_position[], size_t num_location,
			size_t const location[], size_t const lower_index[],
			WorkingData const * const work_data, size_t num_array,
			size_t iarray, YDataType interpolated_data[]) {
		for (size_t k = 1; k < num_location; ++k) {
			InterpolatorImpl::template DoInterpolate<Indexer>(num_base,
					base_position, base_data, num_interpolated,
					interpolated_position, location[k - 1], location[k],
					lower_index[k - 1], work_data, num_array, iarray,
					interpolated_data);
		}
	}
};

/**
 * Nearest interpolation implementation.
 *
 * @tparam XDataType data type for position
 * @tparam YDataType data type for data
 */
template<class XDataType, class YDataType>
struct NearestInterpolatorImpl: public InterpolatorInterface<
		NearestInterpolatorImpl<XDataType, YDataType>,
		NullWorkingData<XDataType, YDataType>, XDataType, YDataType> {
	typedef typename InterpolatorInterface<
			NearestInterpolatorImpl<XDataType, YDataType>,
			NullWorkingData<XDataType, YDataType>, XDataType, YDataType>::WorkingData WorkingData;
	/**
	 * Perform nearest interpolation.
	 * For each @a interpolated_position that locates in between @a base_position[lower_index]
	 * and @a base_position[lower_index+1], perform nearest interpolation and set the result
	 * to @a interpolated_data. Interpolated result will be @a base_data[lower_index] if
	 * @a interpolated_position is closer to @a base_position[lower_index] or @a interpolated_position
	 * is located at the midpoint. Otherwise, the result will be @a base_data[lower_index+1].
	 *
	 * @param[in] num_base number of elements for data points.
	 * @param[in] base_position position of data points. Its length must be @a num_base.
	 * must-be-aligned
	 * @param[in] base_data value of data points. Its length must be @a num_base.
	 * must-be-aligned
	 * @param[in] num_interpolated number of elements for points that wants to get
	 * interpolated value.
	 * @param[in] interpolated_position location of points that wants to get interpolated
	 * value. Its length must be @a num_interpolated.
	 * @param[in] index_from start index for @a interpolated_position and @a interpolated_data
	 * that are located in between @a base_position[lower_index] and @a base_position[lower_index+1]
	 * @param[in] index_to end index for @a interpolated_position and @a interpolated_data
	 * that are located in between @a base_position[lower_index] and @a base_position[lower_index+1]
	 * @param[in] lower_index index for @a base_position. Interpolated area is defined as
	 * @a base_position[lower_index] and @a base_position[lower_index+1],
	 * @param[in] work_data working data for interpolation
	 * @param[out] interpolated_data storage for interpolation result. Its length must be
	 * @a num_interpolated. must-be-aligned
	 */
	template<class Indexer>
	static void DoInterpolate(size_t num_base, XDataType const base_position[],
			YDataType const base_data[], size_t num_interpolated,
			XDataType const interpolated_position[], size_t index_from,
			size_t index_to, size_t lower_index,
			WorkingData const * const work_data, size_t num_array,
			size_t iarray, YDataType interpolated_data[]) {
		auto const midpoint = (base_position[lower_index + 1]
				+ base_position[lower_index])
				/ decltype(base_position[lower_index])(2.0);
		auto const left_value = base_data[lower_index];
		auto const right_value = base_data[lower_index + 1];
		for (size_t i = index_from; i < index_to; ++i) {
			interpolated_data[Indexer::GetIndex(num_interpolated, num_array,
					iarray, i)] =
					(interpolated_position[i] > midpoint) ?
							right_value : left_value;
		}
	}
};

/**
 * Linear interpolation implementation.
 *
 * @tparam XDataType data type for position
 * @tparam YDataType data type for data
 */
template<class XDataType, class YDataType>
struct LinearInterpolatorImpl: public InterpolatorInterface<
		LinearInterpolatorImpl<XDataType, YDataType>,
		NullWorkingData<XDataType, YDataType>, XDataType, YDataType> {
	typedef typename InterpolatorInterface<
			LinearInterpolatorImpl<XDataType, YDataType>,
			NullWorkingData<XDataType, YDataType>, XDataType, YDataType>::WorkingData WorkingData;
	typedef typename std::common_type<XDataType, YDataType>::type WDataType;
	/**
	 * Perform linear interpolation.
	 * For each @a interpolated_position that locates in between @a base_position[lower_index]
	 * and @a base_position[lower_index+1], perform linear interpolation and set the result
	 * to @a interpolated_data. Interpolated curve is a straight line that goes through
	 * two data points (@a base_position[lower_index], @a base_data[lower_index]) and
	 * (@a base_position[lower_index+1], @a base_data[lower_index+1]).
	 *
	 * @param[in] num_base number of elements for data points.
	 * @param[in] base_position position of data points. Its length must be @a num_base.
	 * must-be-aligned
	 * @param[in] base_data value of data points. Its length must be @a num_base.
	 * must-be-aligned
	 * @param[in] num_interpolated number of elements for points that wants to get
	 * interpolated value.
	 * @param[in] interpolated_position location of points that wants to get interpolated
	 * value. Its length must be @a num_interpolated.
	 * @param[in] index_from start index for @a interpolated_position and @a interpolated_data
	 * that are located in between @a base_position[lower_index] and @a base_position[lower_index+1]
	 * @param[in] index_to end index for @a interpolated_position and @a interpolated_data
	 * that are located in between @a base_position[lower_index] and @a base_position[lower_index+1]
	 * @param[in] lower_index index for @a base_position. Interpolated area is defined as
	 * @a base_position[lower_index] and @a base_position[lower_index+1],
	 * @param[in] work_data working data for interpolation
	 * @param[out] interpolated_data storage for interpolation result. Its length must be
	 * @a num_interpolated. must-be-aligned
	 */
	template<class Indexer>
	static void DoInterpolate(size_t num_base, XDataType const base_position[],
			YDataType const base_data[], size_t num_interpolated,
			XDataType const interpolated_position[], size_t index_from,
			size_t index_to, size_t lower_index,
			WorkingData const * const work_data, size_t num_array,
			size_t iarray, YDataType interpolated_data[]) {
		WDataType const dydx = (base_data[lower_index + 1]
				- base_data[lower_index])
				/ (base_position[lower_index + 1] - base_position[lower_index]);
		auto const position_lower = WDataType(base_position[lower_index]);
		auto const data_lower = WDataType(base_data[lower_index]);
		for (size_t i = index_from; i < index_to; ++i) {
			interpolated_data[Indexer::GetIndex(num_interpolated, num_array,
					iarray, i)] = static_cast<YDataType>(data_lower
					+ dydx
							* (static_cast<WDataType>(interpolated_position[i])
									- position_lower));
		}
	}
};

/**
 * Polynomial interpolation implementation.
 *
 * @tparam XDataType data type for position
 * @tparam YDataType data type for data
 */
template<class XDataType, class YDataType>
struct PolynomialInterpolatorImpl: public InterpolatorInterface<
		PolynomialInterpolatorImpl<XDataType, YDataType>,
		PolynomialWorkingData<XDataType, YDataType>, XDataType, YDataType> {
	typedef typename InterpolatorInterface<
			PolynomialInterpolatorImpl<XDataType, YDataType>,
			PolynomialWorkingData<XDataType, YDataType>, XDataType, YDataType>::WorkingData WorkingData;
	typedef typename std::common_type<XDataType, YDataType>::type WDataType;
	/**
	 * Perform linear interpolation.
	 * For each @a interpolated_position that locates in between @a base_position[lower_index]
	 * and @a base_position[lower_index+1], perform linear interpolation and set the result
	 * to @a interpolated_data. Interpolated curve is a polynomial that goes through given
	 * data points. The polynomial is evaluated by the Nevilles algorithm.
	 *
	 * See http://en.wikipedia.org/wiki/Neville's_algorithm for details about the algorithm.
	 *
	 * @param[in] num_base number of elements for data points.
	 * @param[in] base_position position of data points. Its length must be @a num_base.
	 * must-be-aligned
	 * @param[in] base_data value of data points. Its length must be @a num_base.
	 * must-be-aligned
	 * @param[in] num_interpolated number of elements for points that wants to get
	 * interpolated value.
	 * @param[in] interpolated_position location of points that wants to get interpolated
	 * value. Its length must be @a num_interpolated.
	 * @param[in] index_from start index for @a interpolated_position and @a interpolated_data
	 * that are located in between @a base_position[lower_index] and @a base_position[lower_index+1]
	 * @param[in] index_to end index for @a interpolated_position and @a interpolated_data
	 * that are located in between @a base_position[lower_index] and @a base_position[lower_index+1]
	 * @param[in] lower_index index for @a base_position. Interpolated area is defined as
	 * @a base_position[lower_index] and @a base_position[lower_index+1],
	 * @param[in] work_data working data for interpolation
	 * @param[out] interpolated_data storage for interpolation result. Its length must be
	 * @a num_interpolated. must-be-aligned
	 */
	template<class Indexer>
	static void DoInterpolate(size_t num_base, XDataType const base_position[],
			YDataType const base_data[], size_t num_interpolated,
			XDataType const interpolated_position[], size_t index_from,
			size_t index_to, size_t lower_index,
			WorkingData const * const work_data, size_t num_array,
			size_t iarray, YDataType interpolated_data[]) {
		assert(num_base >= work_data->num_elements);
		size_t const left_edge2 = num_base - work_data->num_elements;
		size_t left_edge =
				(lower_index >= work_data->polynomial_order / 2u) ?
						lower_index - work_data->polynomial_order / 2u : 0;
		left_edge = (left_edge > left_edge2) ? left_edge2 : left_edge;
		for (size_t i = index_from; i < index_to; ++i) {
			interpolated_data[Indexer::GetIndex(num_interpolated, num_array,
					iarray, i)] = PerformNevilleAlgorithm(num_base,
					base_position, base_data, left_edge,
					interpolated_position[i], work_data);
		}
	}
private:
	static YDataType PerformNevilleAlgorithm(size_t num_base,
			XDataType const base_position[], YDataType const base_data[],
			size_t lower_index, XDataType interpolated_position,
			WorkingData const * const work_data) {
		// working pointers
		auto const x_ptr = &(base_position[lower_index]);
		auto const y_ptr = &(base_data[lower_index]);

		auto c = work_data->xholder[0].pointer;
		auto d = work_data->xholder[1].pointer;

		for (size_t i = 0; i < work_data->num_elements; ++i) {
			c[i] = WDataType(y_ptr[i]);
			d[i] = c[i];
		}

		// Neville's algorithm
		auto work = c[0];
		for (size_t m = 1; m < work_data->num_elements; ++m) {
			// Evaluate Cm1, Cm2, Cm3, ... Cm[n-m] and Dm1, Dm2, Dm3, ... Dm[n-m].
			// Those are stored to c[0], c[1], ..., c[n-m-1] and d[0], d[1], ...,
			// d[n-m-1].
			for (size_t i = 0; i < work_data->num_elements - m; ++i) {
				auto cd = c[i + 1] - d[i];
				auto const dx = WDataType(x_ptr[i] - x_ptr[i + m]);
				assert(dx != 0);
				cd /= dx;
				c[i] = WDataType(x_ptr[i] - interpolated_position) * cd;
				d[i] = WDataType(x_ptr[i + m] - interpolated_position) * cd;
			}

			// In each step, c[0] holds Cm1 which is a correction between
			// P12...m and P12...[m+1]. Thus, the following repeated update
			// corresponds to the route P1 -> P12 -> P123 ->...-> P123...n.
			work += c[0];
		}

		return static_cast<YDataType>(work);
	}
};

/**
 * Spline interpolation implementation.
 *
 * @tparam XDataType data type for position
 * @tparam YDataType data type for data
 */
template<class XDataType, class YDataType>
struct SplineInterpolatorImpl: public InterpolatorInterface<
		SplineInterpolatorImpl<XDataType, YDataType>,
		SplineWorkingData<XDataType, YDataType>, XDataType, YDataType> {
	typedef typename InterpolatorInterface<
			SplineInterpolatorImpl<XDataType, YDataType>,
			SplineWorkingData<XDataType, YDataType>, XDataType, YDataType>::WorkingData WorkingData;
	typedef typename std::common_type<XDataType, YDataType>::type WDataType;
	/**
	 * Perform spline interpolation.
	 * For each @a interpolated_position that locates in between @a base_position[lower_index]
	 * and @a base_position[lower_index+1], perform linear interpolation and set the result
	 * to @a interpolated_data. Interpolated curve is a spline that goes through all data
	 * points that are defined by @a base_position and @a base_data. Spline curve is evaluated
	 * by using straight line that goes through (@a base_position[lower_index], @a base_data[lower_index])
	 * and (@a base_position[lower_index+1], @a base_data[lower_index+1]) and spline correction term
	 * that is stored in @a work_data.
	 *
	 * @param[in] num_base number of elements for data points.
	 * @param[in] base_position position of data points. Its length must be @a num_base.
	 * must-be-aligned
	 * @param[in] base_data value of data points. Its length must be @a num_base.
	 * must-be-aligned
	 * @param[in] num_interpolated number of elements for points that wants to get
	 * interpolated value.
	 * @param[in] interpolated_position location of points that wants to get interpolated
	 * value. Its length must be @a num_interpolated.
	 * @param[in] index_from start index for @a interpolated_position and @a interpolated_data
	 * that are located in between @a base_position[lower_index] and @a base_position[lower_index+1]
	 * @param[in] index_to end index for @a interpolated_position and @a interpolated_data
	 * that are located in between @a base_position[lower_index] and @a base_position[lower_index+1]
	 * @param[in] lower_index index for @a base_position. Interpolated area is defined as
	 * @a base_position[lower_index] and @a base_position[lower_index+1],
	 * @param[in] work_data working data for interpolation
	 * @param[out] interpolated_data storage for interpolation result. Its length must be
	 * @a num_interpolated. must-be-aligned
	 */
	template<class Indexer>
	static void DoInterpolate(size_t num_base, XDataType const base_position[],
			YDataType const base_data[], size_t num_interpolated,
			XDataType const interpolated_position[], size_t index_from,
			size_t index_to, size_t lower_index,
			WorkingData const * const work_data, size_t num_array,
			size_t iarray, YDataType interpolated_data[]) {
		auto const d2ydx2(work_data->second_derivative.pointer);
		auto const dx = WDataType(
				base_position[lower_index + 1] - base_position[lower_index]);
		auto const dx_factor = dx * dx / decltype(dx)(6.0);
		auto const position_upper = WDataType(base_position[lower_index + 1]);
		auto const data_lower = WDataType(base_data[lower_index]);
		auto const data_upper = WDataType(base_data[lower_index + 1]);
		auto const d2ydx2_lower(d2ydx2[lower_index]);
		auto const d2ydx2_upper(d2ydx2[lower_index + 1]);
		for (size_t i = index_from; i < index_to; ++i) {
			auto const a =
					(position_upper - WDataType(interpolated_position[i])) / dx;
			auto const b = decltype(a)(1.0) - a;
			interpolated_data[Indexer::GetIndex(num_interpolated, num_array,
					iarray, i)] = static_cast<YDataType>(a * data_lower
					+ b * data_upper
					+ ((a * a * a - a) * d2ydx2_lower
							+ (b * b * b - b) * d2ydx2_upper) * dx_factor);
		}
	}
};

/**
 * Pick up valid (mask is true) data and fill out output array with them.
 * If the data is sorted in descending order, the array order will be
 * reversed.
 *
 * @tparam Indexer data accessor
 * @tparam XDataType data type for position
 * @tparam YDataType data type for data
 *
 * @param[in] num_column number of elements for @a position. It also defines a
 * number of columns for the @a data.
 * @param[in] position list of data positions. Length must be @a num_column.
 * must-be-aligned
 * @param[in] num_row number of rows for the @a data.
 * @param[in] the_row target row index for the @a data. It must be less than
 * @a num_row.
 * @param[in] data data array. Its length must be @a num_column times @a num_row.
 * must-be-aligned
 * @param[in] mask data mask. true is valid. Its length must be @a num_column
 * times @a num_row. must-be-aligned
 * @param[out] out_position output valid position list. Its length must be
 * @a num_column. Number of valid data is returned value of this function.
 * @param[out] out_data output valid data list. Its length must be
 * @a num_column. Number of valid data is returned value of this function.
 * @return number of valid data
 */
template<class Indexer, class XDataType, class YDataType>
inline size_t PickValidDataAsAscending(size_t num_column,
		XDataType const position[], size_t num_row, size_t the_row,
		YDataType const data[],
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

/**
 * Copy contents of working data specified as @a in_data into @a out_data.
 * Copied location is defined by @a Indexer.
 *
 * @tparam Indexer data accessor
 * @tparam DataType data type for data
 *
 * @param[in] num_column number of elements for @a in_data. It also defines a
 * number of columns for the @a out_data.
 * @param[in] num_row number of rows for the @a out_data.
 * @param[in] the_row target row index for the @a out_data. It must be less than
 * @a num_row.
 * @param[in] in_data input working data. Its length must be @a num_column.
 * must-be-aligned
 * @param[out] out_data output data. Its length must be @a num_column times
 * @a num_row. must-be-aligned
 */
template<class Indexer, class DataType>
inline void CopyResult(size_t num_column, DataType in_data[], size_t num_row,
		size_t the_row, DataType out_data[]) {
	for (size_t i = 0; i < num_column; ++i) {
		size_t index = Indexer::GetIndex(num_column, num_row, the_row, i);
		assert(index < num_row * num_column);
		out_data[index] = in_data[i];
	}
}

/**
 * It fills specified ranges of output data with specified value.
 *
 * @tparam DataType data type for data
 *
 * @param[in] the_value the value to be filled.
 * @param[in] start start index for @a out.
 * @param[in] end end index for @a out.
 * @param[out] out output data array.
 */
template<class Indexer, class DataType>
inline void FillOutOfRangeAreaWithValue(DataType the_value, size_t start,
		size_t end, size_t num_column, size_t num_row, size_t the_row,
		DataType out[]) {
	for (size_t i = start; i < end; ++i) {
		out[Indexer::GetIndex(num_column, num_row, the_row, i)] = the_value;
	}
}

/**
 * It fills one row with specified value. Data accessing rule is defined
 * by @a Indexer.
 *
 * @tparam Indexer data accessor
 * @tparam DataType data type for data
 *
 * @param[in] num_column number of columns.
 * @param[in] num_row number of rows.
 * @param[in] value the value to be filled.
 * @param[in] the_row target row index. It must be less than @a num_row
 * @param[out] data data array. Its length must be @a num_column times
 * @a num_row. must-be-aligned
 */
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

/**
 * Array accessor for X-axis interpolation with ascending array.
 */
struct XAscendingIndexer {
	static size_t GetIndex(size_t num_column, size_t num_row, size_t the_row,
			size_t index) {
		assert(the_row < num_row);
		assert(index < num_column);
		return num_column * the_row + index;
	}
};

/**
 * Array accessor for X-axis interpolation with descending array.
 */
struct XDescendingIndexer {
	static size_t GetIndex(size_t num_column, size_t num_row, size_t the_row,
			size_t index) {
		assert(the_row < num_row);
		assert(index < num_column);
		return num_column * the_row + num_column - 1 - index;
	}
};

/**
 * Array accessor for Y-axis interpolation with ascending array.
 */
struct YAscendingIndexer {
	static size_t GetIndex(size_t num_column, size_t num_row, size_t the_row,
			size_t index) {
		assert(the_row < num_row);
		assert(index < num_column);
		return the_row + num_row * index;
	}
};

/**
 * Array accessor for Y-axis interpolation with descending array.
 */
struct YDescendingIndexer {
	static size_t GetIndex(size_t num_column, size_t num_row, size_t the_row,
			size_t index) {
		assert(the_row < num_row);
		assert(index < num_column);
		return the_row + num_row * (num_column - 1 - index);
	}
};

/**
 * Helper class for X-axis interpolation.
 */
struct XInterpolatorHelper {
	typedef XAscendingIndexer AscendingIndexer;
	typedef XDescendingIndexer DescendingIndexer;
};

/**
 * Helper class for Y-axis interpolation.
 */
struct YInterpolatorHelper {
	typedef YAscendingIndexer AscendingIndexer;
	typedef YDescendingIndexer DescendingIndexer;
};

/**
 * Core function for interpolation. The procedure is as follows:
 *
 *     -# Check sort order (ascending or descending)
 *     -# Initialization
 *         - Initialize helper function
 *         - Initialize output mask to true
 *         - Prepare working arrays
 *     -# Instantiate working data
 *     -# Main loop on @a num_array:
 *          -# Pick up valid data from @a base_position and @a base_data
 *          -# If there is no valid data, set all @a interpolated_mask elements to false
 *          -# If there is only one valid data, fill @a interpolated_data with the value
 *          -# Otherwise, perform interpolation
 *
 *  Interpolation is done by splitting @a interpolated_data into three regions based on
 *  @a interpolated_position: 1) @a interpolated_position is less than @a base_position[0],
 *  2) @a interpolated_position is in between @a base_position[0] and @a base_position[num_base-1],
 *  and 3) @a interpolated_position is greater than @a base_position[num_base-1].
 *  In regions 1) and 3), no interpolation is done and @a interpolated_data is filled using
 *  nearest value (either @a base_data[0] or @a base_data[num_base-1]). In region 2),
 *  interpolation is performed using @a Interpolator.
 *
 *  See @link sakura_InterpolateXAxisFloat sakura_InterpolateXAxisFloat @endlink for
 *  parameter description.
 *
 * @tparam Interpolator interpolation engine class
 * @tparam Helper class for interpolation
 * @tparam XDataType data type for position
 * @tparam YDataType data type for data
 *
 * @param[in] polynomial_order
 * @param[in] num_base
 * @param[in] base_position
 * @param[in] num_array
 * @param[in] base_data
 * @param[in] base_mask
 * @param[in] num_interpolated
 * @param[in] interpolated_position
 * @param[out] interpolated_data
 * @param[out] interpolated_mask
 */
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
	auto data_picker =
			(is_base_ascending) ?
					PickValidDataAsAscending<AscendingIndexer, XDataType,
							YDataType> :
					PickValidDataAsAscending<DescendingIndexer, XDataType,
							YDataType>;
	auto interpolator =
			(is_interp_ascending) ?
					Interpolator::template Interpolate1D<AscendingIndexer> :
					Interpolator::template Interpolate1D<DescendingIndexer>;
	auto data_filler =
			(is_interp_ascending) ?
					FillOutOfRangeAreaWithValue<AscendingIndexer, YDataType> :
					FillOutOfRangeAreaWithValue<DescendingIndexer, YDataType>;

	// Initialize interpolated_mask to true
	for (size_t i = 0; i < num_interpolated * num_array; ++i) {
		interpolated_mask[i] = true;
	}

	// Working arrays
	StorageAndAlignedPointer<XDataType> x0_storage, x1_storage;
	StorageAndAlignedPointer<YDataType> y0_storage, y1_storage;
	AllocateAndAlign(num_base, &x0_storage);
	AllocateAndAlign(num_base, &y0_storage);
	XDataType *valid_base_position = x0_storage.pointer;
	YDataType *valid_base_data = y0_storage.pointer;
	XDataType const *interpolated_position_ascending = interpolated_position;
	if (!is_interp_ascending) {
		AllocateAndAlign(num_interpolated, &x1_storage);
		XDataType *x_writable = x1_storage.pointer;
		for (size_t i = 0; i < num_interpolated; ++i) {
			x_writable[i] = interpolated_position[num_interpolated - 1 - i];
		}
		interpolated_position_ascending = x_writable;
	}

	std::unique_ptr<WorkingData> wdata_storage(
			WorkingData::Allocate(polynomial_order, num_base));
	WorkingData *work_data = wdata_storage.get();

	StorageAndAlignedPointer<size_t> size_t_holder[2];

	for (size_t iarray = 0; iarray < num_array; ++iarray) {
		// Pick up valid data and associating position from base_position and base_data
		// by referring base_mask. The data and position are reversed if sort order is
		// descending.
		size_t const num_valid_data = data_picker(num_base, base_position,
				num_array, iarray, base_data, base_mask, valid_base_position,
				valid_base_data);
		if (num_valid_data == 0) {
			// cannot perform interpolation, mask all data
			FillOneRowWithValue<AscendingIndexer, bool>(num_interpolated,
					num_array, false, iarray, interpolated_mask);
		} else if (num_valid_data == 1) {
			// no need to interpolate, just substitute single base data
			// to all elements in interpolated data, keep input mask
			FillOneRowWithValue<AscendingIndexer, YDataType>(num_interpolated,
					num_array, valid_base_data[0], iarray, interpolated_data);
		} else {
			// perform interpolation, keep input mask

			// Locate each element in base_position against interpolated_position
			if (size_t_holder[0].pointer == nullptr) {
				AllocateAndAlign<size_t>(num_base, &size_t_holder[0]);
				AllocateAndAlign<size_t>(num_base, &size_t_holder[1]);
			}
			size_t *location_base = size_t_holder[0].pointer;
			size_t *lower_index_base = size_t_holder[1].pointer;
			size_t const num_location_base = Locate<XDataType>(num_interpolated,
					interpolated_position_ascending, num_valid_data,
					valid_base_position, location_base, lower_index_base);

			// interpolated_position is less than base_position[0]
			data_filler(valid_base_data[0], 0, location_base[0],
					num_interpolated, num_array, iarray, interpolated_data);

			// Perform 1-dimensional interpolation
			// interpolated_position is located in between base_position[0]
			// and base_position[num_base-1]
			WorkingData::Initialize(num_valid_data, valid_base_position,
					valid_base_data, work_data);
			interpolator(num_valid_data, valid_base_position, valid_base_data,
					num_interpolated, interpolated_position_ascending,
					num_location_base, location_base, lower_index_base,
					work_data, num_array, iarray, interpolated_data);

			// interpolated_position is greater than base_position[num_base-1]
			data_filler(valid_base_data[num_valid_data - 1],
					location_base[num_location_base - 1], num_interpolated,
					num_interpolated, num_array, iarray, interpolated_data);
		}
	}
}

#ifndef NDEBUG
/**
 * Template utility function for sort check
 *
 * @tparam DataType data type for data
 * @tparam Checker data checker. It must take two DataType arguments and
 * return boolean value.
 *
 * @param[in] checker function to check array elements
 * @param[in] num_elements Length of the @a data
 * @param[in] data data to be inspected. Length must be @a num_elements.
 * must-be-aligned
 * @return true if input array is sorted and no duplicates, otherwise false
 */
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

/**
 * Check if input array is sorted and no duplicates.
 *
 * @tparam DataType data type for data
 *
 * @param[in] num_elements Length of the @a data
 * @param[out] data data to be inspected. Length must be @a num_elements.
 * must-be-aligned
 * @return true if input array is sorted and no duplicates, otherwise false
 */
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

/**
 * Perform basic check of input arguments.
 */
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

/**
 * @a DoInterpolate calls core function @a Interpolate1D with appropriate
 * template arguments according to user-specified @a interpolation method and
 * @a polynomial_order. If any input argument is invalid, it doesn't execute
 * @a Interpolate1D and appropriate error status is returned. Interpolator
 * class are set as the following manner:
 *
 *     - @a NearestInterpolator if interpolation_method is @a sakura_InterpolationMethod_kNearest
 *     - @a LinearInterpolator if interpolation_method is @a sakura_InterpolationMethod_kLinear
 *     - @a SplineInterpolator if interpolation_method is @a sakura_InterpolationMethod_kSpline
 *     - @a PolynomialInterpolator if interpolation_method is @a sakura_InterpolationMethod_kPolynomial
 *       and @a polynomial_order is not 0
 *     - @a NearestInterpolator if interpolation_method is @a sakura_InterpolationMethod_kPolynomial
 *       and @a polynomial_order is 0
 *
 * @tparam XDataType data type for position
 * @tparam YDataType data type for data
 * @tparam InterpolatorHelper Helper class for interpolation
 */
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
 * where N and M correspond to num_base and num_array respectively.
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
