#include <cassert>
#include <sstream>
#include <memory>
#include <cstdalign>
#include <utility>
#include <vector>

#include <libsakura/sakura.h>
#include <libsakura/optimized_implementation_factory_impl.h>
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
	~StorageAndAlignedPointer() {
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

template<class DataType>
class XAxisReordererImpl {
public:
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
};

template<class DataType>
class YAxisReordererImpl {
public:
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
};

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

	assert(LIBSAKURA_SYMBOL(IsAligned)(base_position));

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
			size_t middle_index = (right_index + left_index) / 2;
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
	} else {
		size_t start_position = 0;
		size_t end_position = num_base - 1;
		size_t previous_location = num_base + 1;
		for (size_t i = 0; i < num_located; ++i) {
			size_t location = LocateData(start_position, end_position, num_base,
					base_position, located_position[i]);
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

template<class XDataType, class YDataType>
class NearestXInterpolatorImpl {
public:
	void PrepareForInterpolation(uint8_t polynomial_order, size_t num_base,
			size_t num_array, XDataType const base_position[],
			YDataType const base_data[]) {
	}
	void Interpolate1D(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t num_location,
			size_t const location[]) {
		for (size_t k = 0; k < num_location - 1; ++k) {
			XDataType middle_point = 0.5
					* (base_position[k + 1] + base_position[k]);
			for (size_t j = 0; j < num_array; ++j) {
				YDataType left_value = base_data[j * num_base + k];
				YDataType right_value = base_data[j * num_base + k + 1];
				YDataType *work = &interpolated_data[j * num_interpolated];
				for (size_t i = location[k]; i < location[k + 1]; ++i) {
					// nearest condition
					// - interpolated_position[i] == middle_point
					//     ---> nearest is y_left (left side)
					// - interpolated_position[i] < middle_point and ascending order
					//     ---> nearest is y_left (left side)
					// - interpolated_position[i] > middle_point and ascending order
					//     ---> nearest is y_right (right side)
					// - interpolated_position[i] < middle_point and descending order
					//     ---> nearest is y_right (right side)
					// - interpolated_position[i] > middle_point and descending order
					//     ---> nearest is y_left (left side)
					if ((interpolated_position[i] != middle_point)
							& ((interpolated_position[i] > middle_point))) {
						work[i] = right_value;
					} else {
						work[i] = left_value;
					}
				}
			}
		}
	}
};

template<class XDataType, class YDataType>
class LinearXInterpolatorImpl {
public:
	void PrepareForInterpolation(uint8_t polynomial_order, size_t num_base,
			size_t num_array, XDataType const base_position[],
			YDataType const base_data[]) {
	}
	void Interpolate1D(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t num_location,
			size_t const location[]) {
		for (size_t k = 0; k < num_location - 1; ++k) {
			for (size_t j = 0; j < num_array; ++j) {
				size_t offset_index_left = j * num_base + k;
				XDataType dydx =
						static_cast<XDataType>(base_data[offset_index_left + 1]
								- base_data[offset_index_left])
								/ (base_position[k + 1] - base_position[k]);
				YDataType base_term = base_data[offset_index_left];
				YDataType *y_work = &interpolated_data[j * num_interpolated];
				for (size_t i = location[k]; i < location[k + 1]; ++i) {
					y_work[i] = base_term
							+ static_cast<YDataType>(dydx
									* (interpolated_position[i]
											- base_position[k]));
				}
			}
		}
	}
};

template<class XDataType, class YDataType>
class PolynomialXInterpolatorImpl {
public:
	PolynomialXInterpolatorImpl() :
			polynomial_order_(0), num_elements_(0), holder_(2) {
	}
	void PrepareForInterpolation(uint8_t polynomial_order, size_t num_base,
			size_t num_array, XDataType const base_position[],
			YDataType const base_data[]) {
		polynomial_order_ =
				(polynomial_order + 1u >= num_base) ?
						static_cast<uint8_t>(num_base - 1) : polynomial_order;
		num_elements_ = polynomial_order_ + 1;
		AllocateAndAlign<XDataType>(num_elements_, &holder_[0]);
		AllocateAndAlign<XDataType>(num_elements_, &holder_[1]);
	}
	void Interpolate1D(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t num_location,
			size_t const location[]) {
		for (size_t k = 0; k < num_location - 1; ++k) {
			int left_edge1 = k - polynomial_order_ / 2;
			size_t left_edge2 = num_base - num_elements_;
			size_t left_edge = static_cast<size_t>(
					(left_edge1 > 0) ? left_edge1 : 0);
			left_edge = (left_edge > left_edge2) ? left_edge2 : left_edge;
			for (size_t j = 0; j < num_array; ++j) {
				for (size_t i = location[k]; i < location[k + 1]; ++i) {
					PerformNevilleAlgorithm(num_base, base_position, base_data,
							left_edge, j, interpolated_position[i],
							&interpolated_data[j * num_interpolated + i]);
				}
			}
		}
	}
private:
	void PerformNevilleAlgorithm(size_t num_base,
			XDataType const base_position[], YDataType const base_data[],
			size_t left_index, size_t array_index,
			XDataType interpolated_position, YDataType *interpolated_data) {

		// working pointers
		XDataType const *x_ptr = &(base_position[left_index]);
		YDataType const *y_ptr =
				&(base_data[num_base * array_index + left_index]);

		XDataType *c = holder_[0].pointer;
		XDataType *d = holder_[1].pointer;

		for (size_t i = 0; i < num_elements_; ++i) {
			c[i] = static_cast<XDataType>(y_ptr[i]);
			d[i] = c[i];
		}

		// Neville's algorithm
		XDataType work = c[0];
		for (size_t m = 1; m < num_elements_; ++m) {
			// Evaluate Cm1, Cm2, Cm3, ... Cm[n-m] and Dm1, Dm2, Dm3, ... Dm[n-m].
			// Those are stored to c[0], c[1], ..., c[n-m-1] and d[0], d[1], ...,
			// d[n-m-1].
			for (size_t i = 0; i < num_elements_ - m; ++i) {
				XDataType cd = c[i + 1] - d[i];
				XDataType dx = x_ptr[i] - x_ptr[i + m];
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

		*interpolated_data = static_cast<YDataType>(work);
	}

	uint8_t polynomial_order_;
	size_t num_elements_;
	std::vector<StorageAndAlignedPointer<XDataType> > holder_;
};

template<class XDataType, class YDataType>
class SplineXInterpolatorImpl {
public:
	SplineXInterpolatorImpl() :
			holder_(2) {
	}
	void PrepareForInterpolation(uint8_t polynomial_order, size_t num_base,
			size_t num_array, XDataType const base_position[],
			YDataType const base_data[]) {
		// Derive second derivative at x_base
		AllocateAndAlign<YDataType>(num_base * num_array, &holder_[0]);
		AllocateAndAlign<YDataType>(num_base, &holder_[1]);
		for (size_t i = 0; i < num_array; ++i) {
			size_t start_position = i * num_base;
			YDataType const *base_data_work = &(base_data[start_position]);
			YDataType *d2ydx2 = &(holder_[0].pointer[start_position]);
			DeriveSplineCorrectionTerm(num_base, base_position, base_data_work,
					d2ydx2);
		}
	}
	void Interpolate1D(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t num_location,
			size_t const location[]) {
		YDataType *d2ydx2 = holder_[0].pointer;
		assert(d2ydx2 != nullptr);
		for (size_t k = 0; k < num_location - 1; ++k) {
			XDataType dx = base_position[k + 1] - base_position[k];
			XDataType dx_factor = dx * dx / 6.0;
			for (size_t j = 0; j < num_array; ++j) {
				size_t offset_index_left = j * num_base + k;
				YDataType *work = &interpolated_data[j * num_interpolated];
				for (size_t i = location[k]; i < location[k + 1]; ++i) {
					XDataType a = (base_position[k + 1]
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
	}
private:
	void DeriveSplineCorrectionTerm(size_t num_base,
			XDataType const base_position[], YDataType const base_data[],
			YDataType d2ydx2[]) {
		YDataType *upper_triangular = holder_[1].pointer;

		// This is a condition of natural cubic spline
		d2ydx2[0] = 0.0;
		d2ydx2[num_base - 1] = 0.0;
		upper_triangular[0] = 0.0;

		// Solve tridiagonal system.
		// Here tridiagonal matrix is decomposed to upper triangular matrix.
		// upper_tridiangular stores upper triangular elements, while
		// d2ydx2 stores right-hand-side vector. The diagonal
		// elements are normalized to 1.

		// x_base is ascending order
		XDataType a1 = base_position[1] - base_position[0];
		for (size_t i = 1; i < num_base - 1; ++i) {
			XDataType a2 = base_position[i + 1] - base_position[i];
			XDataType b1 = 1.0 / (base_position[i + 1] - base_position[i - 1]);
			d2ydx2[i] = 3.0 * b1
					* ((base_data[i + 1] - base_data[i]) / a2
							- (base_data[i] - base_data[i - 1]) / a1
							- d2ydx2[i - 1] * 0.5 * a1);
			a1 = 1.0 / (1.0 - upper_triangular[i - 1] * 0.5 * a1 * b1);
			d2ydx2[i] *= a1;
			upper_triangular[i] = 0.5 * a2 * b1 * a1;
			a1 = a2;
		}

		// Solve the system by backsubstitution and store solution to d2ydx2_
		for (size_t k = num_base - 2; k >= 1; --k) {
			d2ydx2[k] -= upper_triangular[k] * d2ydx2[k + 1];
		}
	}
	std::vector<StorageAndAlignedPointer<YDataType> > holder_;
};

template<class XDataType, class YDataType>
class NearestYInterpolatorImpl {
public:
	void PrepareForInterpolation(uint8_t polynomial_order, size_t num_base,
			size_t num_array, XDataType const base_position[],
			YDataType const base_data[]) {
	}
	void Interpolate1D(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t num_location,
			size_t const location[]) {
		for (size_t k = 0; k < num_location - 1; ++k) {
			XDataType middle_point = 0.5
					* (base_position[k + 1] + base_position[k]);
			for (size_t i = location[k]; i < location[k + 1]; ++i) {
				// nearest condition
				// - interpolated_position[i] == middle_point
				//     ---> nearest is y_left (left side)
				// - interpolated_position[i] < middle_point and ascending order
				//     ---> nearest is y_left (left side)
				// - interpolated_position[i] > middle_point and ascending order
				//     ---> nearest is y_right (right side)
				// - interpolated_position[i] < middle_point and descending order
				//     ---> nearest is y_right (right side)
				// - interpolated_position[i] > middle_point and descending order
				//     ---> nearest is y_left (left side)
				size_t offset_index = 0;
				if ((interpolated_position[i] != middle_point)
						& ((interpolated_position[i] > middle_point))) {
					offset_index = 1;
				}
				YDataType *work = &interpolated_data[num_array * i];
				YDataType const *nearest = &base_data[num_array
						* (k + offset_index)];
				for (size_t j = 0; j < num_array; ++j) {
					work[j] = nearest[j];
				}
			}
		}
	}
};

template<class XDataType, class YDataType>
class LinearYInterpolatorImpl {
public:
	void PrepareForInterpolation(uint8_t polynomial_order, size_t num_base,
			size_t num_array, XDataType const base_position[],
			YDataType const base_data[]) {
	}
	void Interpolate1D(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t num_location,
			size_t const location[]) {
		for (size_t k = 0; k < num_location - 1; ++k) {
			for (size_t i = location[k]; i < location[k + 1]; ++i) {
				YDataType fraction =
						static_cast<YDataType>((interpolated_position[i]
								- base_position[k])
								/ (base_position[k + 1] - base_position[k]));
				YDataType *work = &interpolated_data[num_array * i];
				YDataType const *left_value = &base_data[num_array * k];
				YDataType const *right_value = &base_data[num_array * (k + 1)];
				for (size_t j = 0; j < num_array; ++j) {
					work[j] = left_value[j]
							+ fraction * (right_value[j] - left_value[j]);
				}
			}
		}
	}
};

template<class XDataType, class YDataType>
class PolynomialYInterpolatorImpl {
public:
	PolynomialYInterpolatorImpl() :
			polynomial_order_(0), num_elements_(0), holder_(3) {
	}
	void PrepareForInterpolation(uint8_t polynomial_order, size_t num_base,
			size_t num_array, XDataType const base_position[],
			YDataType const base_data[]) {
		polynomial_order_ =
				(polynomial_order + 1u >= num_base) ?
						static_cast<uint8_t>(num_base - 1) : polynomial_order;
		num_elements_ = polynomial_order_ + 1;
		AllocateAndAlign(num_elements_ * num_array, &holder_[0]);
		AllocateAndAlign(num_elements_ * num_array, &holder_[1]);
		AllocateAndAlign(num_array, &holder_[2]);
	}
	void Interpolate1D(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t num_location,
			size_t const location[]) {
		for (size_t k = 0; k < num_location - 1; ++k) {
			int left_edge1 = k - polynomial_order_ / 2;
			size_t left_edge2 = num_base - num_elements_;
			size_t left_edge = static_cast<size_t>(
					(left_edge1 > 0) ? left_edge1 : 0);
			left_edge = (left_edge > left_edge2) ? left_edge2 : left_edge;
			for (size_t i = location[k]; i < location[k + 1]; ++i) {
				PerformNevilleAlgorithm(num_array, base_position, base_data,
						left_edge, interpolated_position[i],
						&interpolated_data[num_array * i]);
			}
		}
	}
private:
	void PerformNevilleAlgorithm(size_t num_array,
			XDataType const base_position[], YDataType const base_data[],
			size_t left_index, XDataType interpolated_position,
			YDataType interpolated_data[]) {
		// working pointers
		XDataType const *x_ptr = &(base_position[left_index]);
		XDataType *c = holder_[0].pointer;
		XDataType *d = holder_[1].pointer;
		XDataType *work = holder_[2].pointer;

		size_t start = left_index * num_array;
		size_t num_elements = num_elements_ * num_array;
		for (size_t i = 0; i < num_elements; ++i) {
			c[i] = static_cast<XDataType>(base_data[start + i]);
			d[i] = c[i];
		}

		// Neville's algorithm
		for (size_t i = 0; i < num_array; ++i) {
			work[i] = c[i];
		}
		for (size_t m = 1; m < num_elements_; ++m) {
			// Evaluate Cm1, Cm2, Cm3, ... Cm[n-m] and Dm1, Dm2, Dm3, ... Dm[n-m].
			// Those are stored to c[0], c[1], ..., c[n-m-1] and d[0], d[1], ...,
			// d[n-m-1].
			for (size_t i = 0; i < num_elements_ - m; ++i) {
				XDataType dx = x_ptr[i] - x_ptr[i + m];
				assert(dx != 0);
				size_t offset = i * num_array;
				for (size_t j = 0; j < num_array; ++j) {
					XDataType cd = (c[offset + num_array + j] - d[offset + j])
							/ dx;
					c[offset + j] = (x_ptr[i] - interpolated_position) * cd;
					d[offset + j] = (x_ptr[i + m] - interpolated_position) * cd;
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

	uint8_t polynomial_order_;
	size_t num_elements_;
	std::vector<StorageAndAlignedPointer<XDataType> > holder_;
};

template<class XDataType, class YDataType>
class SplineYInterpolatorImpl {
public:
	void PrepareForInterpolation(uint8_t polynomial_order, size_t num_base,
			size_t num_array, XDataType const base_position[],
			YDataType const base_data[]) {
		// Derive second derivative at x_base
		AllocateAndAlign<YDataType>(num_base * num_array, &holder_);
		DeriveSplineCorrectionTerm(num_base, base_position, num_array,
				base_data, holder_.pointer);
	}
	void Interpolate1D(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t num_location,
			size_t const location[]) {
		YDataType *d2ydx2 = holder_.pointer;
		assert(d2ydx2 != nullptr);
		for (size_t k = 0; k < num_location - 1; ++k) {
			XDataType dx = base_position[k + 1] - base_position[k];
			YDataType const *left_value = &base_data[k * num_array];
			YDataType const *right_value = &base_data[(k + 1) * num_array];
			YDataType const *d2ydx2_left = &d2ydx2[k * num_array];
			YDataType const *d2ydx2_right = &d2ydx2[(k + 1) * num_array];
			for (size_t i = location[k]; i < location[k + 1]; ++i) {
				XDataType a = (base_position[k + 1] - interpolated_position[i])
						/ dx;
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
	}
private:
	void DeriveSplineCorrectionTerm(size_t num_base,
			XDataType const base_position[], size_t num_array,
			YDataType const base_data[], YDataType d2ydx2[]) {
		StorageAndAlignedPointer<YDataType> holder_for_u;
		AllocateAndAlign<YDataType>(num_base * num_array, &holder_for_u);
		YDataType *upper_triangular = holder_for_u.pointer;

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
		for (size_t i = 1; i < num_base - 1; ++i) {
			XDataType a2 = base_position[i + 1] - base_position[i];
			XDataType b1 = 1.0 / (base_position[i + 1] - base_position[i - 1]);
			for (size_t j = 0; j < num_array; ++j) {
				d2ydx2[num_array * i + j] = 3.0 * b1
						* ((base_data[num_array * (i + 1) + j]
								- base_data[num_array * i + j]) / a2
								- (base_data[num_array * i + j]
										- base_data[num_array * (i - 1) + j])
										/ a1
								- d2ydx2[num_array * (i - 1) + j] * 0.5 * a1);
				XDataType a3 = 1.0
						/ (1.0
								- upper_triangular[num_array * (i - 1) + j]
										* 0.5 * a1 * b1);
				d2ydx2[num_array * i + j] *= a3;
				upper_triangular[num_array * i + j] = 0.5 * a2 * b1 * a3;
			}
			a1 = a2;
		}

		// Solve the system by backsubstitution and store solution to d2ydx2_
		for (size_t k = num_base - 2; k >= 1; --k) {
			for (size_t j = 0; j < num_array; ++j) {
				d2ydx2[k * num_array + j] -= upper_triangular[k * num_array + j]
						* d2ydx2[(k + 1) * num_array + j];
			}
		}
	}
	StorageAndAlignedPointer<YDataType> holder_;
};

template<class InterpolatorImpl, class XDataType, class YDataType>
struct XInterpolator: public InterpolatorImpl {
	typedef XAxisReordererImpl<XDataType> XDataReorderer;
	typedef XAxisReordererImpl<YDataType> YDataReorderer;
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

template<class InterpolatorImpl, class XDataType, class YDataType>
struct YInterpolator: public InterpolatorImpl {
	typedef YAxisReordererImpl<XDataType> XDataReorderer;
	typedef YAxisReordererImpl<YDataType> YDataReorderer;
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

template<class XDataType, class YDataType>
struct XInterpolatorSet {
	typedef XInterpolator<NearestXInterpolatorImpl<XDataType, YDataType>,
			XDataType, YDataType> NearestInterpolator;
	typedef XInterpolator<LinearXInterpolatorImpl<XDataType, YDataType>,
			XDataType, YDataType> LinearInterpolator;
	typedef XInterpolator<PolynomialXInterpolatorImpl<XDataType, YDataType>,
			XDataType, YDataType> PolynomialInterpolator;
	typedef XInterpolator<SplineXInterpolatorImpl<XDataType, YDataType>,
			XDataType, YDataType> SplineInterpolator;
};

template<class XDataType, class YDataType>
struct YInterpolatorSet {
	typedef YInterpolator<NearestYInterpolatorImpl<XDataType, YDataType>,
			XDataType, YDataType> NearestInterpolator;
	typedef YInterpolator<LinearYInterpolatorImpl<XDataType, YDataType>,
			XDataType, YDataType> LinearInterpolator;
	typedef YInterpolator<PolynomialYInterpolatorImpl<XDataType, YDataType>,
			XDataType, YDataType> PolynomialInterpolator;
	typedef YInterpolator<SplineYInterpolatorImpl<XDataType, YDataType>,
			XDataType, YDataType> SplineInterpolator;
};

template<class Interpolator, class XDataType, class YDataType>
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
		Interpolator::SubstituteSingleBaseData(num_base, num_array,
				num_interpolated, base_data, interpolated_data);
		return;
	}

	std::vector<StorageAndAlignedPointer<XDataType> > xdatatype_holder(2);
	StorageAndAlignedPointer<YDataType> ydatatype_holder;
	GetAscendingArray<typename Interpolator::XDataReorderer, XDataType,
			XDataType>(num_base, base_position, 1, base_position,
			&xdatatype_holder[0]);
	GetAscendingArray<typename Interpolator::YDataReorderer, XDataType,
			YDataType>(num_base, base_position, num_array, base_data,
			&ydatatype_holder);
	GetAscendingArray<typename Interpolator::XDataReorderer, XDataType,
			XDataType>(num_interpolated, interpolated_position, 1,
			interpolated_position, &xdatatype_holder[1]);
	XDataType const *base_position_work = xdatatype_holder[0].pointer;
	YDataType const *base_data_work = ydatatype_holder.pointer;
	XDataType const *interpolated_position_work = xdatatype_holder[1].pointer;

	// Generate worker class
	Interpolator interpolator;

	// Perform 1-dimensional interpolation
	// Any preparation for interpolation should be done here
	interpolator.PrepareForInterpolation(polynomial_order, num_base, num_array,
			base_position_work, base_data_work);

	// Locate each element in x_base against x_interpolated
	StorageAndAlignedPointer<size_t> size_t_holder;
	AllocateAndAlign<size_t>(num_base, &size_t_holder);
	size_t *location_base = size_t_holder.pointer;
	size_t num_location_base = Locate<XDataType>(num_interpolated, num_base,
			interpolated_position_work, base_position_work, location_base);

	// Outside of x_base[0]
	Interpolator::SubstituteLeftMostData(location_base[0], num_base, num_array,
			num_interpolated, base_data_work, interpolated_data);

	// Between x_base[0] and x_base[num_x_base-1]
	interpolator.Interpolate1D(num_base, base_position_work, num_array,
			base_data_work, num_interpolated, interpolated_position_work,
			interpolated_data, num_location_base, location_base);

	// Outside of x_base[num_x_base-1]
	Interpolator::SubstituteRightMostData(location_base[num_location_base - 1],
			num_base, num_array, num_interpolated, base_data_work,
			interpolated_data);

	// swap output array
	if (interpolated_position[0]
			> interpolated_position[num_interpolated - 1]) {
		Interpolator::SwapResult(num_array, num_interpolated,
				interpolated_data);
	}
}

template<class Interpolator, class XDataType, class YDataType>
void ExecuteInterpolate1D(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_base,
		XDataType const base_position[], size_t num_array,
		YDataType const base_data[], size_t num_interpolated,
		XDataType const interpolated_position[],
		YDataType interpolated_data[]) {
	typedef void (*Interpolate1DFunc)(uint8_t, size_t, XDataType const *,
			size_t, YDataType const *, size_t, XDataType const *, YDataType *);
	Interpolate1DFunc func = nullptr;
	switch (interpolation_method) {
	case LIBSAKURA_SYMBOL(InterpolationMethod_kNearest):
		func = Interpolate1D<typename Interpolator::NearestInterpolator,
				XDataType, YDataType>;
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kLinear):
		func = Interpolate1D<typename Interpolator::LinearInterpolator,
				XDataType, YDataType>;
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kPolynomial):
		if (polynomial_order == 0) {
			// This is special case: 0-th polynomial interpolation
			// acts like nearest interpolation
			func = Interpolate1D<typename Interpolator::NearestInterpolator,
					XDataType, YDataType>;
		} else {
			func = Interpolate1D<typename Interpolator::PolynomialInterpolator,
					XDataType, YDataType>;
		}
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kSpline):
		func = Interpolate1D<typename Interpolator::SplineInterpolator,
				XDataType, YDataType>;
		break;
	default:
		// invalid interpolation method type
		break;
	}
	if (func != nullptr) {
		(*func)(polynomial_order, num_base, base_position, num_array, base_data,
				num_interpolated, interpolated_position, interpolated_data);
	}
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {

/**
 * Interpolate1DAlongColumn performs 1D interpolation along column based on base_position
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
template<class XDataType, class YDataType>
void ADDSUFFIX(Interpolation, ARCH_SUFFIX)<XDataType, YDataType>::InterpolateXAxis(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_x_base,
		XDataType const x_base[/*num_x_base*/], size_t num_y,
		YDataType const data_base[/*num_x_base*num_y*/],
		size_t num_x_interpolated,
		XDataType const x_interpolated[/*num_x_interpolated*/],
		YDataType data_interpolated[/*num_x_interpolated*num_y*/]) const {
	ExecuteInterpolate1D<XInterpolatorSet<XDataType, YDataType>, XDataType,
			YDataType>(interpolation_method, polynomial_order, num_x_base,
			x_base, num_y, data_base, num_x_interpolated, x_interpolated,
			data_interpolated);
}

/**
 * Interpolate1DAlongColumn performs 1D interpolation along row based on y_base
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
template<class XDataType, class YDataType>
void ADDSUFFIX(Interpolation, ARCH_SUFFIX)<XDataType, YDataType>::InterpolateYAxis(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_y_base,
		XDataType const y_base[/*num_y_base*/], size_t num_x,
		YDataType const data_base[/*num_y_base*num_x*/],
		size_t num_y_interpolated,
		XDataType const y_interpolated[/*num_y_interpolated*/],
		YDataType data_interpolated[/*num_y_interpolated*num_x*/]) const {
	ExecuteInterpolate1D<YInterpolatorSet<XDataType, YDataType>, XDataType,
			YDataType>(interpolation_method, polynomial_order, num_y_base,
			y_base, num_x, data_base, num_y_interpolated, y_interpolated,
			data_interpolated);
}

template class ADDSUFFIX(Interpolation, ARCH_SUFFIX)<double, float> ;
} /* namespace LIBSAKURA_PREFIX */
