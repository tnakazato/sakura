#include <cassert>
#include <sstream>
#include <memory>
#include <cstdalign>
#include <utility>
#include <iostream>

#include <libsakura/sakura.h>
#include <libsakura/optimized_implementation_factory_impl.h>
#include <libsakura/localdef.h>
#include <libsakura/logger.h>
#include <libsakura/memory_manager.h>

namespace {

using namespace LIBSAKURA_PREFIX;

// deleter for interpolation
struct InterpolationDeleter {
	void operator()(void *pointer) const {
		Memory::Free(pointer);
	}
};

// a logger for this module
auto logger = Logger::GetLogger("interpolation");

template<class DataType>
void *AllocateAndAlign(size_t num_array, DataType **array) {
	size_t sakura_alignment = LIBSAKURA_SYMBOL(GetAlignment)();
	size_t elements_in_arena = num_array + sakura_alignment - 1;
	void *storage = Memory::Allocate(sizeof(DataType) * elements_in_arena);
	*array = reinterpret_cast<DataType *>(LIBSAKURA_SYMBOL(AlignAny)(
			sizeof(DataType) * elements_in_arena, storage,
			sizeof(DataType) * num_array));
	return storage;
}

template<class XDataType, class YDataType>
void DeriveSplineCorrectionTerm(bool is_descending, size_t num_base,
		XDataType const x_base[], YDataType const y_base[],
		YDataType y_base_2nd_derivative[]) {
	YDataType *upper_triangular = nullptr;
	std::unique_ptr<void, InterpolationDeleter> storage_for_u(
			AllocateAndAlign<YDataType>(num_base, &upper_triangular));

	// This is a condition of natural cubic spline
	y_base_2nd_derivative[0] = 0.0;
	y_base_2nd_derivative[num_base - 1] = 0.0;
	upper_triangular[0] = 0.0;

	// Solve tridiagonal system.
	// Here tridiagonal matrix is decomposed to upper triangular matrix.
	// upper_tridiangular stores upper triangular elements, while
	// y_base_2nd_derivative stores right-hand-side vector. The diagonal
	// elements are normalized to 1.
	if (is_descending) {
		// x_base is descending order
		XDataType a1 = x_base[num_base - 2] - x_base[num_base - 1];
		for (size_t i = 1; i < num_base - 1; ++i) {
			XDataType a2 = x_base[num_base - i - 2] - x_base[num_base - i - 1];
			XDataType b1 = 1.0
					/ (x_base[num_base - i - 2] - x_base[num_base - i]);
			y_base_2nd_derivative[i] = 3.0 * b1
					* ((y_base[num_base - i - 2] - y_base[num_base - i - 1])
							/ a2
							- (y_base[num_base - i - 1] - y_base[num_base - i])
									/ a1
							- y_base_2nd_derivative[i - 1] * 0.5 * a1);
			a1 = 1.0 / (1.0 - upper_triangular[i - 1] * 0.5 * a1 * b1);
			y_base_2nd_derivative[i] *= a1;
			upper_triangular[i] = 0.5 * a2 * b1 * a1;
			a1 = a2;
		}
	} else {
		// x_base is ascending order
		XDataType a1 = x_base[1] - x_base[0];
		for (size_t i = 1; i < num_base - 1; ++i) {
			XDataType a2 = x_base[i + 1] - x_base[i];
			XDataType b1 = 1.0 / (x_base[i + 1] - x_base[i - 1]);
			y_base_2nd_derivative[i] = 3.0 * b1
					* ((y_base[i + 1] - y_base[i]) / a2
							- (y_base[i] - y_base[i - 1]) / a1
							- y_base_2nd_derivative[i - 1] * 0.5 * a1);
			a1 = 1.0 / (1.0 - upper_triangular[i - 1] * 0.5 * a1 * b1);
			y_base_2nd_derivative[i] *= a1;
			upper_triangular[i] = 0.5 * a2 * b1 * a1;
			a1 = a2;
		}
	}

	// Solve the system by backsubstitution and store solution to
	// y_base_2nd_derivative
	for (size_t k = num_base - 2; k >= 1; --k)
		y_base_2nd_derivative[k] -= upper_triangular[k]
				* y_base_2nd_derivative[k + 1];
}

template<class XDataType, class YDataType>
void DeriveSplineCorrectionTermAlongRow(size_t num_base,
		XDataType const x_base[], size_t num_array, YDataType const y_base[],
		YDataType y_base_2nd_derivative[]) {
	YDataType *upper_triangular = nullptr;
	std::unique_ptr<void, InterpolationDeleter> storage_for_u(
			AllocateAndAlign<YDataType>(num_base * num_array,
					&upper_triangular));

	// This is a condition of natural cubic spline
	for (size_t i = 0; i < num_array; ++i) {
		y_base_2nd_derivative[i] = 0.0;
		y_base_2nd_derivative[(num_base - 1) * num_array + i] = 0.0;
		upper_triangular[i] = 0.0;
	}

	// Solve tridiagonal system.
	// Here tridiagonal matrix is decomposed to upper triangular matrix.
	// upper_tridiangular stores upper triangular elements, while
	// y_base_2nd_derivative stores right-hand-side vector. The diagonal
	// elements are normalized to 1.

	// x_base is ascending order
	XDataType a1 = x_base[1] - x_base[0];
	for (size_t i = 1; i < num_base - 1; ++i) {
		XDataType a2 = x_base[i + 1] - x_base[i];
		XDataType b1 = 1.0 / (x_base[i + 1] - x_base[i - 1]);
		for (size_t j = 0; j < num_array; ++j) {
			y_base_2nd_derivative[num_array * i + j] = 3.0 * b1
					* ((y_base[num_array * (i + 1) + j]
							- y_base[num_array * i + j]) / a2
							- (y_base[num_array * i + j]
									- y_base[num_array * (i - 1) + j]) / a1
							- y_base_2nd_derivative[num_array * (i - 1) + j]
									* 0.5 * a1);
			XDataType a3 = 1.0
					/ (1.0
							- upper_triangular[num_array * (i - 1) + j] * 0.5
									* a1 * b1);
			y_base_2nd_derivative[num_array * i + j] *= a3;
			upper_triangular[num_array * i + j] = 0.5 * a2 * b1 * a3;
		}
		a1 = a2;
	}

	// Solve the system by backsubstitution and store solution to
	// y_base_2nd_derivative
	for (size_t k = num_base - 2; k >= 1; --k) {
		for (size_t j = 0; j < num_array; ++j) {
			y_base_2nd_derivative[k * num_array + j] -= upper_triangular[k
					* num_array + j]
					* y_base_2nd_derivative[(k + 1) * num_array + j];
		}
	}
}

template<class XDataType, class YDataType>
class InterpolationWorker {
public:
	InterpolationWorker(size_t num_base, XDataType const x_base[],
			size_t num_array, YDataType const y_base[],
			uint8_t polynomial_order = 0) :
			num_base_(num_base), num_array_(num_array), x_base_(x_base), y_base_(
					y_base), polynomial_order_(polynomial_order) {
	}
	virtual ~InterpolationWorker() {
	}
	virtual void PrepareForInterpolation() {
	}
	virtual void Interpolate1D(size_t left_index, size_t right_index,
			size_t location, size_t offset, XDataType const x_interpolated[],
			YDataType y_interpolated[]) = 0;
	virtual void Interpolate1DAlongRow(size_t left_index, size_t right_index,
			size_t location, XDataType const x_interpolated[],
			YDataType y_interpolated[]) = 0;
protected:
	size_t const num_base_;
	size_t const num_array_;
	XDataType const *x_base_;
	YDataType const *y_base_;
	int const polynomial_order_;
};

template<class XDataType, class YDataType>
class NearestInterpolationWorker: public InterpolationWorker<XDataType,
		YDataType> {
public:
	NearestInterpolationWorker(size_t num_base, XDataType const x_base[],
			size_t num_array, YDataType const y_base[]) :
			InterpolationWorker<XDataType, YDataType>(num_base, x_base,
					num_array, y_base), is_ascending_(
					this->x_base_[0] < this->x_base_[this->num_base_ - 1]) {
	}
	virtual ~NearestInterpolationWorker() {
	}
	virtual void Interpolate1D(size_t left_index, size_t right_index,
			size_t location, size_t offset, XDataType const x_interpolated[],
			YDataType y_interpolated[]) {
		XDataType middle_point = 0.5
				* (this->x_base_[location] + this->x_base_[location - 1]);
		YDataType y_left =
				this->y_base_[offset * this->num_base_ + location - 1];
		YDataType y_right = this->y_base_[offset * this->num_base_ + location];
		for (size_t i = left_index; i < right_index; ++i) {
			// nearest condition
			// - x_interpolated[i] == middle_point
			//     ---> nearest is y_left (left side)
			// - x_interpolated[i] < middle_point and ascending order
			//     ---> nearest is y_left (left side)
			// - x_interpolated[i] > middle_point and ascending order
			//     ---> nearest is y_right (right side)
			// - x_interpolated[i] < middle_point and descending order
			//     ---> nearest is y_right (right side)
			// - x_interpolated[i] > middle_point and descending order
			//     ---> nearest is y_left (left side)
			if ((x_interpolated[i] != middle_point)
					& ((x_interpolated[i] > middle_point) & is_ascending_)) {
				y_interpolated[i] = y_right;
			} else {
				y_interpolated[i] = y_left;
			}
		}
	}
	virtual void Interpolate1DAlongRow(size_t left_index, size_t right_index,
			size_t location, XDataType const x_interpolated[],
			YDataType y_interpolated[]) {
		XDataType middle_point = 0.5
				* (this->x_base_[location] + this->x_base_[location - 1]);
		for (size_t i = left_index; i < right_index; ++i) {
			// nearest condition
			// - x_interpolated[i] == middle_point
			//     ---> nearest is y_left (left side)
			// - x_interpolated[i] < middle_point and ascending order
			//     ---> nearest is y_left (left side)
			// - x_interpolated[i] > middle_point and ascending order
			//     ---> nearest is y_right (right side)
			// - x_interpolated[i] < middle_point and descending order
			//     ---> nearest is y_right (right side)
			// - x_interpolated[i] > middle_point and descending order
			//     ---> nearest is y_left (left side)
			size_t offset_index = 0;
			if ((x_interpolated[i] != middle_point)
					& ((x_interpolated[i] > middle_point) & is_ascending_)) {
				offset_index = 1;
			}
			YDataType *y_work = &y_interpolated[this->num_array_ * i];
			YDataType const *nearest = &this->y_base_[this->num_array_
					* (location - 1 + offset_index)];
			for (size_t j = 0; j < this->num_array_; ++j) {
				y_work[j] = nearest[j];
			}
		}
	}
private:
	bool const is_ascending_;
};

template<class XDataType, class YDataType>
class LinearInterpolationWorker: public InterpolationWorker<XDataType, YDataType> {
public:
	LinearInterpolationWorker(size_t num_base, XDataType const x_base[],
			size_t num_array, YDataType const y_base[]) :
			InterpolationWorker<XDataType, YDataType>(num_base, x_base,
					num_array, y_base) {
	}
	virtual ~LinearInterpolationWorker() {
	}
	virtual void Interpolate1D(size_t left_index, size_t right_index,
			size_t location, size_t offset, XDataType const x_interpolated[],
			YDataType y_interpolated[]) {
		size_t offset_index_left = offset * this->num_base_ + location - 1;
		XDataType dydx = static_cast<XDataType>(this->y_base_[offset_index_left
				+ 1] - this->y_base_[offset_index_left])
				/ (this->x_base_[location] - this->x_base_[location - 1]);
		YDataType base_term = this->y_base_[offset_index_left];
		for (size_t i = left_index; i < right_index; ++i) {
			y_interpolated[i] =
					base_term
							+ static_cast<YDataType>(dydx
									* (x_interpolated[i]
											- this->x_base_[location - 1]));
		}
	}
	virtual void Interpolate1DAlongRow(size_t left_index, size_t right_index,
			size_t location, XDataType const x_interpolated[],
			YDataType y_interpolated[]) {
		for (size_t i = left_index; i < right_index; ++i) {
			YDataType fraction = static_cast<YDataType>((x_interpolated[i]
					- this->x_base_[location - 1])
					/ (this->x_base_[location] - this->x_base_[location - 1]));
			YDataType *y_work = &y_interpolated[this->num_array_ * i];
			YDataType const *y_left = &this->y_base_[this->num_array_
					* (location - 1)];
			YDataType const *y_right = &this->y_base_[this->num_array_
					* location];
			for (size_t j = 0; j < this->num_array_; ++j) {
				y_work[j] = y_left[j] + fraction * (y_right[j] - y_left[j]);
			}
		}
	}
};

template<class XDataType, class YDataType>
class PolynomialInterpolationWorker: public InterpolationWorker<XDataType,
		YDataType> {
public:
	PolynomialInterpolationWorker(size_t num_base, XDataType const x_base[],
			size_t num_array, YDataType const y_base[],
			uint8_t polynomial_order) :
			InterpolationWorker<XDataType, YDataType>(num_base, x_base,
					num_array, y_base, polynomial_order), kNumElements_(
					static_cast<size_t>(this->polynomial_order_) + 1), storage_for_c_(), storage_for_d_(), c_(
					nullptr), d_(nullptr) {
	}
	virtual ~PolynomialInterpolationWorker() {
	}
	virtual void PrepareForInterpolation() {
//		storage_for_c_.reset(AllocateAndAlign<XDataType>(kNumElements_, &c_));
//		storage_for_d_.reset(AllocateAndAlign<XDataType>(kNumElements_, &d_));
		storage_for_c_.reset(
				AllocateAndAlign<XDataType>(kNumElements_ * this->num_array_,
						&c_));
		storage_for_d_.reset(
				AllocateAndAlign<XDataType>(kNumElements_ * this->num_array_,
						&d_));
	}
	virtual void Interpolate1D(size_t left_index, size_t right_index,
			size_t location, size_t offset, XDataType const x_interpolated[],
			YDataType y_interpolated[]) {
		int j = location - 1 - static_cast<int>(this->polynomial_order_) / 2;
		size_t m = this->num_base_ - 1
				- static_cast<size_t>(this->polynomial_order_);
		size_t k = static_cast<size_t>((j > 0) ? j : 0);
		k = (k > m) ? m : k;
		for (size_t i = left_index; i < right_index; ++i) {
			PerformNevilleAlgorithm(k, offset, x_interpolated[i],
					&y_interpolated[i]);
		}
	}
	virtual void Interpolate1DAlongRow(size_t left_index, size_t right_index,
			size_t location, XDataType const x_interpolated[],
			YDataType y_interpolated[]) {
		int j = location - 1 - static_cast<int>(this->polynomial_order_) / 2;
		size_t m = this->num_base_ - 1
				- static_cast<size_t>(this->polynomial_order_);
		size_t k = static_cast<size_t>((j > 0) ? j : 0);
		k = (k > m) ? m : k;
		for (size_t i = left_index; i < right_index; ++i) {
			PerformNevilleAlgorithmAlongRow(k, x_interpolated[i],
					&y_interpolated[this->num_array_ * i]);
		}
	}
private:
	void PerformNevilleAlgorithm(size_t left_index, size_t array_index,
			XDataType x_interpolated, YDataType y_interpolated[]) {

		// working pointers
		XDataType const *x_ptr = &(this->x_base_[left_index]);
		//YDataType const *y_ptr = &(this->y_base_[left_index]);
		YDataType const *y_ptr = &(this->y_base_[this->num_base_ * array_index
				+ left_index]);

		for (size_t i = 0; i < kNumElements_; ++i) {
			c_[i] = static_cast<XDataType>(y_ptr[i]);
			d_[i] = c_[i];
		}

		// Neville's algorithm
		XDataType y_interpolated_work = c_[0];
		for (size_t m = 1; m < kNumElements_; ++m) {
			// Evaluate Cm1, Cm2, Cm3, ... Cm[n-m] and Dm1, Dm2, Dm3, ... Dm[n-m].
			// Those are stored to c[0], c[1], ..., c[n-m-1] and d[0], d[1], ...,
			// d[n-m-1].
			for (size_t i = 0; i < kNumElements_ - m; ++i) {
				XDataType cd = c_[i + 1] - d_[i];
				XDataType dx = x_ptr[i] - x_ptr[i + m];
				assert(dx != 0);
				cd /= dx;
				c_[i] = (x_ptr[i] - x_interpolated) * cd;
				d_[i] = (x_ptr[i + m] - x_interpolated) * cd;
			}

			// In each step, c[0] holds Cm1 which is a correction between
			// P12...m and P12...[m+1]. Thus, the following repeated update
			// corresponds to the route P1 -> P12 -> P123 ->...-> P123...n.
			y_interpolated_work += c_[0];
		}

		*y_interpolated = static_cast<YDataType>(y_interpolated_work);
	}

	void PerformNevilleAlgorithmAlongRow(size_t left_index,
			XDataType x_interpolated, YDataType y_interpolated[]) {

		// working pointers
		XDataType const *x_ptr = &(this->x_base_[left_index]);
		//YDataType const *y_ptr = &(this->y_base_[left_index]);
//		YDataType const *y_ptr = &(this->y_base_[this->num_base_ * array_index
//				+ left_index]);

		size_t start = left_index * this->num_array_;
		size_t num_elements = kNumElements_ * this->num_array_;
		for (size_t i = 0; i < num_elements; ++i) {
			c_[i] = static_cast<XDataType>(this->y_base_[start + i]);
			d_[i] = c_[i];
		}

		// Neville's algorithm
		XDataType *y_work = nullptr;
		std::unique_ptr<void, InterpolationDeleter> storage_y(
				AllocateAndAlign(this->num_array_, &y_work));
		for (size_t i = 0; i < this->num_array_; ++i) {
			y_work[i] = c_[i];
		}
//		XDataType y_interpolated_work = c_[0];
		for (size_t m = 1; m < kNumElements_; ++m) {
			// Evaluate Cm1, Cm2, Cm3, ... Cm[n-m] and Dm1, Dm2, Dm3, ... Dm[n-m].
			// Those are stored to c[0], c[1], ..., c[n-m-1] and d[0], d[1], ...,
			// d[n-m-1].
			for (size_t i = 0; i < kNumElements_ - m; ++i) {
				XDataType dx = x_ptr[i] - x_ptr[i + m];
				assert(dx != 0);
				size_t offset = i * this->num_array_;
				for (size_t j = 0; j < this->num_array_; ++j) {
					XDataType cd = (c_[offset + this->num_array_ + j]
							- d_[offset + j]) / dx;
					c_[offset + j] = (x_ptr[i] - x_interpolated) * cd;
					d_[offset + j] = (x_ptr[i + m] - x_interpolated) * cd;
				}
			}

			// In each step, c[0] holds Cm1 which is a correction between
			// P12...m and P12...[m+1]. Thus, the following repeated update
			// corresponds to the route P1 -> P12 -> P123 ->...-> P123...n.
//			y_interpolated_work += c_[0];
			for (size_t i = 0; i < this->num_array_; ++i) {
				y_work[i] += c_[i];
			}
		}

		for (size_t i = 0; i < this->num_array_; ++i) {
			y_interpolated[i] = static_cast<YDataType>(y_work[i]);
		}
//		*y_interpolated = static_cast<YDataType>(y_interpolated_work);
	}

	size_t const kNumElements_;
	std::unique_ptr<void, InterpolationDeleter> storage_for_c_;
	std::unique_ptr<void, InterpolationDeleter> storage_for_d_;
	XDataType *c_;
	XDataType *d_;
}
;

template<class XDataType, class YDataType>
class SplineInterpolationWorker: public InterpolationWorker<XDataType, YDataType> {
public:
	SplineInterpolationWorker(size_t num_base, XDataType const x_base[],
			size_t num_array, YDataType const y_base[]) :
			InterpolationWorker<XDataType, YDataType>(num_base, x_base,
					num_array, y_base), storage_for_y_(nullptr), y_base_2nd_derivative_(
					nullptr) {
	}
	virtual ~SplineInterpolationWorker() {
	}
	virtual void PrepareForInterpolation() {
		// Derive second derivative at x_base
		bool const kIsDescending = (this->x_base_[0]
				> this->x_base_[this->num_base_ - 1]);
		storage_for_y_.reset(
				AllocateAndAlign(this->num_base_ * this->num_array_,
						&y_base_2nd_derivative_));
		for (size_t i = 0; i < this->num_array_; ++i) {
			size_t start_position = i * this->num_base_;
			YDataType const *y_base = &(this->y_base_[start_position]);
			YDataType *y_base_2nd_derivative_local =
					&(y_base_2nd_derivative_[start_position]);
			DeriveSplineCorrectionTerm<XDataType, YDataType>(kIsDescending,
					this->num_base_, this->x_base_, y_base,
					y_base_2nd_derivative_local);
		}
	}
protected:
	std::unique_ptr<void, InterpolationDeleter> storage_for_y_;
	YDataType *y_base_2nd_derivative_;
};

template<class XDataType, class YDataType>
class SplineInterpolationWorkerAscending: public SplineInterpolationWorker<
		XDataType, YDataType> {
public:
	SplineInterpolationWorkerAscending(size_t num_base,
			XDataType const x_base[], size_t num_array,
			YDataType const y_base[]) :
			SplineInterpolationWorker<XDataType, YDataType>(num_base, x_base,
					num_array, y_base) {
	}
	virtual ~SplineInterpolationWorkerAscending() {
	}
	virtual void Interpolate1D(size_t left_index, size_t right_index,
			size_t location, size_t offset, XDataType const x_interpolated[],
			YDataType y_interpolated[]) {
		assert(this->y_base_2nd_derivative_ != nullptr);
		size_t offset_index_left = offset * this->num_base_ + location - 1;
		XDataType dx = this->x_base_[location] - this->x_base_[location - 1];
		for (size_t i = left_index; i < right_index; ++i) {
			XDataType a = (this->x_base_[location] - x_interpolated[i]) / dx;
			XDataType b = (x_interpolated[i] - this->x_base_[location - 1])
					/ dx;
			y_interpolated[i] =
					static_cast<YDataType>(a * this->y_base_[offset_index_left]
							+ b * this->y_base_[offset_index_left + 1]
							+ ((a * a * a - a)
									* this->y_base_2nd_derivative_[offset_index_left]
									+ (b * b * b - b)
											* this->y_base_2nd_derivative_[offset_index_left
													+ 1]) * (dx * dx) / 6.0);

		}
	}
	virtual void Interpolate1DAlongRow(size_t left_index, size_t right_index,
			size_t location, XDataType const x_interpolated[],
			YDataType y_interpolated[]) {
		assert(this->y_base_2nd_derivative_ != nullptr);
		XDataType dx = this->x_base_[location] - this->x_base_[location - 1];
		YDataType const *y_left = &this->y_base_[(location - 1)
				* this->num_array_];
		YDataType const *y_right = &this->y_base_[location * this->num_array_];
		YDataType const *d2ydx2_left = &this->y_base_2nd_derivative_[(location
				- 1) * this->num_array_];
		YDataType const *d2ydx2_right = &this->y_base_2nd_derivative_[location
				* this->num_array_];
		for (size_t i = left_index; i < right_index; ++i) {
			XDataType a = (this->x_base_[location] - x_interpolated[i]) / dx;
			XDataType b = (x_interpolated[i] - this->x_base_[location - 1])
					/ dx;
			YDataType *y_work = &y_interpolated[i * this->num_array_];
			for (size_t j = 0; j < this->num_array_; ++j) {
				y_work[j] = static_cast<YDataType>(a * y_left[j]
						+ b * y_right[j]
						+ ((a * a * a - a) * d2ydx2_left[j]
								+ (b * b * b - b) * d2ydx2_right[j]) * (dx * dx)
								/ 6.0);
			}
		}
	}
};

template<class XDataType, class YDataType>
class SplineInterpolationWorkerAscendingAlongRow: public SplineInterpolationWorkerAscending<
		XDataType, YDataType> {
public:
	SplineInterpolationWorkerAscendingAlongRow(size_t num_base,
			XDataType const x_base[], size_t num_array,
			YDataType const y_base[]) :
			SplineInterpolationWorkerAscending<XDataType, YDataType>(num_base,
					x_base, num_array, y_base) {
	}
	virtual ~SplineInterpolationWorkerAscendingAlongRow() {
	}
	virtual void PrepareForInterpolation() {
		// Derive second derivative at x_base
//		bool const kIsDescending = (this->x_base_[0]
//				> this->x_base_[this->num_base_ - 1]);
		this->storage_for_y_.reset(
				AllocateAndAlign(this->num_base_ * this->num_array_,
						&this->y_base_2nd_derivative_));
		DeriveSplineCorrectionTermAlongRow(this->num_base_, this->x_base_,
				this->num_array_, this->y_base_, this->y_base_2nd_derivative_);
	}
};

template<class XDataType, class YDataType>
YDataType const *GetAscendingArray(size_t num_base,
		XDataType const base_array[/*num_base*/], size_t num_array,
		YDataType const unordered_array[/*num_base*num_array*/],
		std::unique_ptr<void, InterpolationDeleter> *storage,
		bool reorder_column) {
	if (base_array[0] < base_array[num_base - 1]) {
		return unordered_array;
	} else {
		YDataType *output_array = nullptr;
		storage->reset(
				AllocateAndAlign<YDataType>(num_base * num_array,
						&output_array));
		YDataType *work_array = const_cast<YDataType *>(output_array);
		if (reorder_column) {
			// reorder column
			for (size_t i = 0; i < num_array; ++i) {
				size_t start_position = num_base * i;
				size_t end_position = start_position + num_base;
				for (size_t j = start_position; j < end_position; ++j) {
					work_array[j] = unordered_array[end_position
							- (j - start_position) - 1];
				}
			}
		} else {
			// reorder row
			for (size_t i = 0; i < num_base; ++i) {
				YDataType *out_storage = &work_array[i * num_array];
				YDataType const *in_storage =
						&unordered_array[(num_base - 1 - i) * num_array];
				for (size_t j = 0; j < num_array; ++j) {
					out_storage[j] = in_storage[j];
				}
			}
		}
		return static_cast<YDataType *>(output_array);
	}
}

template<class XDataType, class YDataType>
InterpolationWorker<XDataType, YDataType> *CreateInterpolationWorker(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method, size_t num_base,
		XDataType const x_base[], size_t num_array, YDataType const y_base[],
		uint8_t polynomial_order,
		bool along_column) {
	switch (interpolation_method) {
	case LIBSAKURA_SYMBOL(InterpolationMethod_kNearest):
		return new NearestInterpolationWorker<XDataType, YDataType>(num_base,
				x_base, num_array, y_base);
	case LIBSAKURA_SYMBOL(InterpolationMethod_kLinear):
		return new LinearInterpolationWorker<XDataType, YDataType>(num_base,
				x_base, num_array, y_base);
	case LIBSAKURA_SYMBOL(InterpolationMethod_kPolynomial):
		if (polynomial_order == 0) {
			// This is special case: 0-th polynomial interpolation
			// acts like nearest interpolation
			return new NearestInterpolationWorker<XDataType, YDataType>(
					num_base, x_base, num_array, y_base);
		} else if (static_cast<size_t>(polynomial_order + 1) >= num_base) {
			// use full region for interpolation
			return new PolynomialInterpolationWorker<XDataType, YDataType>(
					num_base, x_base, num_array, y_base, num_base - 1);
		} else {
			// use sub-region around the nearest points
			return new PolynomialInterpolationWorker<XDataType, YDataType>(
					num_base, x_base, num_array, y_base, polynomial_order);
		}
	case LIBSAKURA_SYMBOL(InterpolationMethod_kSpline):
		if (along_column) {
			return new SplineInterpolationWorkerAscending<XDataType, YDataType>(
					num_base, x_base, num_array, y_base);
		} else {
			return new SplineInterpolationWorkerAscendingAlongRow<XDataType,
					YDataType>(num_base, x_base, num_array, y_base);
		}
	default:
		// invalid interpolation method type
		return nullptr;
	}
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {

//template<class XDataType, class YDataType>
//LIBSAKURA_SYMBOL(Status) ADDSUFFIX(Interpolation, ARCH_SUFFIX)<XDataType,
//		YDataType>::Interpolate1D(
//LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
//		uint8_t polynomial_order, size_t num_base,
//		XDataType const x_base[/* num_base */], size_t num_array,
//		YDataType const y_base[/* num_base * num_array*/],
//		size_t num_interpolated,
//		XDataType const x_interpolated[/* num_interpolated */],
//		YDataType y_interpolated[/*num_interpolated * num_array*/]) const {
//	assert(num_base > 0);
//	assert(num_interpolated > 0);
////	assert(LIBSAKURA_SYMBOL(IsAligned)(x_base));
////	assert(LIBSAKURA_SYMBOL(IsAligned)(y_base));
////	assert(LIBSAKURA_SYMBOL(IsAligned)(x_interpolated));
////	assert(LIBSAKURA_SYMBOL(IsAligned)(y_interpolated));
//
////	// input arrays are not aligned
////	if (!LIBSAKURA_SYMBOL(IsAligned)(x_base)
////			|| !LIBSAKURA_SYMBOL(IsAligned)(y_base)
////			|| !LIBSAKURA_SYMBOL(IsAligned)(x_interpolated)
////			|| !LIBSAKURA_SYMBOL(IsAligned)(y_interpolated)) {
////		if (Logger::IsErrorEnabled(logger)) {
////			std::ostringstream oss;
////			oss << "ERROR: input arrays are not aligned" << std::endl;
////			Logger::Error(logger, oss.str().c_str());
////		}
////		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
////	}
//
//	if (num_base == 1) {
//		// No need to interpolate, just substitute y_base[0]
//		// to all elements in y_interpolated
//		for (size_t i = 0; i < num_array; ++i) {
//			YDataType *y_interpolated_work = &y_interpolated[i
//					* num_interpolated];
//			YDataType const y_base_work = y_base[i * num_base];
//			for (size_t index = 0; index < num_interpolated; ++index) {
//				y_interpolated_work[index] = y_base_work;
//			}
//		}
//		return LIBSAKURA_SYMBOL(Status_kOK);
//	}
//
//	// Perform 1-dimensional interpolation
//
//	// make input arrays ascending order
//	std::unique_ptr<void, InterpolationDeleter> storage_for_x_base_work;
//	std::unique_ptr<void, InterpolationDeleter> storage_for_y_base_work;
//	std::unique_ptr<void, InterpolationDeleter> storage_for_x_interpolated_work;
//	XDataType const *x_base_work = GetAscendingArray<XDataType, XDataType>(
//			num_base, x_base, 1, x_base, &storage_for_x_base_work, true);
//	YDataType const *y_base_work = GetAscendingArray<XDataType, YDataType>(
//			num_base, x_base, num_array, y_base, &storage_for_y_base_work,
//			true);
//	XDataType const *x_interpolated_work = GetAscendingArray<XDataType,
//			XDataType>(num_interpolated, x_interpolated, 1, x_interpolated,
//			&storage_for_x_interpolated_work, true);
//
//	// Generate worker class
//	std::unique_ptr<InterpolationWorker<XDataType, YDataType> > interpolator(
//			CreateInterpolationWorker<XDataType, YDataType>(
//					interpolation_method, num_base, x_base_work, num_array,
//					y_base_work, polynomial_order, true));
//	if (interpolator.get() == nullptr) {
//		// failed to create interpolation worker object
//		if (Logger::IsErrorEnabled(logger)) {
//			std::ostringstream oss;
//			oss << "ERROR: Invalid interpolation method type" << std::endl;
//			Logger::Error(logger, oss.str().c_str());
//		}
//		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
//	}
//
//	// Any preparation for interpolation should be done here
//	interpolator->PrepareForInterpolation();
//
//	// Locate each element in x_base against x_interpolated
//	size_t *location_base = nullptr;
//	std::unique_ptr<void, InterpolationDeleter> storage_for_location_base(
//			AllocateAndAlign<size_t>(num_base, &location_base));
//	if (num_interpolated == 1) {
//		for (size_t i = 0; i < num_base; ++i) {
//			location_base[i] =
//					(x_base_work[i] <= x_interpolated_work[0]) ? 0 : 1;
//		}
//	} else {
//		size_t start_position = 0;
//		size_t end_position = num_interpolated - 1;
//		for (size_t i = 0; i < num_base; ++i) {
//			location_base[i] = this->Locate(start_position, end_position,
//					num_interpolated, x_interpolated_work, x_base_work[i]);
//			start_position = location_base[i];
//		}
//	}
//
//	// Outside of x_base[0]
//	for (size_t j = 0; j < num_array; ++j) {
//		YDataType *y_interpolated_work = &y_interpolated[j * num_interpolated];
//		YDataType const y_base_tmp = y_base_work[j * num_base];
//		for (size_t i = 0; i < location_base[0]; ++i) {
//			y_interpolated_work[i] = y_base_tmp;
//		}
//
//		// Between x_base[0] and x_base[num_base-1]
//		for (size_t i = 0; i < num_base - 1; ++i) {
//			if (location_base[i] != location_base[i + 1]) {
//				interpolator->Interpolate1D(location_base[i],
//						location_base[i + 1], i + 1, j, x_interpolated_work,
//						y_interpolated_work);
//			}
//		}
//
//		// Outside of x_base[num_base-1]
//		YDataType const y_base_tmp2 = y_base_work[(j + 1) * num_base - 1];
//		for (size_t i = location_base[num_base - 1]; i < num_interpolated;
//				++i) {
//			y_interpolated_work[i] = y_base_tmp2;
//		}
//	}
//
//	// swap output array
//	if (x_interpolated[0] >= x_interpolated[num_interpolated - 1]) {
//		size_t num_interpolated_half = num_interpolated / 2;
//		for (size_t j = 0; j < num_array; ++j) {
//			YDataType *y_interpolated_work = &y_interpolated[j
//					* num_interpolated];
//			for (size_t i = 0; i < num_interpolated_half; ++i) {
//				std::swap<YDataType>(y_interpolated_work[i],
//						y_interpolated_work[num_interpolated - 1 - i]);
//			}
//		}
//	}
//	return LIBSAKURA_SYMBOL(Status_kOK);
//}

/**
 * Interpolate1DAlongColumn performs 1D interpolation along column based on x_base
 * and y_base. y_base is a serial array of column-major matrix y. its memory layout
 * is assumed to be:
 *     y_base[0]     = y[0][0]
 *     y_base[1]     = y[1][0]
 *     y_base[2]     = y[2][0]
 *     ...
 *     y_base[N]     = y[N][0]
 *     y_base[N+1]   = y[1][0]
 *     ...
 *     y_base[N*M-1] = y[N][M]
 * where N and M correspond to num_interpolation_axis and num_array respectively.
 */
template<class XDataType, class YDataType>
void ADDSUFFIX(Interpolation, ARCH_SUFFIX)<XDataType, YDataType>::Interpolate1DAlongColumn(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_interpolation_axis,
		XDataType const x_base[/*num_interpolation_axis*/], size_t num_array,
		YDataType const y_base[/*num_interpolation_axis*num_array*/],
		size_t num_interpolated,
		XDataType const x_interpolated[/*num_interpolated*/],
		YDataType y_interpolated[/*num_interpolated*num_array*/]) const {
	assert(num_interpolation_axis > 0);
	assert(num_array > 0);
	assert(num_interpolated > 0);
	assert(x_base != nullptr);
	assert(y_base != nullptr);
	assert(x_interpolated != nullptr);
	assert(y_interpolated != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(x_base));
	assert(LIBSAKURA_SYMBOL(IsAligned)(y_base));
	assert(LIBSAKURA_SYMBOL(IsAligned)(x_interpolated));
	assert(LIBSAKURA_SYMBOL(IsAligned)(y_interpolated));

	if (num_interpolation_axis == 1) {
		// No need to interpolate, just substitute y_base
		// to all elements in y_interpolated
		for (size_t i = 0; i < num_array; ++i) {
			YDataType *y_interpolated_work = &y_interpolated[i
					* num_interpolated];
			YDataType const y = y_base[i * num_interpolation_axis];
			for (size_t j = 0; j < num_interpolated; ++j) {
				y_interpolated_work[j] = y;
			}
		}
		return;
	}

	// Perform 1-dimensional interpolation

	// make input arrays ascending order
	std::unique_ptr<void, InterpolationDeleter> storage_for_x_base_work;
	std::unique_ptr<void, InterpolationDeleter> storage_for_y_base_work;
	std::unique_ptr<void, InterpolationDeleter> storage_for_x_interpolated_work;
	XDataType const *x_base_work = GetAscendingArray<XDataType, XDataType>(
			num_interpolation_axis, x_base, 1, x_base, &storage_for_x_base_work,
			true);
	YDataType const *y_base_work = GetAscendingArray<XDataType, YDataType>(
			num_interpolation_axis, x_base, num_array, y_base,
			&storage_for_y_base_work, true);
	XDataType const *x_interpolated_work = GetAscendingArray<XDataType,
			XDataType>(num_interpolated, x_interpolated, 1, x_interpolated,
			&storage_for_x_interpolated_work, true);

	// Generate worker class
	std::unique_ptr<InterpolationWorker<XDataType, YDataType> > interpolator(
			CreateInterpolationWorker<XDataType, YDataType>(
					interpolation_method, num_interpolation_axis, x_base_work,
					num_array, y_base_work, polynomial_order,
					true));
	assert(interpolator.get() != nullptr);

	// Any preparation for interpolation should be done here
	interpolator->PrepareForInterpolation();

	// Locate each element in x_base against x_interpolated
	size_t *location_base = nullptr;
	std::unique_ptr<void, InterpolationDeleter> storage_for_location_base(
			AllocateAndAlign<size_t>(num_interpolation_axis, &location_base));
	if (num_interpolated == 1) {
		for (size_t i = 0; i < num_interpolation_axis; ++i) {
			location_base[i] =
					(x_base_work[i] <= x_interpolated_work[0]) ? 0 : 1;
		}
	} else {
		size_t start_position = 0;
		size_t end_position = num_interpolated - 1;
		for (size_t i = 0; i < num_interpolation_axis; ++i) {
			location_base[i] = this->Locate(start_position, end_position,
					num_interpolated, x_interpolated_work, x_base_work[i]);
			start_position = location_base[i];
		}
	}

	for (size_t j = 0; j < num_array; ++j) {
		// Outside of x_base[0]
		YDataType *y_interpolated_work = &y_interpolated[j * num_interpolated];
		YDataType const y_base_tmp = y_base_work[j * num_interpolation_axis];
		for (size_t i = 0; i < location_base[0]; ++i) {
			y_interpolated_work[i] = y_base_tmp;
		}

		// Between x_base[0] and x_base[num_interpolation_axis-1]
		for (size_t i = 0; i < num_interpolation_axis - 1; ++i) {
			if (location_base[i] != location_base[i + 1]) {
				interpolator->Interpolate1D(location_base[i],
						location_base[i + 1], i + 1, j, x_interpolated_work,
						y_interpolated_work);
			}
		}

		// Outside of x_base[num_interpolation_axis-1]
		YDataType const y_base_tmp2 = y_base_work[(j + 1)
				* num_interpolation_axis - 1];
		for (size_t i = location_base[num_interpolation_axis - 1];
				i < num_interpolated; ++i) {
			y_interpolated_work[i] = y_base_tmp2;
		}
	}

	// swap output array
	if (x_interpolated[0] >= x_interpolated[num_interpolated - 1]) {
		size_t num_interpolated_half = num_interpolated / 2;
		for (size_t j = 0; j < num_array; ++j) {
			YDataType *y_interpolated_work = &y_interpolated[j
					* num_interpolated];
			for (size_t i = 0; i < num_interpolated_half; ++i) {
				std::swap<YDataType>(y_interpolated_work[i],
						y_interpolated_work[num_interpolated - 1 - i]);
			}
		}
	}
}

/**
 * Interpolate1DAlongColumn performs 1D interpolation along row based on x_base
 * and y_base. y_base is a serial array of column-major matrix y. its memory layout
 * is assumed to be:
 *     y_base[0]     = y[0][0]
 *     y_base[1]     = y[0][1]
 *     y_base[2]     = y[0][2]
 *     ...
 *     y_base[M]     = y[0][M]
 *     y_base[M+1]   = y[1][0]
 *     ...
 *     y_base[M*N-1] = y[N][M]
 * where N and M correspond to num_interpolation_axis and num_array respectively.
 */
template<class XDataType, class YDataType>
void ADDSUFFIX(Interpolation, ARCH_SUFFIX)<XDataType, YDataType>::Interpolate1DAlongRow(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_interpolation_axis,
		XDataType const x_base[/*num_interpolation_axis*/], size_t num_array,
		YDataType const y_base[/*num_array*num_interpolation_axis*/],
		size_t num_interpolated,
		XDataType const x_interpolated[/*num_interpolated*/],
		YDataType y_interpolated[/*num_array*num_interpolated*/]) const {
	assert(num_interpolation_axis > 0);
	assert(num_array > 0);
	assert(num_interpolated > 0);
	assert(x_base != nullptr);
	assert(y_base != nullptr);
	assert(x_interpolated != nullptr);
	assert(y_interpolated != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(x_base));
	assert(LIBSAKURA_SYMBOL(IsAligned)(y_base));
	assert(LIBSAKURA_SYMBOL(IsAligned)(x_interpolated));
	assert(LIBSAKURA_SYMBOL(IsAligned)(y_interpolated));

	if (num_interpolation_axis == 1) {
		// No need to interpolate, just substitute y_base
		// to all elements in y_interpolated
		for (size_t i = 0; i < num_interpolated; ++i) {
			YDataType *y_interpolated_work = &y_interpolated[i * num_array];
			for (size_t j = 0; j < num_array; ++j) {
				y_interpolated_work[j] = y_base[j];
			}
		}
		return;
	}

	// Perform 1-dimensional interpolation

	// make input arrays ascending order
	std::unique_ptr<void, InterpolationDeleter> storage_for_x_base_work;
	std::unique_ptr<void, InterpolationDeleter> storage_for_y_base_work;
	std::unique_ptr<void, InterpolationDeleter> storage_for_x_interpolated_work;
	XDataType const *x_base_work = GetAscendingArray<XDataType, XDataType>(
			num_interpolation_axis, x_base, 1, x_base, &storage_for_x_base_work,
			false);
	YDataType const *y_base_work = GetAscendingArray<XDataType, YDataType>(
			num_interpolation_axis, x_base, num_array, y_base,
			&storage_for_y_base_work, false);
	XDataType const *x_interpolated_work = GetAscendingArray<XDataType,
			XDataType>(num_interpolated, x_interpolated, 1, x_interpolated,
			&storage_for_x_interpolated_work, false);

	// Generate worker class
	std::unique_ptr<InterpolationWorker<XDataType, YDataType> > interpolator(
			CreateInterpolationWorker<XDataType, YDataType>(
					interpolation_method, num_interpolation_axis, x_base_work,
					num_array, y_base_work, polynomial_order,
					false));
	assert(interpolator.get() != nullptr);

	// Any preparation for interpolation should be done here
	interpolator->PrepareForInterpolation();

	// Locate each element in x_base against x_interpolated
	size_t *location_base = nullptr;
	std::unique_ptr<void, InterpolationDeleter> storage_for_location_base(
			AllocateAndAlign<size_t>(num_interpolation_axis, &location_base));

	if (num_interpolated == 1) {
		for (size_t i = 0; i < num_interpolation_axis; ++i) {
			location_base[i] =
					(x_base_work[i] <= x_interpolated_work[0]) ? 0 : 1;
		}
	} else {
		size_t start_position = 0;
		size_t end_position = num_interpolated - 1;
		for (size_t i = 0; i < num_interpolation_axis; ++i) {
			location_base[i] = this->Locate(start_position, end_position,
					num_interpolated, x_interpolated_work, x_base_work[i]);
			start_position = location_base[i];
		}
	}

	// Outside of x_base[0]
	YDataType const *in = &y_base_work[0];
	for (size_t i = 0; i < location_base[0]; ++i) {
		YDataType *out = &y_interpolated[i * num_array];
		for (size_t j = 0; j < num_array; ++j) {
			out[j] = in[j];
		}
	}

	// Between x_base[0] and x_base[num_interpolation_axis-1]
	for (size_t i = 0; i < num_interpolation_axis - 1; ++i) {
		if (location_base[i] != location_base[i + 1]) {
			interpolator->Interpolate1DAlongRow(location_base[i],
					location_base[i + 1], i + 1, x_interpolated_work,
					y_interpolated);
		}
	}

	// Outside of x_base[num_interpolation_axis-1]
	in = &y_base_work[(num_interpolation_axis - 1) * num_array];
	for (size_t i = location_base[num_interpolation_axis - 1];
			i < num_interpolated; ++i) {
		YDataType *out = &y_interpolated[i * num_array];
		for (size_t j = 0; j < num_array; ++j) {
			out[j] = in[j];
		}
	}

	// swap output array
	if (x_interpolated[0] >= x_interpolated[num_interpolated - 1]) {
		size_t num_interpolated_half = num_interpolated / 2;
		for (size_t i = 0; i < num_interpolated_half; ++i) {
			YDataType *a = &y_interpolated[i * num_array];
			YDataType *b = &y_interpolated[(num_interpolated - 1 - i)
					* num_array];
			for (size_t j = 0; j < num_array; ++j) {
				std::swap<YDataType>(a[j], b[j]);
			}
		}
	}
}

template class ADDSUFFIX(Interpolation, ARCH_SUFFIX)<double, float> ;
} /* namespace LIBSAKURA_PREFIX */
