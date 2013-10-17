#include <cassert>
#include <sstream>
#include <memory>
#include <cstdalign>
#include <utility>

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

template<bool x_interpolation, class XDataType, class YDataType>
class InterpolationWorker {
public:
	InterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			polynomial_order_(polynomial_order), num_base_(num_base), num_array_(
					num_array), num_interpolated_(num_interpolated), base_position_(
					base_position), base_data_(base_data), interpolated_position_(
					interpolated_position), interpolated_data_(
					interpolated_data) {

	}
	virtual ~InterpolationWorker() {
	}
	virtual void PrepareForInterpolation() {
	}
	virtual void Interpolate1D(size_t num_location, size_t const location[]) {
	}
protected:
	uint8_t const polynomial_order_;
	size_t const num_base_;
	size_t const num_array_;
	size_t const num_interpolated_;
	XDataType const *base_position_;
	YDataType const *base_data_;
	XDataType const *interpolated_position_;
	YDataType *interpolated_data_;
};

template<bool x_interpolation, class XDataType, class YDataType>
class NearestInterpolationWorker: public InterpolationWorker<x_interpolation,
		XDataType, YDataType> {
public:
	NearestInterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			InterpolationWorker<x_interpolation, XDataType, YDataType>(
					polynomial_order, num_base, base_position, num_array,
					base_data, num_interpolated, interpolated_position,
					interpolated_data) {
	}
	virtual ~NearestInterpolationWorker() {
	}
	virtual void Interpolate1D(size_t num_location, size_t const location[]) {
		for (size_t k = 0; k < num_location - 1; ++k) {
			size_t location_left = location[k];
			size_t location_right = location[k + 1];
			size_t right_index = k + 1;
			XDataType middle_point = 0.5
					* (this->base_position_[right_index]
							+ this->base_position_[right_index - 1]);
			for (size_t j = 0; j < this->num_array_; ++j) {
				YDataType left_value = this->base_data_[j * this->num_base_
						+ right_index - 1];
				YDataType right_value = this->base_data_[j * this->num_base_
						+ right_index];
				YDataType *work = &this->interpolated_data_[j
						* this->num_interpolated_];
				for (size_t i = location_left; i < location_right; ++i) {
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
					if ((this->interpolated_position_[i] != middle_point)
							& ((this->interpolated_position_[i] > middle_point))) {
						work[i] = right_value;
					} else {
						work[i] = left_value;
					}
				}
			}
		}
	}
};

template<bool x_interpolation, class XDataType, class YDataType>
class LinearInterpolationWorker: public InterpolationWorker<x_interpolation,
		XDataType, YDataType> {
public:
	LinearInterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			InterpolationWorker<x_interpolation, XDataType, YDataType>(
					polynomial_order, num_base, base_position, num_array,
					base_data, num_interpolated, interpolated_position,
					interpolated_data) {
	}
	virtual ~LinearInterpolationWorker() {
	}
	virtual void Interpolate1D(size_t num_location, size_t const location[]) {
		for (size_t k = 0; k < num_location - 1; ++k) {
			size_t location_left = location[k];
			size_t location_right = location[k + 1];
			size_t right_index = k + 1;
			for (size_t j = 0; j < this->num_array_; ++j) {
				size_t offset_index_left = j * this->num_base_ + right_index
						- 1;
				XDataType dydx =
						static_cast<XDataType>(this->base_data_[offset_index_left
								+ 1] - this->base_data_[offset_index_left])
								/ (this->base_position_[right_index]
										- this->base_position_[right_index - 1]);
				YDataType base_term = this->base_data_[offset_index_left];
				YDataType *y_work = &this->interpolated_data_[j
						* this->num_interpolated_];
				for (size_t i = location_left; i < location_right; ++i) {
					y_work[i] = base_term
							+ static_cast<YDataType>(dydx
									* (this->interpolated_position_[i]
											- this->base_position_[right_index
													- 1]));
				}
			}
		}
	}
};

template<bool x_interpolation, class XDataType, class YDataType>
class PolynomialInterpolationWorker: public InterpolationWorker<x_interpolation,
		XDataType, YDataType> {
public:
	PolynomialInterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			InterpolationWorker<x_interpolation, XDataType, YDataType>(
					polynomial_order, num_base, base_position, num_array,
					base_data, num_interpolated, interpolated_position,
					interpolated_data), num_elements_(
					static_cast<size_t>(this->polynomial_order_) + 1), c_(
					nullptr), d_(nullptr) {
	}
	virtual ~PolynomialInterpolationWorker() {
	}
	virtual void PrepareForInterpolation() {
		storage_for_c_.reset(AllocateAndAlign<XDataType>(num_elements_, &c_));
		storage_for_d_.reset(AllocateAndAlign<XDataType>(num_elements_, &d_));
	}
	virtual void Interpolate1D(size_t num_location, size_t const location[]) {
		for (size_t k = 0; k < num_location - 1; ++k) {
			size_t location_left = location[k];
			size_t location_right = location[k + 1];
			size_t right_index = k + 1;
			for (size_t jj = 0; jj < this->num_array_; ++jj) {
				int j = right_index - 1
						- static_cast<int>(this->polynomial_order_) / 2;
				size_t m = this->num_base_ - 1
						- static_cast<size_t>(this->polynomial_order_);
				size_t kk = static_cast<size_t>((j > 0) ? j : 0);
				kk = (kk > m) ? m : kk;
				for (size_t i = location_left; i < location_right; ++i) {
					PerformNevilleAlgorithm(kk, jj,
							this->interpolated_position_[i],
							&this->interpolated_data_[jj
									* this->num_interpolated_ + i]);
				}
			}
		}
	}
private:
	void PerformNevilleAlgorithm(size_t left_index, size_t array_index,
			XDataType interpolated_position, YDataType *interpolated_data) {

		// working pointers
		XDataType const *x_ptr = &(this->base_position_[left_index]);
		//YDataType const *y_ptr = &(this->y_base_[left_index]);
		YDataType const *y_ptr = &(this->base_data_[this->num_base_
				* array_index + left_index]);

		for (size_t i = 0; i < num_elements_; ++i) {
			c_[i] = static_cast<XDataType>(y_ptr[i]);
			d_[i] = c_[i];
		}

		// Neville's algorithm
		XDataType work = c_[0];
		for (size_t m = 1; m < num_elements_; ++m) {
			// Evaluate Cm1, Cm2, Cm3, ... Cm[n-m] and Dm1, Dm2, Dm3, ... Dm[n-m].
			// Those are stored to c[0], c[1], ..., c[n-m-1] and d[0], d[1], ...,
			// d[n-m-1].
			for (size_t i = 0; i < num_elements_ - m; ++i) {
				XDataType cd = c_[i + 1] - d_[i];
				XDataType dx = x_ptr[i] - x_ptr[i + m];
				assert(dx != 0);
				cd /= dx;
				c_[i] = (x_ptr[i] - interpolated_position) * cd;
				d_[i] = (x_ptr[i + m] - interpolated_position) * cd;
			}

			// In each step, c[0] holds Cm1 which is a correction between
			// P12...m and P12...[m+1]. Thus, the following repeated update
			// corresponds to the route P1 -> P12 -> P123 ->...-> P123...n.
			work += c_[0];
		}

		*interpolated_data = static_cast<YDataType>(work);
	}

	size_t const num_elements_;
	std::unique_ptr<void, InterpolationDeleter> storage_for_c_;
	std::unique_ptr<void, InterpolationDeleter> storage_for_d_;
	XDataType *c_;
	XDataType *d_;
};

template<bool x_interpolation, class XDataType, class YDataType>
class SplineInterpolationWorker: public InterpolationWorker<x_interpolation,
		XDataType, YDataType> {
public:
	SplineInterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			InterpolationWorker<x_interpolation, XDataType, YDataType>(
					polynomial_order, num_base, base_position, num_array,
					base_data, num_interpolated, interpolated_position,
					interpolated_data), d2ydx2_(nullptr) {
	}
	virtual ~SplineInterpolationWorker() {
	}
	virtual void PrepareForInterpolation() {
		// Derive second derivative at x_base
		storage_for_d2ydx2_.reset(
				AllocateAndAlign(this->num_base_ * this->num_array_, &d2ydx2_));
		for (size_t i = 0; i < this->num_array_; ++i) {
			size_t start_position = i * this->num_base_;
			YDataType const *base_data = &(this->base_data_[start_position]);
			YDataType *d2ydx2_local = &(d2ydx2_[start_position]);
			DeriveSplineCorrectionTerm(this->num_base_, this->base_position_,
					base_data, d2ydx2_local);
		}
	}
	virtual void Interpolate1D(size_t num_location, size_t const location[]) {
		for (size_t k = 0; k < num_location - 1; ++k) {
			size_t location_left = location[k];
			size_t location_right = location[k + 1];
			size_t right_index = k + 1;
			assert(d2ydx2_ != nullptr);
			for (size_t j = 0; j < this->num_array_; ++j) {
				size_t offset_index_left = j * this->num_base_ + right_index
						- 1;
				XDataType dx = this->base_position_[right_index]
						- this->base_position_[right_index - 1];
				YDataType *work = &this->interpolated_data_[j
						* this->num_interpolated_];
				XDataType dx_factor = dx * dx / 6.0;
				for (size_t i = location_left; i < location_right; ++i) {
					XDataType a = (this->base_position_[right_index]
							- this->interpolated_position_[i]) / dx;
					XDataType b = (this->interpolated_position_[i]
							- this->base_position_[right_index - 1]) / dx;
					work[i] = static_cast<YDataType>(a
							* this->base_data_[offset_index_left]
							+ b * this->base_data_[offset_index_left + 1]
							+ ((a * a * a - a) * d2ydx2_[offset_index_left]
									+ (b * b * b - b)
											* d2ydx2_[offset_index_left + 1])
									* dx_factor);
				}
			}
		}
	}
private:
	void DeriveSplineCorrectionTerm(size_t num_base,
			XDataType const base_position[], YDataType const base_data[],
			YDataType d2ydx2[]) {
		YDataType *upper_triangular = nullptr;
		std::unique_ptr<void, InterpolationDeleter> storage_for_u(
				AllocateAndAlign<YDataType>(num_base, &upper_triangular));

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

		// Solve the system by backsubstitution and store solution to
		// d2ydx2
		for (size_t k = num_base - 2; k >= 1; --k) {
			d2ydx2[k] -= upper_triangular[k] * d2ydx2[k + 1];
		}
	}
	std::unique_ptr<void, InterpolationDeleter> storage_for_d2ydx2_;
	YDataType *d2ydx2_;
};

template<class XDataType, class YDataType>
class NearestInterpolationWorker<false, XDataType, YDataType> : public InterpolationWorker<
		false, XDataType, YDataType> {
public:
	NearestInterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			InterpolationWorker<false, XDataType, YDataType>(polynomial_order,
					num_base, base_position, num_array, base_data,
					num_interpolated, interpolated_position, interpolated_data) {
	}
	virtual ~NearestInterpolationWorker() {
	}
	virtual void Interpolate1D(size_t num_location, size_t const location[]) {
		for (size_t k = 0; k < num_location - 1; ++k) {
			size_t location_left = location[k];
			size_t location_right = location[k + 1];
			size_t right_index = k + 1;
			XDataType middle_point = 0.5
					* (this->base_position_[right_index]
							+ this->base_position_[right_index - 1]);
			for (size_t i = location_left; i < location_right; ++i) {
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
				if ((this->interpolated_position_[i] != middle_point)
						& ((this->interpolated_position_[i] > middle_point))) {
					offset_index = 1;
				}
				YDataType *work =
						&this->interpolated_data_[this->num_array_ * i];
				YDataType const *nearest = &this->base_data_[this->num_array_
						* (right_index - 1 + offset_index)];
				for (size_t j = 0; j < this->num_array_; ++j) {
					work[j] = nearest[j];
				}
			}
		}
	}
};

template<class XDataType, class YDataType>
class LinearInterpolationWorker<false, XDataType, YDataType> : public InterpolationWorker<
		false, XDataType, YDataType> {
public:
	LinearInterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			InterpolationWorker<false, XDataType, YDataType>(polynomial_order,
					num_base, base_position, num_array, base_data,
					num_interpolated, interpolated_position, interpolated_data) {
	}
	virtual ~LinearInterpolationWorker() {
	}
	virtual void Interpolate1D(size_t num_location, size_t const location[]) {
		for (size_t k = 0; k < num_location - 1; ++k) {
			size_t location_left = location[k];
			size_t location_right = location[k + 1];
			size_t right_index = k + 1;
			for (size_t i = location_left; i < location_right; ++i) {
				YDataType fraction =
						static_cast<YDataType>((this->interpolated_position_[i]
								- this->base_position_[right_index - 1])
								/ (this->base_position_[right_index]
										- this->base_position_[right_index - 1]));
				YDataType *work =
						&this->interpolated_data_[this->num_array_ * i];
				YDataType const *left_value = &this->base_data_[this->num_array_
						* (right_index - 1)];
				YDataType const *right_value =
						&this->base_data_[this->num_array_ * right_index];
				for (size_t j = 0; j < this->num_array_; ++j) {
					work[j] = left_value[j]
							+ fraction * (right_value[j] - left_value[j]);
				}
			}
		}
	}
};

template<class XDataType, class YDataType>
class PolynomialInterpolationWorker<false, XDataType, YDataType> : public InterpolationWorker<
		false, XDataType, YDataType> {
public:
	PolynomialInterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			InterpolationWorker<false, XDataType, YDataType>(polynomial_order,
					num_base, base_position, num_array, base_data,
					num_interpolated, interpolated_position, interpolated_data), num_elements_(
					static_cast<size_t>(this->polynomial_order_) + 1), c_(
					nullptr), d_(nullptr), work_(nullptr) {
	}
	virtual ~PolynomialInterpolationWorker() {
	}
	virtual void PrepareForInterpolation() {
		storage_for_c_.reset(
				AllocateAndAlign<XDataType>(num_elements_ * this->num_array_,
						&c_));
		storage_for_d_.reset(
				AllocateAndAlign<XDataType>(num_elements_ * this->num_array_,
						&d_));
		storage_for_work_.reset(
				AllocateAndAlign<XDataType>(this->num_array_, &work_));
	}
	virtual void Interpolate1D(size_t num_location, size_t const location[]) {
		for (size_t k = 0; k < num_location - 1; ++k) {
			size_t location_left = location[k];
			size_t location_right = location[k + 1];
			size_t right_index = k + 1;
			int j = right_index - 1
					- static_cast<int>(this->polynomial_order_) / 2;
			size_t m = this->num_base_ - 1
					- static_cast<size_t>(this->polynomial_order_);
			size_t kk = static_cast<size_t>((j > 0) ? j : 0);
			kk = (kk > m) ? m : kk;
			for (size_t i = location_left; i < location_right; ++i) {
				PerformNevilleAlgorithm(kk, this->interpolated_position_[i],
						&this->interpolated_data_[this->num_array_ * i]);
			}
		}
	}
private:
	void PerformNevilleAlgorithm(size_t left_index,
			XDataType interpolated_position, YDataType *interpolated_data) {
		// working pointers
		XDataType const *x_ptr = &(this->base_position_[left_index]);

		size_t start = left_index * this->num_array_;
		size_t num_elements = num_elements_ * this->num_array_;
		for (size_t i = 0; i < num_elements; ++i) {
			c_[i] = static_cast<XDataType>(this->base_data_[start + i]);
			d_[i] = c_[i];
		}

		// Neville's algorithm
		for (size_t i = 0; i < this->num_array_; ++i) {
			work_[i] = c_[i];
		}
		for (size_t m = 1; m < num_elements_; ++m) {
			// Evaluate Cm1, Cm2, Cm3, ... Cm[n-m] and Dm1, Dm2, Dm3, ... Dm[n-m].
			// Those are stored to c[0], c[1], ..., c[n-m-1] and d[0], d[1], ...,
			// d[n-m-1].
			for (size_t i = 0; i < num_elements_ - m; ++i) {
				XDataType dx = x_ptr[i] - x_ptr[i + m];
				assert(dx != 0);
				size_t offset = i * this->num_array_;
				for (size_t j = 0; j < this->num_array_; ++j) {
					XDataType cd = (c_[offset + this->num_array_ + j]
							- d_[offset + j]) / dx;
					c_[offset + j] = (x_ptr[i] - interpolated_position) * cd;
					d_[offset + j] = (x_ptr[i + m] - interpolated_position)
							* cd;
				}
			}

			// In each step, c[0] holds Cm1 which is a correction between
			// P12...m and P12...[m+1]. Thus, the following repeated update
			// corresponds to the route P1 -> P12 -> P123 ->...-> P123...n.
			for (size_t i = 0; i < this->num_array_; ++i) {
				work_[i] += c_[i];
			}
		}

		for (size_t i = 0; i < this->num_array_; ++i) {
			interpolated_data[i] = static_cast<YDataType>(work_[i]);
		}
	}

	size_t const num_elements_;
	std::unique_ptr<void, InterpolationDeleter> storage_for_c_;
	std::unique_ptr<void, InterpolationDeleter> storage_for_d_;
	std::unique_ptr<void, InterpolationDeleter> storage_for_work_;
	XDataType *c_;
	XDataType *d_;
	XDataType *work_;
};

template<class XDataType, class YDataType>
class SplineInterpolationWorker<false, XDataType, YDataType> : public InterpolationWorker<
		false, XDataType, YDataType> {
public:
	SplineInterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			InterpolationWorker<false, XDataType, YDataType>(polynomial_order,
					num_base, base_position, num_array, base_data,
					num_interpolated, interpolated_position, interpolated_data), d2ydx2_(
					nullptr) {
	}
	virtual ~SplineInterpolationWorker() {
	}
	virtual void PrepareForInterpolation() {
		// Derive second derivative at x_base
		storage_for_d2ydx2_.reset(
				AllocateAndAlign(this->num_base_ * this->num_array_, &d2ydx2_));
		DeriveSplineCorrectionTerm(this->num_base_, this->base_position_,
				this->num_array_, this->base_data_, d2ydx2_);
	}
	virtual void Interpolate1D(size_t num_location, size_t const location[]) {
		for (size_t k = 0; k < num_location - 1; ++k) {
			size_t location_left = location[k];
			size_t location_right = location[k + 1];
			size_t right_index = k + 1;
			assert(d2ydx2_ != nullptr);
			XDataType dx = this->base_position_[right_index]
					- this->base_position_[right_index - 1];
			YDataType const *left_value = &this->base_data_[(right_index - 1)
					* this->num_array_];
			YDataType const *right_value = &this->base_data_[right_index
					* this->num_array_];
			YDataType const *d2ydx2_left = &d2ydx2_[(right_index - 1)
					* this->num_array_];
			YDataType const *d2ydx2_right = &d2ydx2_[right_index
					* this->num_array_];
			for (size_t i = location_left; i < location_right; ++i) {
				XDataType a = (this->base_position_[right_index]
						- this->interpolated_position_[i]) / dx;
				XDataType b = (this->interpolated_position_[i]
						- this->base_position_[right_index - 1]) / dx;
				YDataType *work =
						&this->interpolated_data_[i * this->num_array_];
				XDataType aaa = (a * a * a - a) * dx * dx / 6.0;
				XDataType bbb = (b * b * b - b) * dx * dx / 6.0;
				for (size_t j = 0; j < this->num_array_; ++j) {
					work[j] = static_cast<YDataType>(a * left_value[j]
							+ b * right_value[j]
							+ (aaa * d2ydx2_left[j] + bbb * d2ydx2_right[j]));
				}
			}
		}
	}
private:
	void DeriveSplineCorrectionTerm(size_t num_base,
			XDataType const base_position[], size_t num_array,
			YDataType const base_data[], YDataType d2ydx2[]) {
		YDataType *upper_triangular = nullptr;
		std::unique_ptr<void, InterpolationDeleter> storage_for_u(
				AllocateAndAlign<YDataType>(num_base * num_array,
						&upper_triangular));

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

		// Solve the system by backsubstitution and store solution to
		// d2ydx2
		for (size_t k = num_base - 2; k >= 1; --k) {
			for (size_t j = 0; j < num_array; ++j) {
				d2ydx2[k * num_array + j] -= upper_triangular[k * num_array + j]
						* d2ydx2[(k + 1) * num_array + j];
			}
		}
	}
	std::unique_ptr<void, InterpolationDeleter> storage_for_d2ydx2_;
	YDataType *d2ydx2_;
};

template<bool x_interpolation, class XDataType, class YDataType>
InterpolationWorker<x_interpolation, XDataType, YDataType> *CreateInterpolationWorker(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_base,
		XDataType const base_position[], size_t num_array,
		YDataType const base_data[], size_t num_interpolated,
		XDataType const interpolated_position[],
		YDataType interpolated_data[]) {
	switch (interpolation_method) {
	case LIBSAKURA_SYMBOL(InterpolationMethod_kNearest):
		return new NearestInterpolationWorker<x_interpolation, XDataType,
				YDataType>(polynomial_order, num_base, base_position, num_array,
				base_data, num_interpolated, interpolated_position,
				interpolated_data);
	case LIBSAKURA_SYMBOL(InterpolationMethod_kLinear):
		return new LinearInterpolationWorker<x_interpolation, XDataType,
				YDataType>(polynomial_order, num_base, base_position, num_array,
				base_data, num_interpolated, interpolated_position,
				interpolated_data);
	case LIBSAKURA_SYMBOL(InterpolationMethod_kPolynomial):
		if (polynomial_order == 0) {
			// This is special case: 0-th polynomial interpolation
			// acts like nearest interpolation
			return new NearestInterpolationWorker<x_interpolation, XDataType,
					YDataType>(polynomial_order, num_base, base_position,
					num_array, base_data, num_interpolated,
					interpolated_position, interpolated_data);
		} else if (static_cast<size_t>(polynomial_order + 1) >= num_base) {
			// use full region for interpolation
			return new PolynomialInterpolationWorker<x_interpolation, XDataType,
					YDataType>(num_base - 1, num_base, base_position, num_array,
					base_data, num_interpolated, interpolated_position,
					interpolated_data);
		} else {
			// use sub-region around the nearest points
			return new PolynomialInterpolationWorker<x_interpolation, XDataType,
					YDataType>(polynomial_order, num_base, base_position,
					num_array, base_data, num_interpolated,
					interpolated_position, interpolated_data);
		}
	case LIBSAKURA_SYMBOL(InterpolationMethod_kSpline):
		return new SplineInterpolationWorker<x_interpolation, XDataType,
				YDataType>(polynomial_order, num_base, base_position, num_array,
				base_data, num_interpolated, interpolated_position,
				interpolated_data);
	default:
		// invalid interpolation method type
		return nullptr;
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

// default, effectively for x_interpolation is true
template<bool x_interpolation, class YDataType>
class Substitutor {
public:
	static inline void SubstituteSingleBaseData(size_t num_base,
			size_t num_array, size_t num_interpolated,
			YDataType const base_data[], YDataType interpolated_data[]) {
		for (size_t i = 0; i < num_array; ++i) {
			YDataType *work = &interpolated_data[i * num_interpolated];
			YDataType const value = base_data[i * num_base];
			for (size_t j = 0; j < num_interpolated; ++j) {
				work[j] = value;
			}
		}
	}
	static inline void SubstituteLeftMostData(size_t location, size_t num_base,
			size_t num_array, size_t num_interpolated,
			YDataType const base_data[], YDataType interpolated_data[]) {
		for (size_t j = 0; j < num_array; ++j) {
			YDataType *work = &interpolated_data[j * num_interpolated];
			YDataType const value = base_data[j * num_base];
			for (size_t i = 0; i < location; ++i) {
				work[i] = value;
			}
		}
	}
	static inline void SubstituteRightMostData(size_t location, size_t num_base,
			size_t num_array, size_t num_interpolated,
			YDataType const base_data[], YDataType interpolated_data[]) {
		for (size_t j = 0; j < num_array; ++j) {
			YDataType *work = &interpolated_data[j * num_interpolated];
			YDataType const value = base_data[(j + 1) * num_base - 1];
			for (size_t i = location; i < num_interpolated; ++i) {
				work[i] = value;
			}
		}
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
};

// for x_interpolatin is false
template<class YDataType>
class Substitutor<false, YDataType> {
public:
	static inline void SubstituteSingleBaseData(size_t num_base,
			size_t num_array, size_t num_interpolated,
			YDataType const base_data[], YDataType interpolated_data[]) {
		for (size_t i = 0; i < num_interpolated; ++i) {
			YDataType *work = &interpolated_data[i * num_array];
			for (size_t j = 0; j < num_array; ++j) {
				work[j] = base_data[j];
			}
		}
	}
	static inline void SubstituteLeftMostData(size_t location, size_t num_base,
			size_t num_array, size_t num_interpolated,
			YDataType const base_data[], YDataType interpolated_data[]) {
		for (size_t i = 0; i < location; ++i) {
			YDataType *work = &interpolated_data[i * num_array];
			for (size_t j = 0; j < num_array; ++j) {
				work[j] = base_data[j];
			}
		}
	}
	static inline void SubstituteRightMostData(size_t location, size_t num_base,
			size_t num_array, size_t num_interpolated,
			YDataType const base_data_[], YDataType interpolated_data_[]) {
		YDataType const *value = &base_data_[(num_base - 1) * num_array];
		for (size_t i = location; i < num_interpolated; ++i) {
			YDataType *work = &interpolated_data_[i * num_array];
			for (size_t j = 0; j < num_array; ++j) {
				work[j] = value[j];
			}
		}
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
};

template<bool x_interpolation, class XDataType, class YDataType>
void Interpolate1D(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_base,
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
		Substitutor<x_interpolation, YDataType>::SubstituteSingleBaseData(
				num_base, num_array, num_interpolated, base_data,
				interpolated_data);
		return;
	}

	std::unique_ptr<void, InterpolationDeleter> storage_for_base_position;
	std::unique_ptr<void, InterpolationDeleter> storage_for_base_data;
	std::unique_ptr<void, InterpolationDeleter> storage_for_interpolated_position;
	XDataType const *base_position_work =
			GetAscendingArray<XDataType, XDataType>(num_base, base_position, 1,
					base_position, &storage_for_base_position, x_interpolation);
	YDataType const *base_data_work = GetAscendingArray<XDataType, YDataType>(
			num_base, base_position, num_array, base_data,
			&storage_for_base_data, x_interpolation);
	XDataType const *interpolated_position_work = GetAscendingArray<XDataType,
			XDataType>(num_interpolated, interpolated_position, 1,
			interpolated_position, &storage_for_interpolated_position,
			x_interpolation);

	// Generate worker class
	std::unique_ptr<InterpolationWorker<x_interpolation, XDataType, YDataType> > interpolator(
			CreateInterpolationWorker<x_interpolation, XDataType, YDataType>(
					interpolation_method, polynomial_order, num_base,
					base_position_work, num_array, base_data_work,
					num_interpolated, interpolated_position_work,
					interpolated_data));
	assert(interpolator.get() != nullptr);

	// Perform 1-dimensional interpolation
	// Any preparation for interpolation should be done here
	interpolator->PrepareForInterpolation();

	// Locate each element in x_base against x_interpolated
	size_t *location_base = nullptr;
	std::unique_ptr<void, InterpolationDeleter> storage_for_location_base(
			AllocateAndAlign<size_t>(num_base, &location_base));
	size_t num_location_base = Locate<XDataType>(num_interpolated, num_base,
			interpolated_position_work, base_position_work, location_base);

	// Outside of x_base[0]
	Substitutor<x_interpolation, YDataType>::SubstituteLeftMostData(
			location_base[0], num_base, num_array, num_interpolated,
			base_data_work, interpolated_data);

	// Between x_base[0] and x_base[num_x_base-1]
	interpolator->Interpolate1D(num_location_base, location_base);

	// Outside of x_base[num_x_base-1]
	Substitutor<x_interpolation, YDataType>::SubstituteRightMostData(
			location_base[num_location_base - 1], num_base, num_array,
			num_interpolated, base_data_work, interpolated_data);

	// swap output array
	if (interpolated_position[0]
			> interpolated_position[num_interpolated - 1]) {
		Substitutor<x_interpolation, YDataType>::SwapResult(num_array,
				num_interpolated, interpolated_data);
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
void ADDSUFFIX(Interpolation, ARCH_SUFFIX)<XDataType, YDataType>::Interpolate1DX(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_x_base,
		XDataType const x_base[/*num_x_base*/], size_t num_y,
		YDataType const data_base[/*num_x_base*num_y*/],
		size_t num_x_interpolated,
		XDataType const x_interpolated[/*num_x_interpolated*/],
		YDataType data_interpolated[/*num_x_interpolated*num_y*/]) const {
	Interpolate1D<true, XDataType, YDataType>(interpolation_method,
			polynomial_order, num_x_base, x_base, num_y, data_base,
			num_x_interpolated, x_interpolated, data_interpolated);
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
void ADDSUFFIX(Interpolation, ARCH_SUFFIX)<XDataType, YDataType>::Interpolate1DY(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_y_base,
		XDataType const y_base[/*num_y_base*/], size_t num_x,
		YDataType const data_base[/*num_y_base*num_x*/],
		size_t num_y_interpolated,
		XDataType const y_interpolated[/*num_y_interpolated*/],
		YDataType data_interpolated[/*num_y_interpolated*num_x*/]) const {
	Interpolate1D<false, XDataType, YDataType>(interpolation_method,
			polynomial_order, num_y_base, y_base, num_x, data_base,
			num_y_interpolated, y_interpolated, data_interpolated);
}

template class ADDSUFFIX(Interpolation, ARCH_SUFFIX)<double, float> ;
} /* namespace LIBSAKURA_PREFIX */
