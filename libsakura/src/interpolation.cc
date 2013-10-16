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

template<class XDataType, class YDataType>
class InterpolationWorkerTemplate {
public:
	InterpolationWorkerTemplate(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/],
			bool x_interpolation) :
			polynomial_order_(polynomial_order), num_base_(num_base), num_array_(
					num_array), num_interpolated_(num_interpolated), base_position_(
					base_position), base_data_(base_data), interpolated_position_(
					interpolated_position), base_position_work_(nullptr), base_data_work_(
					nullptr), interpolated_position_work_(nullptr), interpolated_data_(
					interpolated_data), x_interpolation_(x_interpolation), is_ascending_(
			false) {
	}
	virtual ~InterpolationWorkerTemplate() {
	}
	void Execute() {
		assert(num_base_ > 0);
		assert(num_array_ > 0);
		assert(num_interpolated_ > 0);
		assert(base_position_ != nullptr);
		assert(base_data_ != nullptr);
		assert(interpolated_position_ != nullptr);
		assert(interpolated_data_ != nullptr);
		assert(LIBSAKURA_SYMBOL(IsAligned)(base_position_));
		assert(LIBSAKURA_SYMBOL(IsAligned)(base_data_));
		assert(LIBSAKURA_SYMBOL(IsAligned)(interpolated_position_));
		assert(LIBSAKURA_SYMBOL(IsAligned)(interpolated_data_));

		if (num_base_ == 1) {
			// No need to interpolate, just substitute base_data
			// to all elements in y_interpolated
			SubstituteSingleBaseData();
			return;
		}

		// Perform 1-dimensional interpolation
		is_ascending_ = interpolated_position_[0]
				<= interpolated_position_[num_interpolated_ - 1];

		// Any preparation for interpolation should be done here
		PrepareForInterpolation();

		// Locate each element in x_base against x_interpolated
		size_t *location_base = nullptr;
		std::unique_ptr<void, InterpolationDeleter> storage_for_location_base(
				AllocateAndAlign<size_t>(num_base_, &location_base));
		size_t num_location_base = Locate(location_base);

		// Outside of x_base[0]
		SubstituteLeftMostData(location_base[0]);

		// Between x_base[0] and x_base[num_x_base-1]
		ExecuteInterpolation(location_base, num_location_base);

		// Outside of x_base[num_x_base-1]
		SubstituteRightMostData(location_base[num_location_base - 1]);

		// swap output array
		if (!is_ascending_) {
			SwapResult();
		}
	}
protected:
	virtual void PrepareForInterpolation() {
		// make input arrays ascending order
		base_position_work_ = GetAscendingArray<XDataType, XDataType>(num_base_,
				base_position_, 1, base_position_, &storage_for_base_position_,
				x_interpolation_);
		base_data_work_ = GetAscendingArray<XDataType, YDataType>(num_base_,
				base_position_, num_array_, base_data_, &storage_for_base_data_,
				x_interpolation_);
		interpolated_position_work_ = GetAscendingArray<XDataType, XDataType>(
				num_interpolated_, interpolated_position_, 1,
				interpolated_position_, &storage_for_interpolated_position_,
				x_interpolation_);
	}
	virtual void SubstituteSingleBaseData() = 0;
	virtual void SubstituteLeftMostData(size_t location) = 0;
	virtual void SubstituteRightMostData(size_t location) = 0;
	virtual void ExecuteInterpolation(size_t const location[],
			size_t num_location) {
		for (size_t i = 0; i < num_location - 1; ++i) {
			Interpolate1D(location[i], location[i + 1], i + 1);
		}
	}
	virtual void Interpolate1D(size_t location_left, size_t location_right,
			size_t right_index) = 0;
	virtual void SwapResult() = 0;

	size_t Locate(size_t location_base[]) {
		size_t num_location_base = 0;
		if (num_interpolated_ == 1) {
			if (base_position_work_[num_base_ - 1]
					<= interpolated_position_work_[0]) {
				num_location_base = 1;
				location_base[0] = 0;
			} else if (base_position_work_[0]
					>= interpolated_position_work_[0]) {
				num_location_base = 1;
				location_base[0] = 1;
			} else {
				num_location_base = 2;
				location_base[0] = 0;
				location_base[1] = 1;
			}
		} else {
			size_t start_position = 0;
			size_t end_position = num_interpolated_ - 1;
			size_t previous_location = num_interpolated_ + 1;
			for (size_t i = 0; i < num_base_; ++i) {
				size_t location = LocateData(start_position, end_position,
						num_interpolated_, interpolated_position_work_,
						base_position_work_[i]);
				if (location != previous_location) {
					location_base[num_location_base] = location;
					num_location_base += 1;
					start_position = location;
					previous_location = location;
				}
			}
		}
		return num_location_base;
	}

protected:
	size_t LocateData(size_t start_position, size_t end_position,
			size_t num_base, XDataType const base_position[],
			XDataType located_position) {
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
			return LocateData(end_position, num_base - 1, num_base,
					base_position, located_position);
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

	uint8_t const polynomial_order_;
	size_t const num_base_;
	size_t const num_array_;
	size_t const num_interpolated_;
	XDataType const *base_position_;
	YDataType const *base_data_;
	XDataType const *interpolated_position_;
	XDataType const *base_position_work_;
	YDataType const *base_data_work_;
	XDataType const *interpolated_position_work_;
	YDataType *interpolated_data_;bool x_interpolation_;bool is_ascending_;
	std::unique_ptr<void, InterpolationDeleter> storage_for_base_position_;
	std::unique_ptr<void, InterpolationDeleter> storage_for_base_data_;
	std::unique_ptr<void, InterpolationDeleter> storage_for_interpolated_position_;
};

template<class XDataType, class YDataType>
class XInterpolationWorker: public InterpolationWorkerTemplate<XDataType,
		YDataType> {
public:
	XInterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			InterpolationWorkerTemplate<XDataType, YDataType>(polynomial_order,
					num_base, base_position, num_array, base_data,
					num_interpolated, interpolated_position, interpolated_data,
					true) {
	}
	virtual ~XInterpolationWorker() {
	}
protected:
	virtual void SubstituteSingleBaseData() {
		for (size_t i = 0; i < this->num_array_; ++i) {
			YDataType *work = &this->interpolated_data_[i
					* this->num_interpolated_];
			YDataType const value = this->base_data_[i * this->num_base_];
			for (size_t j = 0; j < this->num_interpolated_; ++j) {
				work[j] = value;
			}
		}
	}
	virtual void SubstituteLeftMostData(size_t location) {
		for (size_t j = 0; j < this->num_array_; ++j) {
			YDataType *work = &this->interpolated_data_[j
					* this->num_interpolated_];
			YDataType const value = this->base_data_work_[j * this->num_base_];
			for (size_t i = 0; i < location; ++i) {
				work[i] = value;
			}
		}
	}
	virtual void SubstituteRightMostData(size_t location) {
		for (size_t j = 0; j < this->num_array_; ++j) {
			YDataType *work = &this->interpolated_data_[j
					* this->num_interpolated_];
			YDataType const value = this->base_data_work_[(j + 1)
					* this->num_base_ - 1];
			for (size_t i = location; i < this->num_interpolated_; ++i) {
				work[i] = value;
			}
		}
	}
	virtual void SwapResult() {
		size_t middle_point = this->num_interpolated_ / 2;
		size_t right_edge = this->num_interpolated_ - 1;
		for (size_t j = 0; j < this->num_array_; ++j) {
			YDataType *work = &this->interpolated_data_[j
					* this->num_interpolated_];
			for (size_t i = 0; i < middle_point; ++i) {
				std::swap<YDataType>(work[i], work[right_edge - i]);
			}
		}
	}
};

template<class XDataType, class YDataType>
class XNearestInterpolationWorker: public XInterpolationWorker<XDataType,
		YDataType> {
public:
	XNearestInterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			XInterpolationWorker<XDataType, YDataType>(polynomial_order,
					num_base, base_position, num_array, base_data,
					num_interpolated, interpolated_position, interpolated_data) {
	}
	virtual ~XNearestInterpolationWorker() {
	}
protected:
	virtual void Interpolate1D(size_t location_left, size_t location_right,
			size_t right_index) {
		XDataType middle_point = 0.5
				* (this->base_position_work_[right_index]
						+ this->base_position_work_[right_index - 1]);
		for (size_t j = 0; j < this->num_array_; ++j) {
			YDataType left_value = this->base_data_work_[j * this->num_base_
					+ right_index - 1];
			YDataType right_value = this->base_data_work_[j * this->num_base_
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
				if ((this->interpolated_position_work_[i] != middle_point)
						& ((this->interpolated_position_work_[i] > middle_point))) {
					work[i] = right_value;
				} else {
					work[i] = left_value;
				}
			}
		}
	}
};

template<class XDataType, class YDataType>
class XLinearInterpolationWorker: public XInterpolationWorker<XDataType,
		YDataType> {
public:
	XLinearInterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			XInterpolationWorker<XDataType, YDataType>(polynomial_order,
					num_base, base_position, num_array, base_data,
					num_interpolated, interpolated_position, interpolated_data) {
	}
	virtual ~XLinearInterpolationWorker() {
	}
protected:
	virtual void Interpolate1D(size_t location_left, size_t location_right,
			size_t right_index) {
		for (size_t j = 0; j < this->num_array_; ++j) {
			size_t offset_index_left = j * this->num_base_ + right_index - 1;
			XDataType dydx =
					static_cast<XDataType>(this->base_data_work_[offset_index_left
							+ 1] - this->base_data_work_[offset_index_left])
							/ (this->base_position_work_[right_index]
									- this->base_position_work_[right_index - 1]);
			YDataType base_term = this->base_data_work_[offset_index_left];
			YDataType *y_work = &this->interpolated_data_[j
					* this->num_interpolated_];
			for (size_t i = location_left; i < location_right; ++i) {
				y_work[i] = base_term
						+ static_cast<YDataType>(dydx
								* (this->interpolated_position_work_[i]
										- this->base_position_work_[right_index
												- 1]));
			}
		}
	}
};

template<class XDataType, class YDataType>
class XPolynomialInterpolationWorker: public XInterpolationWorker<XDataType,
		YDataType> {
public:
	XPolynomialInterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			XInterpolationWorker<XDataType, YDataType>(polynomial_order,
					num_base, base_position, num_array, base_data,
					num_interpolated, interpolated_position, interpolated_data), num_elements_(
					static_cast<size_t>(this->polynomial_order_) + 1), c_(
					nullptr), d_(nullptr) {
	}
	virtual ~XPolynomialInterpolationWorker() {
	}
	virtual void PrepareForInterpolation() {
		InterpolationWorkerTemplate<XDataType, YDataType>::PrepareForInterpolation();
		storage_for_c_.reset(AllocateAndAlign<XDataType>(num_elements_, &c_));
		storage_for_d_.reset(AllocateAndAlign<XDataType>(num_elements_, &d_));
	}
protected:
	virtual void Interpolate1D(size_t location_left, size_t location_right,
			size_t right_index) {
		for (size_t jj = 0; jj < this->num_array_; ++jj) {
			int j = right_index - 1
					- static_cast<int>(this->polynomial_order_) / 2;
			size_t m = this->num_base_ - 1
					- static_cast<size_t>(this->polynomial_order_);
			size_t k = static_cast<size_t>((j > 0) ? j : 0);
			k = (k > m) ? m : k;
			for (size_t i = location_left; i < location_right; ++i) {
				PerformNevilleAlgorithm(k, jj,
						this->interpolated_position_work_[i],
						&this->interpolated_data_[jj * this->num_interpolated_
								+ i]);
			}
		}
	}
private:
	void PerformNevilleAlgorithm(size_t left_index, size_t array_index,
			XDataType interpolated_position, YDataType *interpolated_data) {

		// working pointers
		XDataType const *x_ptr = &(this->base_position_work_[left_index]);
		//YDataType const *y_ptr = &(this->y_base_[left_index]);
		YDataType const *y_ptr = &(this->base_data_work_[this->num_base_
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

template<class XDataType, class YDataType>
class XSplineInterpolationWorker: public XInterpolationWorker<XDataType,
		YDataType> {
public:
	XSplineInterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			XInterpolationWorker<XDataType, YDataType>(polynomial_order,
					num_base, base_position, num_array, base_data,
					num_interpolated, interpolated_position, interpolated_data), d2ydx2_(
					nullptr) {
	}
	virtual ~XSplineInterpolationWorker() {
	}
	virtual void PrepareForInterpolation() {
		InterpolationWorkerTemplate<XDataType, YDataType>::PrepareForInterpolation();
		// Derive second derivative at x_base
		storage_for_d2ydx2_.reset(
				AllocateAndAlign(this->num_base_ * this->num_array_, &d2ydx2_));
		for (size_t i = 0; i < this->num_array_; ++i) {
			size_t start_position = i * this->num_base_;
			YDataType const *base_data =
					&(this->base_data_work_[start_position]);
			YDataType *d2ydx2_local = &(d2ydx2_[start_position]);
			DeriveSplineCorrectionTerm(this->num_base_,
					this->base_position_work_, base_data, d2ydx2_local);
		}
	}
protected:
	virtual void Interpolate1D(size_t location_left, size_t location_right,
			size_t right_index) {
		assert(d2ydx2_ != nullptr);
		for (size_t j = 0; j < this->num_array_; ++j) {
			size_t offset_index_left = j * this->num_base_ + right_index - 1;
			XDataType dx = this->base_position_work_[right_index]
					- this->base_position_work_[right_index - 1];
			YDataType *work = &this->interpolated_data_[j
					* this->num_interpolated_];
			for (size_t i = location_left; i < location_right; ++i) {
				XDataType a = (this->base_position_work_[right_index]
						- this->interpolated_position_work_[i]) / dx;
				XDataType b = (this->interpolated_position_work_[i]
						- this->base_position_work_[right_index - 1]) / dx;
				work[i] = static_cast<YDataType>(a
						* this->base_data_work_[offset_index_left]
						+ b * this->base_data_work_[offset_index_left + 1]
						+ ((a * a * a - a) * d2ydx2_[offset_index_left]
								+ (b * b * b - b)
										* d2ydx2_[offset_index_left + 1])
								* (dx * dx) / 6.0);

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
}
;

template<class XDataType, class YDataType>
class YInterpolationWorker: public InterpolationWorkerTemplate<XDataType,
		YDataType> {
public:
	YInterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			InterpolationWorkerTemplate<XDataType, YDataType>(polynomial_order,
					num_base, base_position, num_array, base_data,
					num_interpolated, interpolated_position, interpolated_data,
					false) {
	}
	virtual ~YInterpolationWorker() {
	}
protected:
	virtual void SubstituteSingleBaseData() {
		for (size_t i = 0; i < this->num_interpolated_; ++i) {
			YDataType *work = &this->interpolated_data_[i * this->num_array_];
			for (size_t j = 0; j < this->num_array_; ++j) {
				work[j] = this->base_data_[j];
			}
		}
	}
	virtual void SubstituteLeftMostData(size_t location) {
		for (size_t i = 0; i < location; ++i) {
			YDataType *work = &this->interpolated_data_[i * this->num_array_];
			for (size_t j = 0; j < this->num_array_; ++j) {
				work[j] = this->base_data_work_[j];
			}
		}
	}
	virtual void SubstituteRightMostData(size_t location) {
		YDataType const *value = &this->base_data_work_[(this->num_base_ - 1)
				* this->num_array_];
		for (size_t i = location; i < this->num_interpolated_; ++i) {
			YDataType *work = &this->interpolated_data_[i * this->num_array_];
			for (size_t j = 0; j < this->num_array_; ++j) {
				work[j] = value[j];
			}
		}
	}
	virtual void SwapResult() {
		size_t middle_point = this->num_interpolated_ / 2;
		for (size_t i = 0; i < middle_point; ++i) {
			YDataType *a = &this->interpolated_data_[i * this->num_array_];
			YDataType *b = &this->interpolated_data_[(this->num_interpolated_
					- 1 - i) * this->num_array_];
			for (size_t j = 0; j < this->num_array_; ++j) {
				std::swap<YDataType>(a[j], b[j]);
			}
		}
	}
};

template<class XDataType, class YDataType>
class YNearestInterpolationWorker: public YInterpolationWorker<XDataType,
		YDataType> {
public:
	YNearestInterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			YInterpolationWorker<XDataType, YDataType>(polynomial_order,
					num_base, base_position, num_array, base_data,
					num_interpolated, interpolated_position, interpolated_data) {
	}
	virtual ~YNearestInterpolationWorker() {
	}
protected:
	virtual void Interpolate1D(size_t location_left, size_t location_right,
			size_t right_index) {
		XDataType middle_point = 0.5
				* (this->base_position_work_[right_index]
						+ this->base_position_work_[right_index - 1]);
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
			if ((this->interpolated_position_work_[i] != middle_point)
					& ((this->interpolated_position_work_[i] > middle_point))) {
				offset_index = 1;
			}
			YDataType *work = &this->interpolated_data_[this->num_array_ * i];
			YDataType const *nearest = &this->base_data_work_[this->num_array_
					* (right_index - 1 + offset_index)];
			for (size_t j = 0; j < this->num_array_; ++j) {
				work[j] = nearest[j];
			}
		}
	}
};

template<class XDataType, class YDataType>
class YLinearInterpolationWorker: public YInterpolationWorker<XDataType,
		YDataType> {
public:
	YLinearInterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			YInterpolationWorker<XDataType, YDataType>(polynomial_order,
					num_base, base_position, num_array, base_data,
					num_interpolated, interpolated_position, interpolated_data) {
	}
	virtual ~YLinearInterpolationWorker() {
	}
protected:
	virtual void Interpolate1D(size_t location_left, size_t location_right,
			size_t right_index) {
		for (size_t i = location_left; i < location_right; ++i) {
			YDataType fraction =
					static_cast<YDataType>((this->interpolated_position_work_[i]
							- this->base_position_work_[right_index - 1])
							/ (this->base_position_work_[right_index]
									- this->base_position_work_[right_index - 1]));
			YDataType *work = &this->interpolated_data_[this->num_array_ * i];
			YDataType const *left_value =
					&this->base_data_work_[this->num_array_ * (right_index - 1)];
			YDataType const *right_value =
					&this->base_data_work_[this->num_array_ * right_index];
			for (size_t j = 0; j < this->num_array_; ++j) {
				work[j] = left_value[j]
						+ fraction * (right_value[j] - left_value[j]);
			}
		}
	}
};

template<class XDataType, class YDataType>
class YPolynomialInterpolationWorker: public YInterpolationWorker<XDataType,
		YDataType> {
public:
	YPolynomialInterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			YInterpolationWorker<XDataType, YDataType>(polynomial_order,
					num_base, base_position, num_array, base_data,
					num_interpolated, interpolated_position, interpolated_data), num_elements_(
					static_cast<size_t>(this->polynomial_order_) + 1), c_(
					nullptr), d_(nullptr) {
	}
	virtual ~YPolynomialInterpolationWorker() {
	}
	virtual void PrepareForInterpolation() {
		InterpolationWorkerTemplate<XDataType, YDataType>::PrepareForInterpolation();
		storage_for_c_.reset(
				AllocateAndAlign<XDataType>(num_elements_ * this->num_array_,
						&c_));
		storage_for_d_.reset(
				AllocateAndAlign<XDataType>(num_elements_ * this->num_array_,
						&d_));
	}
protected:
	virtual void Interpolate1D(size_t location_left, size_t location_right,
			size_t right_index) {
		int j = right_index - 1 - static_cast<int>(this->polynomial_order_) / 2;
		size_t m = this->num_base_ - 1
				- static_cast<size_t>(this->polynomial_order_);
		size_t k = static_cast<size_t>((j > 0) ? j : 0);
		k = (k > m) ? m : k;
		for (size_t i = location_left; i < location_right; ++i) {
			PerformNevilleAlgorithm(k, this->interpolated_position_work_[i],
					&this->interpolated_data_[this->num_array_ * i]);
		}
	}
private:
	void PerformNevilleAlgorithm(size_t left_index,
			XDataType interpolated_position, YDataType *interpolated_data) {
		// working pointers
		XDataType const *x_ptr = &(this->base_position_work_[left_index]);

		size_t start = left_index * this->num_array_;
		size_t num_elements = num_elements_ * this->num_array_;
		for (size_t i = 0; i < num_elements; ++i) {
			c_[i] = static_cast<XDataType>(this->base_data_work_[start + i]);
			d_[i] = c_[i];
		}

		// Neville's algorithm
		XDataType *work = nullptr;
		std::unique_ptr<void, InterpolationDeleter> storage_y(
				AllocateAndAlign(this->num_array_, &work));
		for (size_t i = 0; i < this->num_array_; ++i) {
			work[i] = c_[i];
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
				work[i] += c_[i];
			}
		}

		for (size_t i = 0; i < this->num_array_; ++i) {
			interpolated_data[i] = static_cast<YDataType>(work[i]);
		}
	}

	size_t const num_elements_;
	std::unique_ptr<void, InterpolationDeleter> storage_for_c_;
	std::unique_ptr<void, InterpolationDeleter> storage_for_d_;
	XDataType *c_;
	XDataType *d_;
};

template<class XDataType, class YDataType>
class YSplineInterpolationWorker: public YInterpolationWorker<XDataType,
		YDataType> {
public:
	YSplineInterpolationWorker(uint8_t polynomial_order, size_t num_base,
			XDataType const base_position[/*num_base*/], size_t num_array,
			YDataType const base_data[/*num_base*num_array*/],
			size_t num_interpolated,
			XDataType const interpolated_position[/*num_interpolated*/],
			YDataType interpolated_data[/*num_interpolated*num_array*/]) :
			YInterpolationWorker<XDataType, YDataType>(polynomial_order,
					num_base, base_position, num_array, base_data,
					num_interpolated, interpolated_position, interpolated_data), d2ydx2_(
					nullptr) {
	}
	virtual ~YSplineInterpolationWorker() {
	}
	virtual void PrepareForInterpolation() {
		InterpolationWorkerTemplate<XDataType, YDataType>::PrepareForInterpolation();
		// Derive second derivative at x_base
		storage_for_d2ydx2_.reset(
				AllocateAndAlign(this->num_base_ * this->num_array_, &d2ydx2_));
		DeriveSplineCorrectionTerm(this->num_base_, this->base_position_work_,
				this->num_array_, this->base_data_work_, d2ydx2_);
	}
protected:
	virtual void Interpolate1D(size_t location_left, size_t location_right,
			size_t right_index) {
		assert(d2ydx2_ != nullptr);
		XDataType dx = this->base_position_work_[right_index]
				- this->base_position_work_[right_index - 1];
		YDataType const *left_value = &this->base_data_work_[(right_index - 1)
				* this->num_array_];
		YDataType const *right_value = &this->base_data_work_[right_index
				* this->num_array_];
		YDataType const *d2ydx2_left = &d2ydx2_[(right_index - 1)
				* this->num_array_];
		YDataType const *d2ydx2_right = &d2ydx2_[right_index * this->num_array_];
		for (size_t i = location_left; i < location_right; ++i) {
			XDataType a = (this->base_position_work_[right_index]
					- this->interpolated_position_work_[i]) / dx;
			XDataType b = (this->interpolated_position_work_[i]
					- this->base_position_work_[right_index - 1]) / dx;
			YDataType *work = &this->interpolated_data_[i * this->num_array_];
			for (size_t j = 0; j < this->num_array_; ++j) {
				work[j] = static_cast<YDataType>(a * left_value[j]
						+ b * right_value[j]
						+ ((a * a * a - a) * d2ydx2_left[j]
								+ (b * b * b - b) * d2ydx2_right[j]) * (dx * dx)
								/ 6.0);
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

template<class XDataType, class YDataType>
InterpolationWorkerTemplate<XDataType, YDataType> *CreateInterpolationWorker(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_base,
		XDataType const base_position[], size_t num_array,
		YDataType const base_data[], size_t num_interpolated,
		XDataType const interpolated_position[], YDataType interpolated_data[],
		bool x_interpolation) {
	if (x_interpolation) {
		switch (interpolation_method) {
		case LIBSAKURA_SYMBOL(InterpolationMethod_kNearest):
			return new XNearestInterpolationWorker<XDataType, YDataType>(
					polynomial_order, num_base, base_position, num_array,
					base_data, num_interpolated, interpolated_position,
					interpolated_data);
		case LIBSAKURA_SYMBOL(InterpolationMethod_kLinear):
			return new XLinearInterpolationWorker<XDataType, YDataType>(
					polynomial_order, num_base, base_position, num_array,
					base_data, num_interpolated, interpolated_position,
					interpolated_data);
		case LIBSAKURA_SYMBOL(InterpolationMethod_kPolynomial):
			if (polynomial_order == 0) {
				// This is special case: 0-th polynomial interpolation
				// acts like nearest interpolation
				return new XNearestInterpolationWorker<XDataType, YDataType>(
						polynomial_order, num_base, base_position, num_array,
						base_data, num_interpolated, interpolated_position,
						interpolated_data);
			} else if (static_cast<size_t>(polynomial_order + 1) >= num_base) {
				// use full region for interpolation
				return new XPolynomialInterpolationWorker<XDataType, YDataType>(
						num_base - 1, num_base, base_position, num_array,
						base_data, num_interpolated, interpolated_position,
						interpolated_data);
			} else {
				// use sub-region around the nearest points
				return new XPolynomialInterpolationWorker<XDataType, YDataType>(
						polynomial_order, num_base, base_position, num_array,
						base_data, num_interpolated, interpolated_position,
						interpolated_data);
			}
		case LIBSAKURA_SYMBOL(InterpolationMethod_kSpline):
			return new XSplineInterpolationWorker<XDataType, YDataType>(
					polynomial_order, num_base, base_position, num_array,
					base_data, num_interpolated, interpolated_position,
					interpolated_data);
		default:
			// invalid interpolation method type
			return nullptr;
		}
	} else {
		switch (interpolation_method) {
		case LIBSAKURA_SYMBOL(InterpolationMethod_kNearest):
			return new YNearestInterpolationWorker<XDataType, YDataType>(
					polynomial_order, num_base, base_position, num_array,
					base_data, num_interpolated, interpolated_position,
					interpolated_data);
		case LIBSAKURA_SYMBOL(InterpolationMethod_kLinear):
			return new YLinearInterpolationWorker<XDataType, YDataType>(
					polynomial_order, num_base, base_position, num_array,
					base_data, num_interpolated, interpolated_position,
					interpolated_data);
		case LIBSAKURA_SYMBOL(InterpolationMethod_kPolynomial):
			if (polynomial_order == 0) {
				// This is special case: 0-th polynomial interpolation
				// acts like nearest interpolation
				return new YNearestInterpolationWorker<XDataType, YDataType>(
						polynomial_order, num_base, base_position, num_array,
						base_data, num_interpolated, interpolated_position,
						interpolated_data);
			} else if (static_cast<size_t>(polynomial_order + 1) >= num_base) {
				// use full region for interpolation
				return new YPolynomialInterpolationWorker<XDataType, YDataType>(
						num_base - 1, num_base, base_position, num_array,
						base_data, num_interpolated, interpolated_position,
						interpolated_data);
			} else {
				// use sub-region around the nearest points
				return new YPolynomialInterpolationWorker<XDataType, YDataType>(
						polynomial_order, num_base, base_position, num_array,
						base_data, num_interpolated, interpolated_position,
						interpolated_data);
			}
		case LIBSAKURA_SYMBOL(InterpolationMethod_kSpline):
			return new YSplineInterpolationWorker<XDataType, YDataType>(
					polynomial_order, num_base, base_position, num_array,
					base_data, num_interpolated, interpolated_position,
					interpolated_data);
		default:
			// invalid interpolation method type
			return nullptr;
		}
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
	// Generate worker class
	std::unique_ptr<InterpolationWorkerTemplate<XDataType, YDataType> > interpolator(
			CreateInterpolationWorker<XDataType, YDataType>(
					interpolation_method, polynomial_order, num_x_base, x_base,
					num_y, data_base, num_x_interpolated, x_interpolated,
					data_interpolated,
					true));
	assert(interpolator.get() != nullptr);
	interpolator->Execute();
}

/**
 * Interpolate1DAlongColumn performs 1D interpolation along row based on y_base
 * and data_base. data_base is a serial array of column-major matrix data.
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
		XDataType const x_interpolated[/*num_y_interpolated*/],
		YDataType data_interpolated[/*num_y_interpolated*num_x*/]) const {
	// Generate worker class
	std::unique_ptr<InterpolationWorkerTemplate<XDataType, YDataType> > interpolator(
			CreateInterpolationWorker<XDataType, YDataType>(
					interpolation_method, polynomial_order, num_y_base, y_base,
					num_x, data_base, num_y_interpolated, x_interpolated,
					data_interpolated,
					false));
	assert(interpolator.get() != nullptr);
	interpolator->Execute();
}

template class ADDSUFFIX(Interpolation, ARCH_SUFFIX)<double, float> ;
} /* namespace LIBSAKURA_PREFIX */
