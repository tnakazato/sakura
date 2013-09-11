#include <cassert>
#include <iostream>
#include <memory>
#include <cstdalign>

#include <libsakura/sakura.h>
#include <libsakura/optimized_implementation_factory_impl.h>
#include <libsakura/localdef.h>

namespace {

using namespace LIBSAKURA_PREFIX;

// Polynomial interpolation using Neville's algorithm
template<class XDataType, class YDataType>
void PerformNevilleAlgorithm(XDataType const x_base[], YDataType const y_base[],
		size_t left_index, size_t num_elements, XDataType x_interpolated,
		YDataType *y_interpolated) {

	// working pointers
	XDataType const *x_ptr = &x_base[left_index];
	YDataType const *y_ptr = &y_base[left_index];

	// storage for C and D in Neville's algorithm
	size_t sakura_alignment = LIBSAKURA_SYMBOL(GetAlignment)();
	size_t elements_in_arena = num_elements + sakura_alignment - 1;
	std::unique_ptr<YDataType[]> storage_for_c(new YDataType[elements_in_arena]);
	std::unique_ptr<YDataType[]> storage_for_d(new YDataType[elements_in_arena]);
	YDataType *c = reinterpret_cast<YDataType*>(LIBSAKURA_SYMBOL(AlignAny)(
			sizeof(YDataType) * elements_in_arena, storage_for_c.get(),
			sizeof(YDataType) * num_elements));
	YDataType *d = reinterpret_cast<YDataType*>(LIBSAKURA_SYMBOL(AlignAny)(
			sizeof(YDataType) * elements_in_arena, storage_for_d.get(),
			sizeof(YDataType) * num_elements));

	for (size_t i = 0; i < num_elements; ++i) {
		c[i] = y_ptr[i];
		d[i] = y_ptr[i];
	}

	// Neville's algorithm
	YDataType y_interpolated_work = c[0];
	for (size_t m = 1; m < num_elements; ++m) {
		// Evaluate Cm1, Cm2, Cm3, ... Cm[n-m] and Dm1, Dm2, Dm3, ... Dm[n-m].
		// Those are stored to c[0], c[1], ..., c[n-m-1] and d[0], d[1], ...,
		// d[n-m-1].
		for (size_t i = 0; i < num_elements - m; ++i) {
			XDataType cd = static_cast<XDataType>(c[i + 1] - d[i]);
			XDataType dx = x_ptr[i] - x_ptr[i + m];
			assert(dx != 0);
			cd /= dx;
			c[i] = (x_ptr[i] - x_interpolated) * cd;
			d[i] = (x_ptr[i + m] - x_interpolated) * cd;
		}

		// In each step, c[0] holds Cm1 which is a correction between
		// P12...m and P12...[m+1]. Thus, the following repeated update
		// corresponds to the route P1 -> P12 -> P123 ->...-> P123...n.
		y_interpolated_work += c[0];
	}

	*y_interpolated = y_interpolated_work;
}

template<class XDataType, class YDataType>
void DeriveSplineCorrectionTerm(bool is_descending, size_t num_base,
		XDataType const x_base[], YDataType const y_base[],
		YDataType y_base_2nd_derivative[]) {
	size_t sakura_alignment = LIBSAKURA_SYMBOL(GetAlignment)();
	size_t elements_in_arena = num_base + sakura_alignment - 2;
	std::unique_ptr<YDataType[]> storage_for_u(
			new YDataType[elements_in_arena]);
	YDataType *upper_triangular =
			const_cast<YDataType *>(LIBSAKURA_SYMBOL(AlignFloat)(
					elements_in_arena, storage_for_u.get(), num_base - 1));

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
YDataType DoSplineInterpolation(size_t lower_index, size_t upper_index,
		size_t lower_index_correct, size_t upper_index_correct,
		XDataType x_interpolated, XDataType const x_base[],
		YDataType const y_base[], float const y_base_2nd_derivative[]) {
	XDataType dx = x_base[upper_index] - x_base[lower_index];
	XDataType a = (x_base[upper_index] - x_interpolated) / dx;
	XDataType b = (x_interpolated - x_base[lower_index]) / dx;
	return static_cast<YDataType>(a * y_base[lower_index]
			+ b * y_base[upper_index]
			+ ((a * a * a - a) * y_base_2nd_derivative[lower_index_correct]
					+ (b * b * b - b)
							* y_base_2nd_derivative[upper_index_correct])
					* (dx * dx) / 6.0);
}

template<class XDataType, class YDataType>
class InterpolationWorker {
public:
	InterpolationWorker(size_t num_base, XDataType const x_base[],
			YDataType const y_base[], int polynomial_order = -1) :
			num_base_(num_base), x_base_(x_base), y_base_(y_base), polynomial_order_(
					polynomial_order) {
	}
	virtual ~InterpolationWorker() {
	}
	virtual void PrepareForInterpolation() {
	}
	virtual void Interpolate(size_t location, XDataType x_interpolated,
			YDataType *y_interpolated) = 0;
protected:
	size_t num_base_;
	XDataType const *x_base_;
	YDataType const *y_base_;
	int const polynomial_order_;
};

template<class XDataType, class YDataType>
class NearestInterpolationWorker: public InterpolationWorker<XDataType,
		YDataType> {
public:
	NearestInterpolationWorker(size_t num_base, XDataType const x_base[],
			YDataType const y_base[]) :
			InterpolationWorker<XDataType, YDataType>(num_base, x_base, y_base) {
	}
	virtual ~NearestInterpolationWorker() {
	}
	virtual void Interpolate(size_t location, XDataType x_interpolated,
			YDataType *y_interpolated) {
		XDataType dx = static_cast<XDataType>(fabs(
				this->x_base_[location] - this->x_base_[location - 1]));
		XDataType dx_right = static_cast<XDataType>(fabs(
				x_interpolated - this->x_base_[location]));
		// nearest condition
		// if dx_right / dx > 0.5, index will be (location - 1) which means nearest
		// is left side. Otherwise, index will be (location), right side.
		*y_interpolated = this->y_base_[location
				- static_cast<int>((dx_right) / dx + 0.5)];
	}
};

template<class XDataType, class YDataType>
class LinearInterpolationWorker: public InterpolationWorker<XDataType, YDataType> {
public:
	LinearInterpolationWorker(size_t num_base, XDataType const x_base[],
			YDataType const y_base[]) :
			InterpolationWorker<XDataType, YDataType>(num_base, x_base, y_base) {
	}
	virtual ~LinearInterpolationWorker() {
	}
	virtual void Interpolate(size_t location, XDataType x_interpolated,
			YDataType *y_interpolated) {
		*y_interpolated =
				this->y_base_[location - 1]
						+ (this->y_base_[location] - this->y_base_[location - 1])
								* (x_interpolated - this->x_base_[location - 1])
								/ (this->x_base_[location]
										- this->x_base_[location - 1]);
	}
};

template<class XDataType, class YDataType>
class PolynomialInterpolationWorker: public InterpolationWorker<XDataType,
		YDataType> {
public:
	PolynomialInterpolationWorker(size_t num_base, XDataType const x_base[],
			YDataType const y_base[], int polynomial_order) :
			InterpolationWorker<XDataType, YDataType>(num_base, x_base, y_base,
					polynomial_order) {
	}
	virtual ~PolynomialInterpolationWorker() {
	}
	virtual void Interpolate(size_t location, XDataType x_interpolated,
			YDataType *y_interpolated) {
		int j = location - 1 - this->polynomial_order_ / 2;
		size_t m = this->num_base_ - 1
				- static_cast<size_t>(this->polynomial_order_);
		size_t k = static_cast<size_t>((j > 0) ? j : 0);
		k = (k > m) ? m : k;
		size_t num_elements = static_cast<size_t>(this->polynomial_order_ + 1);
		// call polynomial interpolation
		PerformNevilleAlgorithm<XDataType, YDataType>(this->x_base_, this->y_base_, k,
				num_elements, x_interpolated, y_interpolated);
	}
};

template<class XDataType, class YDataType>
class SplineInterpolationWorker: public InterpolationWorker<XDataType, YDataType> {
public:
	SplineInterpolationWorker(size_t num_base, XDataType const x_base[],
			YDataType const y_base[]) :
			InterpolationWorker<XDataType, YDataType>(num_base, x_base, y_base), storage_for_y_(
					nullptr), y_base_2nd_derivative_(nullptr) {
	}
	virtual ~SplineInterpolationWorker() {
	}
	virtual void PrepareForInterpolation() {
		// Derive second derivative at x_base
		const bool kIsDescending = (this->x_base_[0]
				> this->x_base_[this->num_base_ - 1]);
		size_t sakura_alignment = LIBSAKURA_SYMBOL(GetAlignment)();
		size_t elements_in_arena = this->num_base_ + sakura_alignment - 1;
		storage_for_y_.reset(new YDataType[elements_in_arena]);
		y_base_2nd_derivative_ =
				const_cast<YDataType *>(reinterpret_cast<YDataType const *>(LIBSAKURA_SYMBOL(AlignAny)(
						sizeof(YDataType) * elements_in_arena,
						storage_for_y_.get(),
						sizeof(YDataType) * this->num_base_)));
		DeriveSplineCorrectionTerm<XDataType, YDataType>(kIsDescending,
				this->num_base_, this->x_base_, this->y_base_,
				y_base_2nd_derivative_);
	}
protected:
	std::unique_ptr<YDataType[]> storage_for_y_;
	YDataType *y_base_2nd_derivative_;
};

template<class XDataType, class YDataType>
class SplineInterpolationWorkerDescending: public SplineInterpolationWorker<
		XDataType, YDataType> {
public:
	SplineInterpolationWorkerDescending(size_t num_base,
			XDataType const x_base[], YDataType const y_base[]) :
			SplineInterpolationWorker<XDataType, YDataType>(num_base, x_base,
					y_base) {
	}
	virtual ~SplineInterpolationWorkerDescending() {
	}
	virtual void Interpolate(size_t location, XDataType x_interpolated,
			YDataType *y_interpolated) {
		assert(this->y_base_2nd_derivative_ != nullptr);
		*y_interpolated = DoSplineInterpolation<XDataType, YDataType>(location,
				location - 1, this->num_base_ - location - 1,
				this->num_base_ - location - 2, x_interpolated, this->x_base_,
				this->y_base_, this->y_base_2nd_derivative_);
	}
};

template<class XDataType, class YDataType>
class SplineInterpolationWorkerAscending: public SplineInterpolationWorker<
		XDataType, YDataType> {
public:
	SplineInterpolationWorkerAscending(size_t num_base,
			XDataType const x_base[], YDataType const y_base[]) :
			SplineInterpolationWorker<XDataType, YDataType>(num_base, x_base,
					y_base) {
	}
	virtual ~SplineInterpolationWorkerAscending() {
	}
	virtual void Interpolate(size_t location, XDataType x_interpolated,
			YDataType *y_interpolated) {
		assert(this->y_base_2nd_derivative_ != nullptr);
		*y_interpolated = DoSplineInterpolation<XDataType, YDataType>(
				location - 1, location, location - 1, location, x_interpolated,
				this->x_base_, this->y_base_, this->y_base_2nd_derivative_);
	}
};

template<class XDataType, class YDataType>
InterpolationWorker<XDataType, YDataType> *CreateInterpolationWorker(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method, size_t num_base,
		XDataType const x_base[], YDataType const y_base[],
		int polynomial_order) {
	switch (interpolation_method) {
	case LIBSAKURA_SYMBOL(InterpolationMethod_kNearest):
		return new NearestInterpolationWorker<XDataType, YDataType>(num_base,
				x_base, y_base);
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kLinear):
		return new LinearInterpolationWorker<XDataType, YDataType>(num_base,
				x_base, y_base);
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kPolynomial):
		if (polynomial_order == 0) {
			// This is special case: 0-th polynomial interpolation
			// acts like nearest interpolation
			return new NearestInterpolationWorker<XDataType, YDataType>(
					num_base, x_base, y_base);
		} else if (static_cast<size_t>(polynomial_order + 1) >= num_base) {
			// use full region for interpolation
			return new PolynomialInterpolationWorker<XDataType, YDataType>(
					num_base, x_base, y_base, num_base - 1);
		} else {
			// use sub-region around the nearest points
			return new PolynomialInterpolationWorker<XDataType, YDataType>(
					num_base, x_base, y_base, polynomial_order);
		}
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kSpline):
		return new SplineInterpolationWorkerAscending<XDataType, YDataType>(
				num_base, x_base, y_base);
		break;
	default:
		// invalid interpolation method type
		return nullptr;
		break;
	}
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {

template<class XDataType, class YDataType>
LIBSAKURA_SYMBOL(Status) ADDSUFFIX(Interpolation, ARCH_SUFFIX)<XDataType, YDataType>::Interpolate1d(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		int polynomial_order, size_t num_base,
		XDataType const x_base[/* num_base */],
		YDataType const y_base[/* num_base */], size_t num_interpolated,
		XDataType const x_interpolated[/* num_interpolated */],
		YDataType y_interpolated[/*num_interpolated*/]) const {
	assert(num_base > 0);
	assert(num_interpolated > 0);
	assert(LIBSAKURA_SYMBOL(IsAligned)(x_base));
	assert(LIBSAKURA_SYMBOL(IsAligned)(y_base));
	assert(LIBSAKURA_SYMBOL(IsAligned)(x_interpolated));
	assert(LIBSAKURA_SYMBOL(IsAligned)(y_interpolated));

	if (num_base == 1) {
		// No need to interpolate, just substitute y_base[0]
		// to all elements in y_interpolated
		for (size_t index = 0; index < num_interpolated; ++index) {
			y_interpolated[index] = y_base[0];
		}
		return LIBSAKURA_SYMBOL(Status_kOK);
	} else {
		size_t sakura_alignment = LIBSAKURA_SYMBOL(GetAlignment)();

		// make input arrays ascending order
		size_t elements_in_arena = num_base + sakura_alignment - 1;
		std::unique_ptr<XDataType[]> storage_for_x_base_work;
		std::unique_ptr<YDataType[]> storage_for_y_base_work;
		const XDataType *x_base_work;
		const YDataType *y_base_work;
		if (x_base[0] < x_base[num_base - 1]) {
			x_base_work = const_cast<XDataType *>(x_base);
			y_base_work = const_cast<YDataType *>(y_base);
		} else {
			storage_for_x_base_work.reset(new XDataType[elements_in_arena]);
			x_base_work = reinterpret_cast<XDataType *>(LIBSAKURA_SYMBOL(AlignAny)(
					sizeof(XDataType) * elements_in_arena,
					storage_for_x_base_work.get(), sizeof(XDataType) * num_base));
			storage_for_y_base_work.reset(new YDataType[elements_in_arena]);
			y_base_work =
					reinterpret_cast<YDataType *>(LIBSAKURA_SYMBOL(AlignAny)(
							sizeof(YDataType) * elements_in_arena,
							storage_for_y_base_work.get(),
							sizeof(YDataType) * num_base));
			XDataType *x_base_work2 = const_cast<XDataType *>(x_base_work);
			YDataType *y_base_work2 = const_cast<YDataType *>(y_base_work);
			for (size_t i = 0; i < num_base; ++i) {
				x_base_work2[i] = x_base[num_base - 1 - i];
				y_base_work2[i] = y_base[num_base - 1 - i];
			}
		}
		elements_in_arena = num_interpolated + sakura_alignment - 1;
		std::unique_ptr<XDataType[]> storage_for_x_interpolated_work;
		std::unique_ptr<YDataType[]> storage_for_y_interpolated_work;
		const XDataType *x_interpolated_work;
		YDataType *y_interpolated_work;
		bool interpolated_is_ascending = (x_interpolated[0]
				< x_interpolated[num_interpolated - 1]);
		if (interpolated_is_ascending) {
			x_interpolated_work = const_cast<XDataType *>(x_interpolated);
			y_interpolated_work = y_interpolated;
		} else {
			storage_for_x_interpolated_work.reset(
					new XDataType[elements_in_arena]);
			x_interpolated_work =
					reinterpret_cast<XDataType *>(LIBSAKURA_SYMBOL(AlignAny)(
							sizeof(XDataType) * elements_in_arena,
							storage_for_x_interpolated_work.get(),
							sizeof(XDataType) * num_interpolated));
			storage_for_y_interpolated_work.reset(
					new YDataType[elements_in_arena]);
			y_interpolated_work =
					reinterpret_cast<YDataType *>(LIBSAKURA_SYMBOL(AlignAny)(
							sizeof(YDataType) * elements_in_arena,
							storage_for_y_interpolated_work.get(),
							sizeof(YDataType) * num_interpolated));
			XDataType *x_interpolated_work2 =
					const_cast<XDataType *>(x_interpolated_work);
			for (size_t i = 0; i < num_interpolated; ++i) {
				x_interpolated_work2[i] = x_interpolated[num_interpolated - 1
						- i];
			}
		}

		// Generate worker class
		std::unique_ptr<InterpolationWorker<XDataType, YDataType> > interpolator(
				CreateInterpolationWorker<XDataType, YDataType>(
						interpolation_method, num_base, x_base_work,
						y_base_work, polynomial_order));
		if (interpolator.get() == nullptr) {
			// failed to create interpolation worker object
			// probably due to the invalid method type
			std::cerr << "ERROR: Invalid interpolation method type"
					<< std::endl;
			return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		}

		// Any preparation for interpolation should be done here
		interpolator->PrepareForInterpolation();

		elements_in_arena = num_base + sakura_alignment - 1;
		std::unique_ptr<size_t[]> storage_for_location_base(
				new size_t[elements_in_arena]);
		size_t *location_base =
				reinterpret_cast<size_t *>(LIBSAKURA_SYMBOL(AlignAny)(
						sizeof(size_t) * elements_in_arena,
						storage_for_location_base.get(),
						sizeof(size_t) * num_base));

		// Locate each element in x_base against x_interpolated
		size_t start_position = 0;
		size_t end_position = num_interpolated - 1;
		for (size_t i = 0; i < num_base; ++i) {
			location_base[i] = this->Locate(start_position, end_position,
					num_interpolated, x_interpolated_work, x_base_work[i]);
			start_position = location_base[i];
		}

		// Outside of x_base[0]
		size_t left_index = 0;
		size_t right_index = location_base[0];
		for (size_t i = left_index; i < right_index; ++i) {
			y_interpolated_work[i] = y_base_work[0];
		}

		// Between x_base[0] and x_base[num_base-1]
		for (size_t i = 0; i < num_base - 1; ++i) {
			left_index = location_base[i];
			right_index = location_base[i + 1];
			for (size_t j = left_index; j < right_index; ++j) {
				interpolator->Interpolate(i + 1, x_interpolated_work[j],
						&y_interpolated_work[j]);
			}
		}

		// Outside of x_base[num_base-1]
		left_index = location_base[num_base - 1];
		right_index = num_interpolated;
		for (size_t i = left_index; i < right_index; ++i) {
			y_interpolated_work[i] = y_base_work[num_base - 1];
		}

		// copy work array to output array
		if (!interpolated_is_ascending) {
			for (size_t i = 0; i < num_interpolated; ++i) {
				y_interpolated[i] =
						y_interpolated_work[num_interpolated - 1 - i];
			}
		}
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

template class ADDSUFFIX(Interpolation, ARCH_SUFFIX)<double, float> ;
} /* namespace LIBSAKURA_PREFIX */
