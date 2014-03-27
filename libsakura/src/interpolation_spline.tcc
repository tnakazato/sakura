#ifndef _LIBSAKURA_INTERPOLATION_SPLINE_TCC_
#define _LIBSAKURA_INTERPOLATION_SPLINE_TCC_

#include "interpolation_utils.tcc"

namespace {

// TODO: documentation must be added
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
		AllocateAndAlign < YDataType > (num_base * num_array, &holder_[0]);
		AllocateAndAlign < YDataType > (num_base, &holder_[1]);
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
			size_t const location[], size_t offset) {
		YDataType *d2ydx2 = holder_[0].pointer;
		assert(d2ydx2 != nullptr);
		for (size_t k = 1; k < num_location; ++k) {
			size_t left_index = k + offset - 1;
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
class SplineYInterpolatorImpl {
public:
	void PrepareForInterpolation(uint8_t polynomial_order, size_t num_base,
			size_t num_array, XDataType const base_position[],
			YDataType const base_data[]) {
		// Derive second derivative at x_base
		AllocateAndAlign < YDataType > (num_base * num_array, &holder_);
		DeriveSplineCorrectionTerm(num_base, base_position, num_array,
				base_data, holder_.pointer);
	}
	void Interpolate1D(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t num_location,
			size_t const location[], size_t offset) {
		YDataType *d2ydx2 = holder_.pointer;
		assert(d2ydx2 != nullptr);
		for (size_t k = 1; k < num_location; ++k) {
			size_t left_index = k + offset - 1;
			XDataType dx = base_position[left_index + 1]
					- base_position[left_index];
			YDataType const *left_value = &base_data[left_index * num_array];
			YDataType const *right_value = &base_data[(left_index + 1)
					* num_array];
			YDataType const *d2ydx2_left = &d2ydx2[left_index * num_array];
			YDataType const *d2ydx2_right =
					&d2ydx2[(left_index + 1) * num_array];
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
	}
private:
	void DeriveSplineCorrectionTerm(size_t num_base,
			XDataType const base_position[], size_t num_array,
			YDataType const base_data[], YDataType d2ydx2[]) {
		StorageAndAlignedPointer<YDataType> holder_for_u;
		AllocateAndAlign < YDataType > (num_base * num_array, &holder_for_u);
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

}

#endif /* _LIBSAKURA_INTERPOLATION_SPLINE_TCC_ */
