#ifndef _LIBSAKURA_INTERPOLATION_LINEAR_H_
#define _LIBSAKURA_INTERPOLATION_LINEAR_H_

namespace {

// TODO: documentation must be added
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

}

#endif /* _LIBSAKURA_INTERPOLATION_LINEAR_H_ */
