#ifndef _LIBSAKURA_INTERPOLATION_LINEAR_TCC_
#define _LIBSAKURA_INTERPOLATION_LINEAR_TCC_

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
			size_t const location[], size_t offset) {
		for (size_t k = 0; k < num_location - 1; ++k) {
			size_t left_index = k + offset;
			for (size_t j = 0; j < num_array; ++j) {
				size_t offset_index_left = j * num_base + left_index;
				XDataType dydx =
						static_cast<XDataType>(base_data[offset_index_left + 1]
								- base_data[offset_index_left])
								/ (base_position[left_index + 1] - base_position[left_index]);
				YDataType base_term = base_data[offset_index_left];
				YDataType *y_work = &interpolated_data[j * num_interpolated];
				for (size_t i = location[k]; i < location[k + 1]; ++i) {
					y_work[i] = base_term
							+ static_cast<YDataType>(dydx
									* (interpolated_position[i]
											- base_position[left_index]));
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
			size_t const location[], size_t offset) {
		for (size_t k = 0; k < num_location - 1; ++k) {
			size_t left_index = k + offset;
			for (size_t i = location[k]; i < location[k + 1]; ++i) {
				YDataType fraction =
						static_cast<YDataType>((interpolated_position[i]
								- base_position[left_index])
								/ (base_position[left_index + 1] - base_position[left_index]));
				YDataType *work = &interpolated_data[num_array * i];
				YDataType const *left_value = &base_data[num_array * left_index];
				YDataType const *right_value = &base_data[num_array * (left_index + 1)];
				for (size_t j = 0; j < num_array; ++j) {
					work[j] = left_value[j]
							+ fraction * (right_value[j] - left_value[j]);
				}
			}
		}
	}
};

}

#endif /* _LIBSAKURA_INTERPOLATION_LINEAR_TCC_ */
