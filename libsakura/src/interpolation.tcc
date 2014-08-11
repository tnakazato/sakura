/*
 * @SAKURA_LICENSE_HEADER_START@
 * @SAKURA_LICENSE_HEADER_END@
 */
#ifndef _LIBSAKURA_INTERPOLATION_TCC_
#define _LIBSAKURA_INTERPOLATION_TCC_

#include "interpolation_utils.tcc"
#include "interpolation_nearest.tcc"
#include "interpolation_linear.tcc"
#include "interpolation_polynomial.tcc"
#include "interpolation_spline.tcc"

namespace {

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

}

#endif /* _LIBSAKURA_INTERPOLATION_TCC_ */
