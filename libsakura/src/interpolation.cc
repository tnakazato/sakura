#include <cassert>
#include <sstream>
#include <memory>
#include <cstdalign>
#include <utility>
#include <vector>

#include <libsakura/sakura.h>
#include <libsakura/optimized_implementation_factory_impl.h>
#include <libsakura/localdef.h>
//#include <libsakura/logger.h>
#include <libsakura/memory_manager.h>

#include "locator.tcc"
#include "interpolation_utils.tcc"
#include "interpolation.tcc"

namespace {
// a logger for this module
//auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("interpolation");

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
	size_t offset = 0;
	if (base_position_work[0] < interpolated_position_work[0]) {
		for (size_t i = 0; i < num_base - 1; ++i) {
			if (base_position_work[offset + 1]
					< interpolated_position_work[0]) {
				offset++;
			} else {
				break;
			}
		}
	}
	interpolator.Interpolate1D(num_base, base_position_work, num_array,
			base_data_work, num_interpolated, interpolated_position_work,
			interpolated_data, num_location_base, location_base, offset);

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

template<class InterpolatorSet, class XDataType, class YDataType>
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
		func = Interpolate1D<typename InterpolatorSet::NearestInterpolator,
				XDataType, YDataType>;
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kLinear):
		func = Interpolate1D<typename InterpolatorSet::LinearInterpolator,
				XDataType, YDataType>;
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kPolynomial):
		if (polynomial_order == 0) {
			// This is special case: 0-th polynomial interpolation
			// acts like nearest interpolation
			func = Interpolate1D<typename InterpolatorSet::NearestInterpolator,
					XDataType, YDataType>;
		} else {
			func = Interpolate1D<
					typename InterpolatorSet::PolynomialInterpolator, XDataType,
					YDataType>;
		}
		break;
	case LIBSAKURA_SYMBOL(InterpolationMethod_kSpline):
		func = Interpolate1D<typename InterpolatorSet::SplineInterpolator,
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
