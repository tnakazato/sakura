#ifndef LIBSAKURA_LIBSAKURA_OPTIMIZED_IMPLEMENTATION_FACTORY_H_
#define LIBSAKURA_LIBSAKURA_OPTIMIZED_IMPLEMENTATION_FACTORY_H_

// Author: Kohji Nakamura <k.nakamura@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution

#include <complex>
#include <cstdint>
#include <cstddef>
#include <libsakura/sakura.h>

namespace LIBSAKURA_PREFIX {

template<class DataType>
class ApplyCalibration {
public:
	virtual ~ApplyCalibration() {
	}
	virtual void ApplyPositionSwitchCalibration(size_t num_scaling_factor,
			DataType const scaling_factor[/*num_scaling_factor*/],
			size_t num_data, DataType const target[/*num_data*/],
			DataType const reference[/*num_data*/],
			DataType result[/*num_data*/]) const = 0;
};

class Baseline {
public:
	virtual ~Baseline() {
	}
	virtual void SubtractBaselinePolynomial(size_t num_chan,
			float const in_data[/*num_chan*/], bool const in_mask[/*num_chan*/],
			int order, float clipping_threshold_sigma,
			int clipping_max_iteration, bool get_residual,
			float out[/*num_chan*/]) const = 0;
	virtual void GetBaselineModel(size_t num_chan, int order,
			double out[/*(order+1)*num_chan*/]) const = 0;
	virtual void DoSubtractBaseline(size_t num_chan,
			float const in_data[/*num_chan*/], bool const in_mask[/*num_chan*/],
			size_t num_model, double model_data[/*num_model * num_chan*/],
			float clipping_threshold_sigma, int clipping_max_iteration,
			bool get_residual, float out[/*num_chan*/]) const = 0;
};

template<typename DataType>
class BitOperation {
public:
	virtual ~BitOperation() {
	}

	virtual void OperateBitsAnd(DataType bit_mask, size_t num_data,
			DataType const data[/*num_data*/],
			bool const edit_mask[/*num_data*/],
			DataType result[/*num_data*/]) const = 0;
	virtual void OperateBitsConverseNonImplication(DataType bit_mask,
			size_t num_data, DataType const data[/*num_data*/],
			bool const edit_mask[/*num_data*/],
			DataType result[/*num_data*/]) const = 0;
	virtual void OperateBitsNot(size_t num_data,
			DataType const data[/*num_data*/],
			bool const edit_mask[/*num_data*/],
			DataType result[/*num_data*/]) const = 0;
	virtual void OperateBitsOr(DataType bit_mask, size_t num_data,
			DataType const data[/*num_data*/],
			bool const edit_mask[/*num_data*/],
			DataType result[/*num_data*/]) const = 0;
};

template<typename DataType>
class BoolFilterCollection {
public:
	virtual ~BoolFilterCollection() {
	}

	virtual void SetTrueInRangesInclusive(size_t num_data,
			DataType const data[/*num_data*/], size_t num_condition,
			DataType const lower_bounds[/*num_condition*/],
			DataType const upper_bounds[/*num_condition*/],
			bool result[/*num_data*/]) const = 0;
	virtual void ToBool(size_t num_data, DataType const data[/*num_data*/],
	bool result[/*num_data*/]) const = 0;
	virtual void InvertBool(size_t num_data, bool const data[/*num_data*/],
	bool result[/*num_data*/]) const = 0;
};

class Convolution {
public:
	virtual ~Convolution() {
	}
	virtual void CreateConvolve1DContext(size_t num_channel,
	LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type, size_t kernel_width,
	bool use_fft, LIBSAKURA_SYMBOL(Convolve1DContext) **context) const = 0;
	virtual void Convolve1D(LIBSAKURA_SYMBOL(Convolve1DContext) **context,
			float input_spectrum[/*num_channels*/],
			bool const input_flag[/*num_channels*/],
			float output_spectrum[/*num_channels*/]) const = 0;
	virtual void DestroyConvolve1DContext(
	LIBSAKURA_SYMBOL(Convolve1DContext) *context) const = 0;
};

class Gridding {
public:
	typedef int32_t integer;
	typedef unsigned char flag_t;
	//typedef std::complex<float> value_t;
	typedef float value_t;

	virtual ~Gridding() {
	}

	virtual void GridConvolving(size_t num_spectra, size_t start_spectrum,
			size_t end_spectrum, bool const spectrum_mask[/*num_spectra*/],
			double const x[/*num_spectra*/], double const y[/*num_spectra*/],
			size_t support, size_t sampling, size_t num_polarization,
			uint32_t const polarization_map[/*num_polarization*/],
			size_t num_channels, uint32_t const channel_map[/*num_channels*/],
			bool const mask/*[num_spectra][num_polarization]*/[/*num_channels*/],
			float const value/*[num_spectra][num_polarization]*/[/*num_channels*/],
			float const weight/*[num_spectra]*/[/*num_channels*/],
			bool do_weight,
			size_t num_convolution_table/*= ceil(sqrt(2.)*(support+1)*sampling)*/,
			float const convolution_table[/*num_convolution_table*/],
			size_t num_polarization_for_grid, size_t num_channels_for_grid,
			size_t width, size_t height,
			double weight_sum/*[num_polarization_for_grid]*/[/*num_channels_for_grid*/],
			float weight_of_grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/],
			float grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/]) const = 0;

#if 0
	virtual void Transform(integer ny, integer nx, integer nchan, integer npol,
			value_t grid_from/*[ny][nx][nchan]*/[/*npol*/],
			float wgrid_from/*[ny][nx][nchan]*/[/*npol*/],
			value_t gridTo/*[nchan][npol][ny]*/[/*nx*/],
			float wgridTo/*[nchan][npol][ny]*/[/*nx*/]) const = 0;
#endif
};

template<class XDataType, class YDataType>
class Interpolation {
public:
	virtual ~Interpolation() {
	}
	virtual LIBSAKURA_SYMBOL(Status) Interpolate1D(
	LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
			uint8_t polynomial_order, size_t num_base,
			XDataType const x_base[/*num_base*/], size_t num_array,
			YDataType const y_base[/*num_base * num_base_array*/],
			size_t num_interpolated,
			XDataType const x_interpolated[/*num_interpolated*/],
			YDataType y_interpolated[/*num_interpolated * num_base_array*/]) const = 0;
};

class LogicalOperation {
public:
	virtual ~LogicalOperation() {
	}

	virtual void OperateLogicalAnd(size_t num_in, bool const in1[/*num_in*/],
	bool const in2[/*num_in*/], bool out[/*num_in*/]) const = 0;
};

class NumericOperation {
public:
	virtual ~NumericOperation() {
	}

	virtual void OperateFloatSubtraction(size_t num_in,
			float const in1[/*num_in*/], float const in2[/*num_in*/],
			float out[/*num_in*/]) const = 0;
	virtual void GetBestFitModel(size_t num_in, float const in_data[/*num_in*/],
	bool const in_mask[/*num_in*/], size_t num_model,
			double const model[/*num_model * num_in*/],
			float out[/*num_in*/]) const = 0;
	virtual void GetLeastSquareMatrix(size_t num_in,
			float const in_data[/*num_in*/], bool const in_mask[/*num_in*/],
			size_t num_model, double const model[/*num_model * num_in*/],
			double out[/*num_model * num_model*/],
			double out_vector[/*num_model*/]) const = 0;
	virtual void SolveSimultaneousEquationsByLU(size_t num_eqn,
			double const lsq_matrix0[/*num_eqn * num_eqn*/],
			double const lsq_vector0[/*num_eqn*/],
			double out[/*num_eqn*/]) const = 0;
	virtual void DoGetBestFitModel(size_t num_chan, size_t num_eqn,
			double const model[/*num_eqn * num_chan*/],
			double const coeff[/*num_eqn*/], float out[/*num_in*/]) const = 0;
};

class Statistics {
public:
	virtual ~Statistics() {
	}
	virtual void ComputeStatistics(float const data[], bool const mask[],
			size_t elements,
			LIBSAKURA_SYMBOL(StatisticsResult) *result) const = 0;
};

class OptimizedImplementationFactory {
public:
	/**
	 * @~
	 * MT-unsafe
	 */
	static void InitializeFactory(char const *simd_spec);
	/**
	 * @~
	 * MT-unsafe
	 */
	static void CleanUpFactory();
	/**
	 * @~
	 * MT-safe
	 */
	static OptimizedImplementationFactory const *GetFactory();
	/**
	 * @~
	 * MT-unsafe
	 */
	virtual ~OptimizedImplementationFactory() {
	}
	virtual char const *GetName() const = 0;
	virtual ApplyCalibration<float> const *GetApplyCalibrationImpl() const = 0;
	virtual Baseline const *GetBaselineImpl() const = 0;
	virtual BitOperation<uint8_t> const *GetBitOperationImplUint8() const = 0;
	virtual BitOperation<uint32_t> const *GetBitOperationImplUint32() const = 0;
	virtual BoolFilterCollection<float> const *GetBoolFilterCollectionImplFloat() const = 0;
	virtual BoolFilterCollection<int> const *GetBoolFilterCollectionImplInt() const = 0;
	virtual BoolFilterCollection<uint8_t> const *GetBoolFilterCollectionImplUint8() const = 0;
	virtual BoolFilterCollection<uint32_t> const *GetBoolFilterCollectionImplUint32() const = 0;
	virtual Convolution const *GetConvolutionImpl() const = 0;
	virtual Gridding const *GetGriddingImpl() const = 0;
	virtual Interpolation<double, float> const *GetInterpolationImpl() const = 0;
	virtual LogicalOperation const *GetLogicalOperationImpl() const = 0;
	virtual NumericOperation const *GetNumericOperationImpl() const = 0;
	virtual Statistics const *GetStatisticsImpl() const = 0;
protected:
	static OptimizedImplementationFactory const *factory_;
	OptimizedImplementationFactory() {
	}
};

} // namespace LIBSAKURA_PREFIX

#endif /* LIBSAKURA_LIBSAKURA_OPTIMIZED_IMPLEMENTATION_FACTORY_H_ */
