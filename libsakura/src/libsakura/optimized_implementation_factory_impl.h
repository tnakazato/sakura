#ifndef LIBSAKURA_LIBSAKURA_OPTIMIZED_IMPLEMENTATION_FACTORY_IMPL_H_
#define LIBSAKURA_LIBSAKURA_OPTIMIZED_IMPLEMENTATION_FACTORY_IMPL_H_

// Author: Kohji Nakamura <k.nakamura@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution

#include <libsakura/optimized_implementation_factory.h>

namespace LIBSAKURA_PREFIX {

class BaselineDefault: public Baseline {
public:
	virtual ~BaselineDefault() {
	}
	virtual void SubtractBaselinePolynomial(size_t num_chan,
			float const in_data[/*num_chan*/],
			bool const in_mask[/*num_chan*/], int order,
			float clipping_threshold_sigma, int clipping_max_iteration,
			bool get_residual, float out[/*num_chan*/]) const;
	virtual void GetBaselineModel(size_t num_chan, int order,
			double out[/*(order+1)*num_chan*/]) const;
	virtual void DoSubtractBaseline(size_t num_chan,
			float const in_data[/*num_chan*/],
			bool const in_mask[/*num_chan*/], size_t num_model,
			double model_data[/*num_model * num_chan*/],
			float clipping_threshold_sigma, int clipping_max_iteration,
			bool get_residual, float out[/*num_chan*/]) const;
};

class BaselineAfterSandyBridge: public Baseline {
public:
	virtual ~BaselineAfterSandyBridge() {
	}
	virtual void SubtractBaselinePolynomial(size_t num_chan,
			float const in_data[/*num_chan*/],
			bool const in_mask[/*num_chan*/], int order,
			float clipping_threshold_sigma, int clipping_max_iteration,
			bool get_residual, float out[/*num_chan*/]) const;
	virtual void GetBaselineModel(size_t num_chan, int order,
			double out[/*(order+1)*num_chan*/]) const;
	virtual void DoSubtractBaseline(size_t num_chan,
			float const in_data[/*num_chan*/],
			bool const in_mask[/*num_chan*/], size_t num_model,
			double model_data[/*num_model * num_chan*/],
			float clipping_threshold_sigma, int clipping_max_iteration,
			bool get_residual, float out[/*num_chan*/]) const;
};

template<typename DataType>
class BitOperationDefault: public BitOperation<DataType> {
public:
	virtual ~BitOperationDefault() {
	}
	virtual void OperateBitsAnd(DataType bit_mask, size_t num_data,
			DataType const data[/*num_data*/], bool const edit_mask[/*num_data*/],
			DataType result[/*num_data*/]) const;
};

template<typename DataType>
class BitOperationAfterSandyBridge: public BitOperation<DataType> {
public:
	virtual ~BitOperationAfterSandyBridge() {
	}
	virtual void OperateBitsAnd(DataType bit_mask, size_t num_data,
			DataType const data[/*num_data*/], bool const edit_mask[/*num_data*/],
			DataType result[/*num_data*/]) const;
};

template<typename DataType>
class BoolFilterCollectionDefault: public BoolFilterCollection<DataType> {
public:
	virtual ~BoolFilterCollectionDefault() {
	}
	virtual void SetTrueInRangesInclusive(size_t num_data,
			DataType const data[/*num_data*/], size_t num_condition,
			DataType const lower_bounds[/*num_condition*/],
			DataType const upper_bounds[/*num_condition*/],
			bool result[/*num_data*/]) const;
	virtual void ToBool(size_t num_data, DataType const data[/*num_data*/],
			bool result[/*num_data*/]) const;
	virtual void InvertBool(size_t num_data, bool const data[/*num_data*/],
			bool result[/*num_data*/]) const;
};

template<typename DataType>
class BoolFilterCollectionAfterSandyBridge: public BoolFilterCollection<DataType> {
public:
	virtual ~BoolFilterCollectionAfterSandyBridge() {
	}
	virtual void SetTrueInRangesInclusive(size_t num_data,
			DataType const data[/*num_data*/], size_t num_condition,
			DataType const lower_bounds[/*num_condition*/],
			DataType const upper_bounds[/*num_condition*/],
			bool result[/*num_data*/]) const;
	virtual void ToBool(size_t num_data, DataType const data[/*num_data*/],
			bool result[/*num_data*/]) const;
	virtual void InvertBool(size_t num_data, bool const data[/*num_data*/],
			bool result[/*num_data*/]) const;
};

class ConvolutionDefault: public Convolution {
public:
	virtual ~ConvolutionDefault() {
	}
	virtual void CreateConvolve1DContext(size_t num_channel,
			LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type, size_t kernel_width,
			bool use_fft,LIBSAKURA_SYMBOL(Convolve1DContext) **context) const;
	virtual void Convolve1D(LIBSAKURA_SYMBOL(Convolve1DContext) **context,
			float input_spectrum[/*num_channels*/],bool const input_flag[/*num_channels*/],
			float output_spectrum[/*num_channels*/]) const;
	virtual void DestroyConvolve1DContext(
			LIBSAKURA_SYMBOL(Convolve1DContext) *context) const;
};

class ConvolutionAfterSandyBridge: public Convolution {
public:
	virtual ~ConvolutionAfterSandyBridge() {
	}
	virtual void CreateConvolve1DContext(size_t num_channel,
			LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type, size_t kernel_width,
			bool use_fft,LIBSAKURA_SYMBOL(Convolve1DContext) **context) const;
	virtual void Convolve1D(LIBSAKURA_SYMBOL(Convolve1DContext) **context,
			float input_spectrum[/*num_channels*/],bool const input_flag[/*num_channels*/],
			float output_spectrum[/*num_channels*/]) const;
	virtual void DestroyConvolve1DContext(
			LIBSAKURA_SYMBOL(Convolve1DContext) *context) const;
};

class GriddingDefault: public Gridding {
public:
	virtual void GridConvolving(size_t num_spectra, size_t start_spectrum,
			size_t end_spectrum,
			bool const spectrum_mask[/*num_spectra*/],
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
			float grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/]) const;
};

class GriddingAfterSandyBridge: public Gridding {
public:
	virtual void GridConvolving(size_t num_spectra, size_t start_spectrum,
			size_t end_spectrum,
			bool const spectrum_mask[/*num_spectra*/],
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
			float grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/]) const;
};

template<class XDataType, class YDataType>
class InterpolationImpl: public Interpolation<XDataType, YDataType> {
public:
	virtual ~InterpolationImpl() {
	}
protected:
	virtual size_t Locate(size_t start_position, size_t end_position,
			size_t num_base, XDataType const x_base[/*num_base*/],
			XDataType x_located) const;
};

template<class XDataType, class YDataType>
class InterpolationDefault: public InterpolationImpl<XDataType, YDataType> {
public:
	virtual ~InterpolationDefault() {
	}
	virtual LIBSAKURA_SYMBOL(Status) Interpolate1D(
	LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
			int polynomial_order, size_t num_base,
			XDataType const x_base[/*num_base*/], size_t num_array,
			YDataType const y_base[/*num_base * num_base_array*/],
			size_t num_interpolated,
			XDataType const x_interpolated[/*num_interpolated*/],
			YDataType y_interpolated[/*num_interpolated * num_base_array*/]) const;
};

template<class XDataType, class YDataType>
class InterpolationAfterSandyBridge: public InterpolationImpl<XDataType, YDataType> {
public:
	virtual ~InterpolationAfterSandyBridge() {
	}
	virtual LIBSAKURA_SYMBOL(Status) Interpolate1D(
	LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
			int polynomial_order, size_t num_base,
			XDataType const x_base[/*num_base*/], size_t num_array,
			YDataType const y_base[/*num_base * num_base_array*/],
			size_t num_interpolated,
			XDataType const x_interpolated[/*num_interpolated*/],
			YDataType y_interpolated[/*num_interpolated * num_base_array*/]) const;
};

class LogicalOperationDefault: public LogicalOperation {
public:
	virtual ~LogicalOperationDefault() {
	}
	virtual void OperateLogicalAnd(size_t num_in, bool const in1[/*num_in*/],
	bool const in2[/*num_in*/], bool out[/*num_in*/]) const;
};

class LogicalOperationAfterSandyBridge: public LogicalOperation {
public:
	virtual ~LogicalOperationAfterSandyBridge() {
	}
	virtual void OperateLogicalAnd(size_t num_in, bool const in1[/*num_in*/],
	bool const in2[/*num_in*/], bool out[/*num_in*/]) const;
};

class NumericOperationDefault: public NumericOperation {
public:
	virtual ~NumericOperationDefault() {
	}
	virtual void OperateFloatSubtraction(size_t num_in,
			float const in1[/*num_in*/], float const in2[/*num_in*/],
			float out[/*num_in*/]) const;
	virtual void GetBestFitModel(size_t num_in, float const in_data[/*num_in*/],
	bool const in_mask[/*num_in*/], size_t num_model,
			double const model[/*num_model * num_in*/],
			float out[/*num_in*/]) const;
	virtual void GetLeastSquareMatrix(size_t num_in,
			float const in_data[/*num_in*/],
			bool const in_mask[/*num_in*/], size_t num_model,
			double const model[/*num_model * num_in*/],
			double out[/*num_model * num_model*/],
			double out_vector[/*num_model*/]) const;
	virtual void SolveSimultaneousEquationsByLU(size_t num_eqn,
			double const lsq_matrix0[/*num_eqn * num_eqn*/],
			double const lsq_vector0[/*num_eqn*/],
			double out[/*num_eqn*/]) const;
	virtual void DoGetBestFitModel(size_t num_chan, size_t num_eqn,
			double const model[/*num_eqn * num_in*/],
			double const coeff[/*num_eqn*/], float out[/*num_in*/]) const;
};

class NumericOperationAfterSandyBridge: public NumericOperation {
public:
	virtual ~NumericOperationAfterSandyBridge() {
	}
	virtual void OperateFloatSubtraction(size_t num_in,
			float const in1[/*num_in*/], float const in2[/*num_in*/],
			float out[/*num_in*/]) const;
	virtual void GetBestFitModel(size_t num_in, float const in_data[/*num_in*/],
	bool const in_mask[/*num_in*/], size_t num_model,
			double const model[/*num_model * num_in*/],
			float out[/*num_in*/]) const;
protected:
	virtual void GetLeastSquareMatrix(size_t num_in,
			float const in_data[/*num_in*/],
			bool const in_mask[/*num_in*/], size_t num_model,
			double const model[/*num_model * num_in*/],
			double out[/*num_model * num_model*/],
			double out_vector[/*num_model*/]) const;
	virtual void SolveSimultaneousEquationsByLU(size_t num_eqn,
			double const lsq_matrix0[/*num_eqn * num_eqn*/],
			double const lsq_vector0[/*num_eqn*/],
			double out[/*num_eqn*/]) const;
	virtual void DoGetBestFitModel(size_t num_chan, size_t num_eqn,
			double const model[/*num_eqn * num_chan*/],
			double const coeff[/*num_eqn*/], float out[/*num_in*/]) const;
};

class StatisticsDefault: public Statistics {
public:
	virtual void ComputeStatistics(float const data[], bool const mask[],
			size_t elements,
			LIBSAKURA_SYMBOL(StatisticsResult) *result) const;
};

class StatisticsAfterSandyBridge: public Statistics {
public:
	virtual void ComputeStatistics(float const data[], bool const mask[],
			size_t elements,
			LIBSAKURA_SYMBOL(StatisticsResult) *result) const;
};

} // namespace LIBSAKURA_PREFIX

#endif /* LIBSAKURA_LIBSAKURA_OPTIMIZED_IMPLEMENTATION_FACTORY_IMPL_H_ */
