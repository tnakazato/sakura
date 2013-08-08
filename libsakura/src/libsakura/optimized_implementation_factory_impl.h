#ifndef LIBSAKURA_LIBSAKURA_OPTIMIZED_IMPLEMENTATION_FACTORY_IMPL_H_
#define LIBSAKURA_LIBSAKURA_OPTIMIZED_IMPLEMENTATION_FACTORY_IMPL_H_

// Author: Kohji Nakamura <k.nakamura@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution

#include <libsakura/optimized_implementation_factory.h>

namespace LIBSAKURA_PREFIX {
class GriddingDefault: public Gridding {
public:
	virtual void GridConvolving(size_t num_spectra,
			size_t start_spectrum, size_t end_spectrum,
			bool const spectrum_mask[/*num_spectra*/],
			double const x[/*num_spectra*/],
			double const y[/*num_spectra*/],
			size_t support, size_t sampling,
			size_t num_polarization,
			uint32_t const polarization_map[/*num_polarization*/],
			size_t num_channels,
			uint32_t const channel_map[/*num_channels*/],
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
			float grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/]
			) const;
};

class GriddingAfterSandyBridge: public Gridding {
public:
	virtual void GridConvolving(size_t num_spectra,
			size_t start_spectrum, size_t end_spectrum,
			bool const spectrum_mask[/*num_spectra*/],
			double const x[/*num_spectra*/],
			double const y[/*num_spectra*/],
			size_t support, size_t sampling,
			size_t num_polarization,
			uint32_t const polarization_map[/*num_polarization*/],
			size_t num_channels,
			uint32_t const channel_map[/*num_channels*/],
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
			float grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/]
			) const;
};

class StatisticsDefault: public Statistics {
public:
	virtual void ComputeStatistics(float const data[], bool const mask[], size_t elements,
			LIBSAKURA_SYMBOL(StatisticsResult) *result) const;
};

class StatisticsAfterSandyBridge: public Statistics {
public:
	virtual void ComputeStatistics(float const data[], bool const mask[], size_t elements,
			LIBSAKURA_SYMBOL(StatisticsResult) *result) const;
};

template<typename DataType>
class BitOperationDefault: public BitOperation<DataType> {
public:
	virtual ~BitOperationDefault(){}
	virtual void OperateBitsAnd(DataType bit_mask, size_t num_in,
			DataType const in[/*num_in*/], bool const edit_mask[/*num_in*/],
			DataType out[]) const;
//	virtual void OperateBitsAnd(uint8_t bit_mask, size_t num_in,
//			uint8_t const in[/*num_in*/], bool const edit_mask[/*num_in*/],
//			uint8_t out[/*num_in*/]) const;
//	virtual void OperateBitsAnd(uint32_t bit_mask, size_t num_in,
//			uint32_t const in[/*num_in*/], bool const edit_mask[/*num_in*/],
//			uint32_t out[/*num_in*/]) const;
};

template<typename DataType>
class BitOperationAfterSandyBridge: public BitOperation<DataType> {
public:
	virtual ~BitOperationAfterSandyBridge(){}
	virtual void OperateBitsAnd(DataType bit_mask, size_t num_in,
			DataType const in[/*num_in*/], bool const edit_mask[/*num_in*/],
			DataType out[/*num_in*/]) const;
//	virtual void OperateBitsAnd(uint8_t bit_mask, size_t num_in,
//			uint8_t const in[/*num_in*/], bool const edit_mask[/*num_in*/],
//			uint8_t out[/*num_in*/]) const;
//	virtual void OperateBitsAnd(uint32_t bit_mask, size_t num_in,
//			uint32_t const in[/*num_in*/], bool const edit_mask[/*num_in*/],
//			uint32_t out[/*num_in*/]) const;
};

class InterpolationDefault: public Interpolation {
public:
	virtual void Interpolate1dFloatNearest(size_t num_base,
			double const x_base[/*num_base*/], float const y_base[/*num_base*/],
			size_t num_interpolated,
			double const x_interpolated[/*num_interpolated*/],
			float y_interpolated[/*num_interpolated*/]) const;
	virtual void Interpolate1dFloatLinear(size_t num_base,
			double const x_base[/*num_base*/], float const y_base[/*num_base*/],
			size_t num_interpolated,
			double const x_interpolated[/*num_interpolated*/],
			float y_interpolated[/*num_interpolated*/]) const;
	virtual void Interpolate1dFloatPolynomial(int polynomial_order, size_t num_base,
			double const x_base[/*num_base*/], float const y_base[/*num_base*/],
			size_t num_interpolated,
			double const x_interpolated[/*num_interpolated*/],
			float y_interpolated[/*num_interpolated*/]) const;
	virtual void Interpolate1dFloatSpline(size_t num_base,
			double const x_base[/*num_base*/], float const y_base[/*num_base*/],
			size_t num_interpolated,
			double const x_interpolated[/*num_interpolated*/],
			float y_interpolated[/*num_interpolated*/]) const;
protected:
	virtual int Locate(int start_position, int end_position, size_t num_base,
			double const x_base[/*num_base*/], double x_located) const;
};

class InterpolationAfterSandyBridge: public Interpolation {
public:
	virtual void Interpolate1dFloatNearest(size_t num_base,
			double const x_base[/*num_base*/], float const y_base[/*num_base*/],
			size_t num_interpolated,
			double const x_interpolated[/*num_interpolated*/],
			float y_interpolated[/*num_interpolated*/]) const;
	virtual void Interpolate1dFloatLinear(size_t num_base,
			double const x_base[/*num_base*/], float const y_base[/*num_base*/],
			size_t num_interpolated,
			double const x_interpolated[/*num_interpolated*/],
			float y_interpolated[/*num_interpolated*/]) const;
	virtual void Interpolate1dFloatPolynomial(int polynomial_order, size_t num_base,
			double const x_base[/*num_base*/], float const y_base[/*num_base*/],
			size_t num_interpolated,
			double const x_interpolated[/*num_interpolated*/],
			float y_interpolated[/*num_interpolated*/]) const;
	virtual void Interpolate1dFloatSpline(size_t num_base,
			double const x_base[/*num_base*/], float const y_base[/*num_base*/],
			size_t num_interpolated,
			double const x_interpolated[/*num_interpolated*/],
			float y_interpolated[/*num_interpolated*/]) const;
protected:
	virtual int Locate(int start_position, int end_position, size_t num_base,
			double const x_base[/*num_base*/], double x_located) const;
};

} // namespace LIBSAKURA_PREFIX

#endif /* LIBSAKURA_LIBSAKURA_OPTIMIZED_IMPLEMENTATION_FACTORY_IMPL_H_ */
