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

class BasicOperation {
public:
	virtual ~BasicOperation() {
	}

	virtual void OperateBoolsAnd(size_t num_in, bool const in1[/*num_in*/],
			bool const in2[/*num_in*/], bool out[/*num_in*/]) const = 0;
};

template<typename DataType>
class BitOperation {
public:
	virtual ~BitOperation() {
	}

	virtual void OperateBitsAnd(DataType bit_mask, size_t num_in,
			DataType const in[/*num_in*/], bool const edit_mask[/*num_in*/],
			DataType out[/*num_in*/]) const = 0;
};

class Convolution {
public:
	virtual ~Convolution() {}
	virtual void CreateConvolve1DContext(size_t num_channel,LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
            size_t kernel_width,bool use_fft,LIBSAKURA_SYMBOL(Convole1DContext) **context) const = 0;
};

class Gridding {
public:
	typedef int32_t integer;
	typedef unsigned char flag_t;
	//typedef std::complex<float> value_t;
	typedef float value_t;

	virtual ~Gridding() {
	}

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
			) const = 0;

#if 0
	virtual void Transform(integer ny, integer nx, integer nchan, integer npol,
			value_t grid_from/*[ny][nx][nchan]*/[/*npol*/],
			float wgrid_from/*[ny][nx][nchan]*/[/*npol*/],
			value_t gridTo/*[nchan][npol][ny]*/[/*nx*/],
			float wgridTo/*[nchan][npol][ny]*/[/*nx*/]) const = 0;
#endif
};

class Interpolation {
public:
	virtual ~Interpolation() {
	}
	virtual void Interpolate1dFloatNearest(size_t num_base,
			double const x_base[/*num_base*/], float const y_base[/*num_base*/],
			size_t num_interpolated,
			double const x_interpolated[/*num_interpolated*/],
			float y_interpolated[/*num_interpolated*/]) const = 0;
	virtual void Interpolate1dFloatLinear(size_t num_base,
			double const x_base[/*num_base*/], float const y_base[/*num_base*/],
			size_t num_interpolated,
			double const x_interpolated[/*num_interpolated*/],
			float y_interpolated[/*num_interpolated*/]) const = 0;
	virtual void Interpolate1dFloatPolynomial(int polynomial_order, size_t num_base,
			double const x_base[/*num_base*/], float const y_base[/*num_base*/],
			size_t num_interpolated,
			double const x_interpolated[/*num_interpolated*/],
			float y_interpolated[/*num_interpolated*/]) const = 0;
	virtual void Interpolate1dFloatSpline(size_t num_base,
			double const x_base[/*num_base*/], float const y_base[/*num_base*/],
			size_t num_interpolated,
			double const x_interpolated[/*num_interpolated*/],
			float y_interpolated[/*num_interpolated*/]) const = 0;
protected:
	virtual int Locate(int start_position, int end_position, size_t num_base,
			double const x_base[/*num_base*/], double x_located) const = 0;
};

class Statistics {
public:
	virtual ~Statistics() {
	}
	virtual void ComputeStatistics(float const data[], bool const mask[], size_t elements,
			LIBSAKURA_SYMBOL(StatisticsResult) *result) const = 0;
};

class OptimizedImplementationFactory {
public:
	static OptimizedImplementationFactory const *GetFactory();
	virtual ~OptimizedImplementationFactory() {
	}
	virtual BasicOperation const *GetBasicOperationImpl() const = 0;
	virtual BitOperation<uint8_t> const *GetBitOperationImplUint8() const = 0;
	virtual BitOperation<uint32_t> const *GetBitOperationImplUint32() const = 0;
	virtual Convolution const *GetConvolutionImpl() const = 0;
	virtual Gridding const *GetGriddingImpl() const = 0;
	virtual Interpolation const *GetInterpolationImpl() const = 0;
	virtual Statistics const *GetStatisticsImpl() const = 0;
protected:
	OptimizedImplementationFactory() {
	}
};

} // namespace LIBSAKURA_PREFIX

#endif /* LIBSAKURA_LIBSAKURA_OPTIMIZED_IMPLEMENTATION_FACTORY_H_ */
