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

class Gridding {
public:
	typedef int32_t integer;
	typedef unsigned char flag_t;
	//typedef std::complex<float> value_t;
	typedef float value_t;

	virtual ~Gridding() {
	}

	virtual void GridSd(double const xy[/*nrow*/][2],
			value_t const values/*[nrow][nvischan]*/[/*nvispol*/],
			integer nvispol, integer nvischan, bool dowt,
			flag_t const flag/*[nrow][nvischan]*/[/*nvispol*/],
			integer const rflag[/*nrow*/],
			float const weight/*[nrow]*/[/*nvischan*/], integer nrow,
			integer irow, value_t grid/*[nchan][npol][ny]*/[/*nx*/],
			float wgrid/*[nchan][npol][ny]*/[/*nx*/], integer nx, integer ny,
			integer npol, integer nchan, integer support, integer sampling,
			float const conv_table[], integer const chanmap[/*nvischan*/],
			integer const polmap[/*nvispol*/],
			double sumwt/*[nchan]*/[/*npol*/]) const = 0;
	virtual void GridSdForSpeed(double const xy[/*nrow*/][2],
			value_t const values/*[nrow][nvischan]*/[/*nvispol*/],
			integer nvispol, integer nvischan, bool dowt,
			flag_t const flag/*[nrow][nvischan]*/[/*nvispol*/],
			integer const rflag[/*nrow*/],
			float const weight/*[nrow]*/[/*nvischan*/], integer nrow,
			integer irow, value_t grid/*[ny][nx][nchan]*/[/*npol*/],
			float wgrid/*[ny][nx][nchan]*/[/*npol*/], integer nx, integer ny,
			integer npol, integer nchan, integer support, integer sampling,
			float const conv_table[], integer const chanmap[/*nvischan*/],
			integer const polmap[/*nvispol*/],
			double sumwt/*[nchan]*/[/*npol*/]) const = 0;
	virtual void Transform(integer ny, integer nx, integer nchan, integer npol,
			value_t grid_from/*[ny][nx][nchan]*/[/*npol*/],
			float wgrid_from/*[ny][nx][nchan]*/[/*npol*/],
			value_t gridTo/*[nchan][npol][ny]*/[/*nx*/],
			float wgridTo/*[nchan][npol][ny]*/[/*nx*/]) const = 0;
};

class Statistics {
public:
	virtual ~Statistics() {
	}
	virtual void ComputeStatistics(float const data[], bool const mask[], size_t elements,
			LIBSAKURA_SYMBOL(StatisticsResult) &result) const = 0;
};

class OptimizedImplementationFactory {
public:
	virtual ~OptimizedImplementationFactory() {
	}
	virtual Gridding const *GetGriddingImpl() const = 0;
	virtual Statistics const *GetStatisticsImpl() const = 0;
	static OptimizedImplementationFactory const *GetFactory();
protected:
	OptimizedImplementationFactory() {
	}
};

} // namespace LIBSAKURA_PREFIX

#endif /* LIBSAKURA_LIBSAKURA_OPTIMIZED_IMPLEMENTATION_FACTORY_H_ */
