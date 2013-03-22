#ifndef LIBSAKURA_LIBSAKURA_OPTIMIZED_IMPLEMENTATION_FACTORY_IMPL_H_
#define LIBSAKURA_LIBSAKURA_OPTIMIZED_IMPLEMENTATION_FACTORY_IMPL_H_

// Author: Kohji Nakamura <k.nakamura@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution

#include <libsakura/optimized_implementation_factory.h>

namespace LIBSAKURA_PREFIX {
class GriddingDefault: public Gridding {
public:
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
			double sumwt/*[nchan]*/[/*npol*/]) const;
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
			double sumwt/*[nchan]*/[/*npol*/]) const;
	virtual void Transform(integer ny, integer nx, integer nchan, integer npol,
			value_t grid_from/*[ny][nx][nchan]*/[/*npol*/],
			float wgrid_from/*[ny][nx][nchan]*/[/*npol*/],
			value_t gridTo/*[nchan][npol][ny]*/[/*nx*/],
			float wgridTo/*[nchan][npol][ny]*/[/*nx*/]) const;
};

class GriddingAfterSandyBridge: public Gridding {
public:
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
			double sumwt/*[nchan]*/[/*npol*/]) const;
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
			double sumwt/*[nchan]*/[/*npol*/]) const;
	virtual void Transform(integer ny, integer nx, integer nchan, integer npol,
			value_t grid_from/*[ny][nx][nchan]*/[/*npol*/],
			float wgrid_from/*[ny][nx][nchan]*/[/*npol*/],
			value_t gridTo/*[nchan][npol][ny]*/[/*nx*/],
			float wgridTo/*[nchan][npol][ny]*/[/*nx*/]) const;
};

class StatisticsDefault: public Statistics {
public:
	virtual void ComputeStatistics(float const data[], bool const mask[], size_t elements,
			LIBSAKURA_SYMBOL(StatisticsResult) &result) const;
};

class StatisticsAfterSandyBridge: public Statistics {
public:
	virtual void ComputeStatistics(float const data[], bool const mask[], size_t elements,
			LIBSAKURA_SYMBOL(StatisticsResult) &result) const;
};

} // namespace LIBSAKURA_PREFIX

#endif /* LIBSAKURA_LIBSAKURA_OPTIMIZED_IMPLEMENTATION_FACTORY_IMPL_H_ */
