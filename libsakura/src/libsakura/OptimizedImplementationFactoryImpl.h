#ifndef _OPTIMIZED_IMPLEMENTATION_FACTORY_IMPL_H
#define _OPTIMIZED_IMPLEMENTATION_FACTORY_IMPL_H

// Author: Kohji Nakamura <k.nakamura@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution

#include "libsakura/OptimizedImplementationFactory.h"

namespace sakura {
class GriddingDefault: public Gridding {
public:
	virtual void gridsd(double const xy[/*nrow*/][2],
			value_t const values/*[nrow][nvischan]*/[/*nvispol*/],
			integer nvispol, integer nvischan, bool dowt,
			flag_t const flag/*[nrow][nvischan]*/[/*nvispol*/],
			integer const rflag[/*nrow*/],
			float const weight/*[nrow]*/[/*nvischan*/], integer nrow,
			integer irow, value_t grid/*[nchan][npol][ny]*/[/*nx*/],
			float wgrid/*[nchan][npol][ny]*/[/*nx*/], integer nx, integer ny,
			integer npol, integer nchan, integer support, integer sampling,
			float const convTable[], integer const chanmap[/*nvischan*/],
			integer const polmap[/*nvispol*/],
			double sumwt/*[nchan]*/[/*npol*/]) const;
	virtual void gridsdForSpeed(double const xy[/*nrow*/][2],
			value_t const values/*[nrow][nvischan]*/[/*nvispol*/],
			integer nvispol, integer nvischan, bool dowt,
			flag_t const flag/*[nrow][nvischan]*/[/*nvispol*/],
			integer const rflag[/*nrow*/],
			float const weight/*[nrow]*/[/*nvischan*/], integer nrow,
			integer irow, value_t grid/*[ny][nx][nchan]*/[/*npol*/],
			float wgrid/*[ny][nx][nchan]*/[/*npol*/], integer nx, integer ny,
			integer npol, integer nchan, integer support, integer sampling,
			float const convTable[], integer const chanmap[/*nvischan*/],
			integer const polmap[/*nvispol*/],
			double sumwt/*[nchan]*/[/*npol*/]) const;
	virtual void transform(integer ny, integer nx, integer nchan, integer npol,
			value_t gridFrom/*[ny][nx][nchan]*/[/*npol*/],
			float wgridFrom/*[ny][nx][nchan]*/[/*npol*/],
			value_t gridTo/*[nchan][npol][ny]*/[/*nx*/],
			float wgridTo/*[nchan][npol][ny]*/[/*nx*/]) const;
};

class GriddingAfterSandyBridge: public Gridding {
public:
	virtual void gridsd(double const xy[/*nrow*/][2],
			value_t const values/*[nrow][nvischan]*/[/*nvispol*/],
			integer nvispol, integer nvischan, bool dowt,
			flag_t const flag/*[nrow][nvischan]*/[/*nvispol*/],
			integer const rflag[/*nrow*/],
			float const weight/*[nrow]*/[/*nvischan*/], integer nrow,
			integer irow, value_t grid/*[nchan][npol][ny]*/[/*nx*/],
			float wgrid/*[nchan][npol][ny]*/[/*nx*/], integer nx, integer ny,
			integer npol, integer nchan, integer support, integer sampling,
			float const convTable[], integer const chanmap[/*nvischan*/],
			integer const polmap[/*nvispol*/],
			double sumwt/*[nchan]*/[/*npol*/]) const;
	virtual void gridsdForSpeed(double const xy[/*nrow*/][2],
			value_t const values/*[nrow][nvischan]*/[/*nvispol*/],
			integer nvispol, integer nvischan, bool dowt,
			flag_t const flag/*[nrow][nvischan]*/[/*nvispol*/],
			integer const rflag[/*nrow*/],
			float const weight/*[nrow]*/[/*nvischan*/], integer nrow,
			integer irow, value_t grid/*[ny][nx][nchan]*/[/*npol*/],
			float wgrid/*[ny][nx][nchan]*/[/*npol*/], integer nx, integer ny,
			integer npol, integer nchan, integer support, integer sampling,
			float const convTable[], integer const chanmap[/*nvischan*/],
			integer const polmap[/*nvispol*/],
			double sumwt/*[nchan]*/[/*npol*/]) const;
	virtual void transform(integer ny, integer nx, integer nchan, integer npol,
			value_t gridFrom/*[ny][nx][nchan]*/[/*npol*/],
			float wgridFrom/*[ny][nx][nchan]*/[/*npol*/],
			value_t gridTo/*[nchan][npol][ny]*/[/*nx*/],
			float wgridTo/*[nchan][npol][ny]*/[/*nx*/]) const;
};

class StatisticsDefault: public Statistics {
public:
	virtual void reduce(libsakura_symbol(statistics_result) &result,
			float const *data, bool const *mask, size_t elements) const;
};

class StatisticsAfterSandyBridge: public Statistics {
public:
	virtual void reduce(libsakura_symbol(statistics_result) &result,
			float const *data, bool const *mask, size_t elements) const;
};

}

#endif /* _OPTIMIZED_IMPLEMENTATION_FACTORY_IMPL_H */
