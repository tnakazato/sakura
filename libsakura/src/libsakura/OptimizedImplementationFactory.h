#ifndef _OPTIMIZED_IMPLEMENTATION_FACTORY_H
#define _OPTIMIZED_IMPLEMENTATION_FACTORY_H

// Author: Kohji Nakamura <k.nakamura@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution

#include <complex>
#include <stdint.h>

namespace sakura {
  class Gridding {
  public:
    typedef int32_t integer;
    typedef unsigned char flag_t;
    //typedef std::complex<float> value_t;
    typedef float value_t;
    virtual void gridsd(double const xy[/*nrow*/][2],
			value_t const values/*[nrow][nvischan]*/[/*nvispol*/],
			integer nvispol, integer nvischan, 
			bool dowt,
			flag_t const flag/*[nrow][nvischan]*/[/*nvispol*/],
			integer const rflag[/*nrow*/],
			float const weight/*[nrow]*/[/*nvischan*/],
			integer nrow, integer irow,
			value_t grid/*[nchan][npol][ny]*/[/*nx*/],
			float wgrid/*[nchan][npol][ny]*/[/*nx*/],
			integer nx, integer ny, integer npol, integer nchan, 
			integer support, integer sampling,
			float const convTable[],
			integer const chanmap[/*nvischan*/],
			integer const polmap[/*nvispol*/],
			double sumwt/*[nchan]*/[/*npol*/]) = 0;
    virtual void gridsdForSpeed(double const xy[/*nrow*/][2],
			value_t const values/*[nrow][nvischan]*/[/*nvispol*/],
			integer nvispol, integer nvischan, 
			bool dowt,
			flag_t const flag/*[nrow][nvischan]*/[/*nvispol*/],
			integer const rflag[/*nrow*/],
			float const weight/*[nrow]*/[/*nvischan*/],
			integer nrow, integer irow,
			value_t grid/*[ny][nx][nchan]*/[/*npol*/],
			float wgrid/*[ny][nx][nchan]*/[/*npol*/],
			integer nx, integer ny, integer npol, integer nchan, 
			integer support, integer sampling,
			float const convTable[],
			integer const chanmap[/*nvischan*/],
			integer const polmap[/*nvispol*/],
			double sumwt/*[nchan]*/[/*npol*/]) = 0;
    virtual void transform(
			integer ny, integer nx, integer nchan,  integer npol,
			value_t gridFrom/*[ny][nx][nchan]*/[/*npol*/],
			float wgridFrom/*[ny][nx][nchan]*/[/*npol*/],
			value_t gridTo/*[nchan][npol][ny]*/[/*nx*/],
			float wgridTo/*[nchan][npol][ny]*/[/*nx*/]) = 0;
    virtual ~Gridding() {}
  };

  class OptimizedImplementationFactory {
  protected:
    OptimizedImplementationFactory() {}
  public:
    virtual ~OptimizedImplementationFactory() {}
    virtual Gridding *getGriddingImpl() = 0;
    static OptimizedImplementationFactory *getFactory();
  };
}

#endif /* _OPTIMIZED_IMPLEMENTATION_FACTORY_H */
