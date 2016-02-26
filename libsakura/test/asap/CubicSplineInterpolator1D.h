/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2016
 * National Astronomical Observatory of Japan
 * 2-21-1, Osawa, Mitaka, Tokyo, 181-8588, Japan.
 * 
 * This file is part of Sakura.
 * 
 * Sakura is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the 
 * Free Software Foundation, either version 3 of the License, or (at your 
 * option) any later version.
 * 
 * Sakura is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License 
 * along with Sakura.  If not, see <http://www.gnu.org/licenses/>.
 * @SAKURA_LICENSE_HEADER_END@
 */
//
// C++ Interface: CubicSplineInterpolator1D
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_CUBIC_SPLINE_INTERPOLATOR_1D_H
#define ASAP_CUBIC_SPLINE_INTERPOLATOR_1D_H

#include "Interpolator1D.h"

namespace asap {

/**
 * Implementation of (natural) cubic spline interpolation.
 * @author TakeshiNakazato
 */
template <class T, class U>
class CubicSplineInterpolator1D : public Interpolator1D<T, U> {
public:
  // Default constructor.
  CubicSplineInterpolator1D();

  // Destructor.
  virtual ~CubicSplineInterpolator1D();

  // Override Interpolator1D::setData.
  // @see Interpolator1D::setData
  void setData(T *x, U *y, unsigned int n);

  // Override Interpolate1D::setX.
  // @see Interpolate1D::setX()
  void setX(T *x, unsigned int n);

  // Override Interpolator1D::setY.
  // @see Interpolator1D::setY()
  void setY(U *y, unsigned int n);

  // Perform interpolation.
  // @param[in] x horizontal location where the value is evaluated 
  //              by interpolation.
  // @return interpolated value at x.
  U interpolate(T x);
private:
  // Determine second derivatives of each point based on 
  // natural cubic spline condition (second derivative at each 
  // end is zero).
  void evaly2();

  // Do interpolation using second derivatives determined by evaly2().
  // @param[in] x horizontal location where the value is evaluated 
  //              by interpolation.
  // @param[in] i location index for x.
  // @return interpolated value at x.
  U dospline(T x, unsigned int i);
  
  // Array to store second derivatives on the data points.
  U *y2_;

  // number of data points for second derivatives
  unsigned int ny2_;

  // Boolean parameter whether buffered values are effective or not.
  bool reusable_;
};

}

#include "CubicSplineInterpolator1D.tcc"

#endif
