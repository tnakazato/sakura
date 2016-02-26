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
// C++ Implementation: BisectionLocator
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <assert.h>

#include "BisectionLocator.h"

namespace asap {
template <class T> BisectionLocator<T>::BisectionLocator()
  : Locator<T>()
{
}


template <class T> BisectionLocator<T>::BisectionLocator(T *v, unsigned int n, bool copystorage)
  : Locator<T>(v, n, copystorage)
{}

template <class T> BisectionLocator<T>::~BisectionLocator()
{}

template <class T> unsigned int BisectionLocator<T>::locate(T x)
{
  if (this->n_ == 1)
    return 0;

  return this->bisection(x, 0, this->n_);
}

}
