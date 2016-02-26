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
/*
 * loginit.h
 *
 *  Created on: 2013/11/29
 *      Author: kohji
 */

#ifndef LIBSAKURA_LOGINIT_H_
#define LIBSAKURA_LOGINIT_H_

#include <libsakura/config.h>

#if LIBSAKURA_HAS_LOG4CXX
#include <log4cxx/propertyconfigurator.h>
#endif

namespace {
struct LogInitializer {
	LogInitializer() {
#if LIBSAKURA_HAS_LOG4CXX
		::log4cxx::PropertyConfigurator::configure("libsakura.log4j");
#endif
	}
};
static LogInitializer initializer;
}

#endif /* LIBSAKURA_LOGINIT_H_ */
