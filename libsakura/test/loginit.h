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
