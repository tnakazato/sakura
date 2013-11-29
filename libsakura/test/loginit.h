/*
 * loginit.h
 *
 *  Created on: 2013/11/29
 *      Author: kohji
 */

#ifndef LIBSAKURA_LOGINIT_H_
#define LIBSAKURA_LOGINIT_H_

#include <log4cxx/propertyconfigurator.h>

namespace {
struct LogInitializer {
	LogInitializer() {
		::log4cxx::PropertyConfigurator::configure("libsakura.log4j");
	}
};
static LogInitializer initializer;
}

#endif /* LIBSAKURA_LOGINIT_H_ */
