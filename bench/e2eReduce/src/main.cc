/*
 * main.cc
 *
 *  Created on: 2013/10/21
 *      Author: kohji
 */

#include <log4cxx/logger.h>
#include <log4cxx/propertyconfigurator.h>
#include <libsakura/sakura.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>

namespace {
auto logger = log4cxx::Logger::getLogger("app");

void e2eReduce(int argc, char const* const argv[]) {
	LOG4CXX_INFO(logger, "e2eReduce");
}

}

int main(int argc, char const* const argv[]) {
	::log4cxx::PropertyConfigurator::configure("config.log4j");
	LOG4CXX_INFO(logger, "start");
	sakura_Status result = sakura_Initialize(nullptr, nullptr);
	if (result == sakura_Status_kOK) {
		try {
			e2eReduce(argc, argv);
		} catch(...) {
			LOG4CXX_ERROR(logger, "Exception raised");
			return 1;
		}
		sakura_CleanUp();
	} else {
		LOG4CXX_ERROR(logger, "Failed to initialize libsakura.");
		return 1;
	}
	return 0;
}
