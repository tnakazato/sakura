/*
 * main.cc
 *
 *  Created on: 2013/10/21
 *      Author: kohji
 */

#include <iostream>
#include <unistd.h>
#include <vector>
#include <memory>
#include <string>
#include <map>
#include <log4cxx/logger.h>
#include <log4cxx/propertyconfigurator.h>
#include <xdispatch/dispatch>

#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>

#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <libsakura/sakura.h>

#include "option_parser.h"
#include "utils.h"

namespace {
auto logger = log4cxx::Logger::getLogger("app");

inline void ExecuteCalibration() {

}

inline void ExecuteBaseline() {

}

inline void ExecuteSmoothing() {

}

inline void ExecuteNanOrInfFlag(size_t num_data, float const data[],
		bool result[]) {
	//sakura_SetFalseFloatIfNanOrInf(num_data, data, result);
}

inline void ExecuteChannelFlagging() {

}

inline void CalculateStatistics() {

}

void JobFinished(int i) {
	std::cout << "Job ";
	std::cout << i;
	std::cout << " have done.\n";
}

void ParallelJob(int job_id, unsigned int num_v, float const v[],
		uint8_t const f[], E2EOptions const option_list) {
	float sum = 0.0;
	for (unsigned int i = 0; i < num_v; ++i)
		sum += v[i];
	std::cout << "Job " << job_id << ": v[0]=" << v[0] << ", sum(v)=" << sum
			<< std::endl;
	// sleep(job_id % 3); // my job is sleeping.
	xdispatch::main_queue().sync([=]() {
		JobFinished(job_id);
	});
}

void E2eReduce(int argc, char const* const argv[]) {
	LOG4CXX_INFO(logger, "Enter: E2eReduce");
	bool const serialize = false;

	// Read configuration file
	// Default configuration file is "e2etest.config"
	std::string configuration_file = "e2etest.config";
	if (argc > 1) {
		configuration_file = argv[1];
	}

	E2EOptions options;
	OptionParser::Parse(configuration_file, &options);

	// config file summary
	LOG4CXX_INFO(logger, OptionParser::GetSummary(options));

	double start_time = sakura_GetCurrentTime();
	if (options.input_file.size() > 0) {
		AlignedArrayGenerator array_generator;

		casa::Table table = GetSelectedTable(options.input_file, options.ifno);

		casa::ROArrayColumn<float> spectra_column(table, "SPECTRA");
		casa::ROArrayColumn<unsigned char> flagtra_column(table, "FLAGTRA");
		unsigned int num_chan = spectra_column(0).nelements();

		// sky table
		float *sky_spectra;
		double *sky_time;
		size_t num_chan_sky, num_row_sky;
		GetFromCalTable(options.calibration.sky_table, options.ifno, "SPECTRA",
				&array_generator, &sky_spectra, &sky_time, &num_chan_sky,
				&num_row_sky);
		if (num_chan != num_chan_sky) {
			throw "";
		}

		// tsys table
		float *tsys;
		double *tsys_time;
		size_t num_chan_tsys, num_row_tsys;
		GetFromCalTable(options.calibration.tsys_table,
				options.calibration.tsys_ifno, "TSYS", &array_generator, &tsys,
				&tsys_time, &num_chan_tsys, &num_row_tsys);

		auto group = xdispatch::group();

		for (int i = 0; i < 20; ++i) {
			// allocate and align
			float *spectrum = array_generator.GetAlignedArray<float>(num_chan);
			unsigned char *flag =
					array_generator.GetAlignedArray<unsigned char>(num_chan);

			// get data from the table
			GetArrayCell(spectrum, i, spectra_column,
					casa::IPosition(1, num_chan));
			GetArrayCell(flag, i, flagtra_column, casa::IPosition(1, num_chan));

			// commit job
			float const *v = spectrum;
			uint8_t const *f = reinterpret_cast<uint8_t const *>(flag);
			if (serialize) {
				ParallelJob(i, num_chan, v, f, options);
			} else {
				group.async([=] {
					ParallelJob(i, num_chan, v, f, options);
				});
			}
		}
		group.wait(xdispatch::time_forever);
	}
	xdispatch::main_queue().async([=] {
		double end_time = sakura_GetCurrentTime();
		std::cout << "finished " << end_time - start_time << "secs\n";
	});
	LOG4CXX_INFO(logger, "Leave: E2eReduce");
}

void main_(int argc, char const* const argv[]) {
	LOG4CXX_INFO(logger, "start");

	sakura_Status result;
	xdispatch::main_queue().sync([&] {
		result = sakura_Initialize(nullptr, nullptr);
	});
	if (result == sakura_Status_kOK) {
		try {
			E2eReduce(argc, argv);
		} catch (casa::AipsError &e) {
			LOG4CXX_ERROR(logger, "Exception raised: " << e.getMesg());
		} catch (...) {
			LOG4CXX_ERROR(logger, "Exception raised");
		}
		xdispatch::main_queue().sync([] {
			LOG4CXX_INFO(logger, "Cleaning up libsakura");
			sakura_CleanUp();
			exit(0);
		});
	} else {
		LOG4CXX_ERROR(logger, "Failed to initialize libsakura.");
	}
	xdispatch::main_queue().sync([] {
		exit(1);
	});
}

}

int main(int argc, char const* const argv[]) {
	::log4cxx::PropertyConfigurator::configure("config.log4j");
	xdispatch::global_queue().async([=] {
		main_(argc, argv);
	});
	xdispatch::exec(); // never returns
	return 0;
}
