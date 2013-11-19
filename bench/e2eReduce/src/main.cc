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

#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/ArrayIO.h>

#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/ExprNode.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <libsakura/sakura.h>

#include "config_file_reader.h"
#include "option_parser.h"

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

void JobFinished(int i) {
	std::cout << "Job ";
	std::cout << i;
	std::cout << " have done.\n";
}

void ParallelJob(int job_id, unsigned int num_v, float const *v) {
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
	OptionList options;
	::ConfigFileReader::read(configuration_file, &options);

	std::string input_file;
	std::string output_file;
	OptionParser::ParseE2e(options, &input_file, &output_file);
	LOG4CXX_INFO(logger,
			"input filename=" << input_file << "\n\toutput filename=" << output_file << "\n");

	std::string sky_table;
	std::string tsys_table;
	OptionParser::ParseCalibration(options, &sky_table, &tsys_table);
	LOG4CXX_INFO(logger,
			"sky filename=" << sky_table << "\n\ttsys filename=" << tsys_table << "\n");

	int edge_channels;
	float clipping_threshold;
	bool do_clipping;
	OptionParser::ParseFlagging(options, &edge_channels, &clipping_threshold,
			&do_clipping);
	LOG4CXX_INFO(logger,
			"edge channels=" << edge_channels << "\n\tclipping_threshold=" << clipping_threshold << "\n\tdo_clipping=" << do_clipping);

	std::vector<uint64_t> line_mask;
	OptionParser::ParseBaseline(options, &line_mask);
	std::ostringstream oss;
	char separator = '[';
	oss << "line_mask=";
	for (auto i = line_mask.begin() ; i != line_mask.end() ; ++i) {
		oss << separator << *i;
		separator = ',';
	}
	oss << "]";
	LOG4CXX_INFO(logger, oss.str());

	double start_time = sakura_GetCurrentTime();
	if (input_file.size() > 0) {
		casa::Table table(casa::String(input_file), casa::Table::Old);

		casa::ROScalarColumn<unsigned int> ifno_column(table, "IFNO");
		unsigned int ifno = ifno_column(0);
		{
			std::ostringstream oss;
			oss << "Processing IFNO " << ifno << std::endl;
			LOG4CXX_INFO(logger, oss.str().c_str());
		}
		casa::Table selected_table = table(table.col("IFNO") == ifno);

		casa::ROArrayColumn<float> spectra_column(selected_table, "SPECTRA");
		unsigned int num_chan = spectra_column(0).nelements();
		auto group = xdispatch::group();
		size_t alignment = sakura_GetAlignment();
		size_t num_arena = num_chan + ((alignment - 1) / sizeof(float) + 1);
		std::vector<std::unique_ptr<float[]> > pointer_holder(20);
		for (int i = 0; i < 20; ++i) {
			pointer_holder[i].reset(new float[num_arena]);
			float *spectrum = sakura_AlignFloat(num_arena,
					pointer_holder[i].get(), num_chan);
			casa::Vector<float> spectrum_vector(casa::IPosition(1, num_chan),
					spectrum, casa::SHARE);
			spectra_column.get(i, spectrum_vector);
			float const *v = spectrum;
			if (serialize) {
				ParallelJob(i, num_chan, v);
			} else {
				group.async([=] {
					ParallelJob(i, num_chan, v);
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
