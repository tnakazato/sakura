/*
 * main.cc
 *
 *  Created on: 2013/10/21
 *      Author: kohji
 */

#include <getopt.h>
#include <iostream>
#include <unistd.h>
#include <vector>
#include <memory>
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

namespace {
auto logger = log4cxx::Logger::getLogger("app");

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

void PrintUsage(char const *name) {
	std::ostringstream oss;
	oss << name << "\t-i|--input-file input_filename\n"
			<< "\t[-o|--output-file output_filename]\n"
			<< "\t-s|--sky-file skycal_table\n"
			<< "\t[-t|--tsys-file tsyscal_table]\n";
	std::cout << oss.str();
}

void AnalyseArguments(int argc, char const* const argv[],
		casa::String *input_file, casa::String *output_file,
		casa::String *calsky_file, casa::String *caltsys_file) {
	// argument analysis
	int opt = -1;
	int option_index = -1;
	struct option long_options[] = { { "help", no_argument, nullptr, 'h' }, {
			"input-file", required_argument, nullptr, 'i' }, { "calsky-file",
	required_argument, nullptr, 's' }, { "caltsys-file",
	required_argument, nullptr, 't' }, { "output-file",
	required_argument, nullptr, 'o' }, { 0, 0, 0, 0 } };

	while ((opt = getopt_long(argc, const_cast<char * const *>(argv),
			"hi:s:t:o:", long_options, &option_index)) != -1) {
		switch (opt) {
		case 'i':
			*input_file = optarg;
			break;
		case 's':
			*calsky_file = optarg;
			break;
		case 't':
			*caltsys_file = optarg;
			break;
		case 'o':
			*output_file = optarg;
			break;
		case 'h':
			PrintUsage(argv[0]);
			return;
		case '?':
			std::cout << "Invalid options were specified" << std::endl;
			PrintUsage(argv[0]);
			exit(1);
		}
	}

	if (output_file->length() == 0) {
		if (input_file->rfind("/") == input_file->length() - 1) {
			std::cout << "match" << std::endl;
		}
		*output_file = (
				(input_file->rfind("/") == input_file->length() - 1) ?
						input_file->substr(0, input_file->length() - 1) :
						*input_file) + "_out";
	}
	LOG4CXX_INFO(logger,
			"Input parameter summary:\n" << "\tinput-file=\"" << *input_file << "\"\n" << "\tcalsky-file=\"" << *calsky_file << "\"\n" << "\tcaltsys-file=\"" << *caltsys_file << "\"\n" << "\toutput-file=\"" << *output_file << "\"");
}

void E2eReduce(int argc, char const* const argv[]) {
	LOG4CXX_INFO(logger, "Enter: E2eReduce");
	bool const serialize = false;

	casa::String input_file;
	casa::String output_file;
	casa::String calsky_file;
	casa::String caltsys_file;
	double start_time = sakura_GetCurrentTime();
	AnalyseArguments(argc, argv, &input_file, &output_file, &calsky_file,
			&caltsys_file);
	if (input_file.length() > 0) {
		casa::Table table(input_file, casa::Table::Old);

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
