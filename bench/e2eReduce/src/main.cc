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
#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/ArrayIO.h>

#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/ExprNode.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/TableParse.h>

#include <libsakura/sakura.h>

#include "config_file_reader.h"
#include "option_parser.h"

namespace {
auto logger = log4cxx::Logger::getLogger("app");

struct Deleter {
	inline void operator()(void *ptr) const {
		free(ptr);
	}
};

class AligendArrayGenerator {
public:
	AligendArrayGenerator() :
			alignment_(sakura_GetAlignment()), index_(0), pointer_holder_(128) {
	}
	template<class T> inline T *GetAlignedArray(size_t num_elements) {
		size_t num_arena = (num_elements + alignment_ - 1) * sizeof(T);
		if (index_ >= pointer_holder_.size()) {
			pointer_holder_.resize(pointer_holder_.size() + 128);
		}
		void *ptr = malloc(num_arena);
		pointer_holder_[index_++].reset(ptr);
		return reinterpret_cast<T *>(sakura_AlignAny(num_arena, ptr,
				num_elements * sizeof(T)));
	}
private:
	size_t alignment_;
	size_t index_;
	std::vector<std::unique_ptr<void, Deleter> > pointer_holder_;
};

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
		uint8_t const f[], int edge_channels, bool do_clipping,
		float clipping_threshold, std::vector<uint64_t> const line_mask) {
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

std::string GetTaqlString(std::string table_name, unsigned int ifno) {
	std::ostringstream oss;
	oss << "SELECT FROM \"" << table_name << "\" WHERE IFNO == " << ifno
			<< " ORDER BY TIME";
	return oss.str();
}

casa::Table GetSelectedTable(std::string table_name, unsigned int ifno) {
	casa::String taql(GetTaqlString(table_name, ifno));
	return tableCommand(taql);
}

template<class T>
void GetArrayCell(T *array, size_t row_index,
		casa::ROArrayColumn<T> const column, casa::IPosition const cell_shape) {
	casa::Array < T > casa_array(cell_shape, array, casa::SHARE);
	column.get(row_index, casa_array);
}

template<class T>
void GetArrayColumn(T *array, casa::ROArrayColumn<T> const column,
		casa::IPosition const column_shape) {
	casa::Array < T > casa_array(column_shape, array, casa::SHARE);
	column.getColumn(casa_array);
}

template<class T>
void GetScalarColumn(T *array, casa::ROScalarColumn<T> const column,
		size_t num_rows) {
	casa::Vector < T
			> casa_array(casa::IPosition(1, num_rows), array, casa::SHARE);
	column.getColumn(casa_array);
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
	unsigned int ifno;
	OptionParser::ParseE2e(options, &input_file, &output_file, &ifno);

	std::string sky_table_name;
	std::string tsys_table_name;
	unsigned int tsys_ifno;
	OptionParser::ParseCalibration(options, &sky_table_name, &tsys_table_name,
			&tsys_ifno);

	int edge_channels;
	float clipping_threshold;
	bool do_clipping;
	OptionParser::ParseFlagging(options, &edge_channels, &clipping_threshold,
			&do_clipping);

	std::vector<uint64_t> line_mask;
	OptionParser::ParseBaseline(options, &line_mask);

	// config file summary
	{
		std::ostringstream oss;
		oss << "config file (" << configuration_file << ") summary:\n";
		oss << "\tinput filename=" << input_file << "\n\toutput filename="
				<< output_file << "\n\tspw=" << ifno << "\n";
		oss << "\tsky filename=" << sky_table_name << "\n\ttsys filename="
				<< tsys_table_name << "\n\ttsys_spw=" << tsys_ifno << "\n";
		oss << "\tedge channels=" << edge_channels << "\n\tclipping threshold=";
		if (do_clipping) {
			oss << clipping_threshold << "\n";
		} else {
			oss << "\n";
		}
		char separator = '[';
		oss << "\tline mask=";
		for (auto i = line_mask.begin(); i != line_mask.end(); ++i) {
			oss << separator << *i;
			separator = ',';
		}
		oss << "]";
		LOG4CXX_INFO(logger, oss.str());
	}

	double start_time = sakura_GetCurrentTime();
	if (input_file.size() > 0) {
		AligendArrayGenerator array_generator;

		casa::Table table = GetSelectedTable(input_file, ifno);

		casa::ROArrayColumn<float> spectra_column(table, "SPECTRA");
		casa::ROArrayColumn<unsigned char> flagtra_column(table, "FLAGTRA");
		unsigned int num_chan = spectra_column(0).nelements();

		// sky table
		casa::Table sky_table = GetSelectedTable(sky_table_name, ifno);
		size_t num_row_sky = sky_table.nrow();
		if (num_row_sky == 0) {
			throw "";
		}
		casa::ROArrayColumn<float> caldata_column(sky_table, "SPECTRA");
		casa::ROScalarColumn<double> time_column(sky_table, "TIME");
		float *sky_spectra = array_generator.GetAlignedArray<float>(
				num_chan * num_row_sky);
		double *sky_time = array_generator.GetAlignedArray<double>(
				num_row_sky);
		GetArrayColumn(sky_spectra, caldata_column,
				casa::IPosition(2, num_chan, num_row_sky));
		GetScalarColumn(sky_time, time_column, num_row_sky);

		// tsys table
		casa::Table tsys_table = GetSelectedTable(tsys_table_name, tsys_ifno);
		size_t num_row_tsys = tsys_table.nrow();
		if (num_row_tsys == 0) {
			throw "";
		}
		caldata_column.attach(tsys_table, "TSYS");
		time_column.attach(tsys_table, "TIME");
		size_t num_chan_tsys = caldata_column(0).nelements();
		float *tsys = array_generator.GetAlignedArray<float>(
				num_chan_tsys * num_row_tsys);
		double *tsys_time = array_generator.GetAlignedArray<double>(
				num_row_tsys);
		GetArrayColumn(tsys, caldata_column,
				casa::IPosition(2, num_chan_tsys, num_row_tsys));
		GetScalarColumn(tsys_time, time_column, num_row_tsys);

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
				ParallelJob(i, num_chan, v, f, edge_channels, do_clipping,
						clipping_threshold, line_mask);
			} else {
				group.async([=] {
					ParallelJob(i, num_chan, v, f, edge_channels, do_clipping,
							clipping_threshold, line_mask);
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
