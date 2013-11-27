/*
 * main.cc
 *
 *  Created on: 2013/10/21
 *      Author: kohji
 */

#include <cassert>
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
#include "context_handler.h"
#include "utils.h"

namespace {
auto logger = log4cxx::Logger::getLogger("app");

inline void ExecuteBitFlagToMask(size_t num_data, uint8_t const input_flag[],
		bool output_mask[]) {
//	std::cout << "ExecuteBitMaskToFlag" << std::endl;
	sakura_Uint8ToBool(num_data, input_flag, output_mask);
}

inline void ExecuteCalibration(size_t num_data, float const input_data[],
		CalibrationContext const &calibration_context, float out_data[]) {
//	std::cout << "ExecuteCalibration" << std::endl;

}

inline void ExecuteBaseline(sakura_BaselineContext const *context,
		size_t num_data, float const input_data[], float output_data[],
		bool mask[]) {
//	std::cout << "ExecuteBaseline" << std::endl;
}

inline void ExecuteSmoothing(sakura_Convolve1DContext const *context,
		size_t num_data, float const input_data[], float output_data[],
		bool mask[]) {
//	std::cout << "ExecuteSmoothing" << std::endl;
}

inline void ExecuteNanOrInfFlag(size_t num_data, float const input_data[],
		bool output_mask[]) {
//	std::cout << "ExecuteFlagNanOrInf" << std::endl;
	//sakura_SetFalseFloatIfNanOrInf(num_data, data, result);
}

inline void ExecuteFlagEdge(size_t num_edge, size_t num_data, bool mask[]) {
//	std::cout << "ExecuteFlagEdge" << std::endl;
}

inline void ExecuteClipping(float threshold, size_t num_data,
		float const input_data[], bool mask[]) {
//	std::cout << "ExecuteClipping" << std::endl;
}

inline void CalculateStatistics(size_t num_data, float const input_data[],
		bool const input_mask[]) {
//	std::cout << "CalculateStatistics" << std::endl;
}

void JobFinished(int i) {
	std::cout << "Job ";
	std::cout << i;
	std::cout << " have done.\n";
}

void ParallelJob(int job_id, size_t num_v, float const v[], uint8_t const f[],
		E2EOptions const &option_list,
		CalibrationContext const &calibration_context,
		sakura_BaselineContext const *baseline_context,
		sakura_Convolve1DContext const *convolve1d_context) {
	// generate temporary storage
	AlignedArrayGenerator generator;
	float *out_data = generator.GetAlignedArray<float>(num_v);
	bool *mask = generator.GetAlignedArray<bool>(num_v);

	// Execute Convert bit flag to boolean mask
	ExecuteBitFlagToMask(num_v, f, mask);

	// Execute Calibration
	ExecuteCalibration(num_v, v, calibration_context, out_data);

	// Execute Flag NaN/Inf
	ExecuteNanOrInfFlag(num_v, out_data, mask);

	// Execute Flag edge
	ExecuteFlagEdge(option_list.flagging.edge_channels, num_v, mask);

	// Execute Baseline
	ExecuteBaseline(baseline_context, num_v, out_data, out_data, mask);

	// Execute Clipping
	ExecuteClipping(option_list.flagging.clipping_threshold, num_v, out_data,
			mask);

	// Execute Smoothing
	ExecuteSmoothing(convolve1d_context, num_v, out_data, out_data, mask);

	// Calculate Statistics
	CalculateStatistics(num_v, out_data, mask);

//	float sum = 0.0;
//	for (unsigned int i = 0; i < num_v; ++i)
//		sum += v[i];
//	std::cout << "Job " << job_id << ": v[0]=" << v[0] << ", sum(v)=" << sum
//			<< std::endl;
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
		float *sky_spectra = nullptr;
		double *sky_time = nullptr;
		size_t num_chan_sky, num_row_sky;
		GetFromCalTable(options.calibration.sky_table, options.ifno, "SPECTRA",
				&array_generator, &sky_spectra, &sky_time, &num_chan_sky,
				&num_row_sky);
		if (num_chan != num_chan_sky) {
			throw "";
		}
		assert(sky_spectra != nullptr && sky_time != nullptr);

		// tsys table
		float *tsys = nullptr;
		double *tsys_time = nullptr;
		size_t num_chan_tsys = 0, num_row_tsys = 0;
		GetFromCalTable(options.calibration.tsys_table,
				options.calibration.tsys_ifno, "TSYS", &array_generator, &tsys,
				&tsys_time, &num_chan_tsys, &num_row_tsys);
		assert(tsys != nullptr && tsys_time != nullptr);

		// get frequency label from the table
		// for spectral data
		double *frequency_label_target =
				array_generator.GetAlignedArray<double>(num_chan);
		GetFrequencyLabelFromScantable(table, options.ifno, num_chan,
				frequency_label_target);

		// for Tsys
		double *frequency_label_tsys = array_generator.GetAlignedArray<double>(
				num_chan_tsys);
		GetFrequencyLabelFromScantable(options.calibration.tsys_table,
				options.calibration.tsys_ifno, num_chan_tsys,
				frequency_label_tsys);

		// Create Context and struct for calibration
		// calibration context
		CalibrationContext calibration_context;
		calibration_context.num_channel_sky = num_chan;
		calibration_context.num_channel_tsys = num_chan_tsys;
		calibration_context.num_data_sky = num_row_sky;
		calibration_context.num_data_tsys = num_row_tsys;
		calibration_context.timestamp_sky = sky_time;
		calibration_context.timestamp_tsys = tsys_time;
		calibration_context.sky_spectra = sky_spectra;
		calibration_context.tsys = tsys;

		// baseline context
		sakura_BaselineContext *baseline_context = nullptr;
		if (sakura_CreateBaselineContext(options.baseline.baseline_type,
				options.baseline.order, num_chan, &baseline_context)
				!= sakura_Status_kOK) {
			throw "";
		}
		assert(baseline_context != nullptr);
		std::unique_ptr<sakura_BaselineContext, BaselineContextDeleter> baseline_context_holder(
				baseline_context);

		// convolve1d context
		sakura_Convolve1DContext *convolve1d_context = nullptr;
		if (sakura_CreateConvolve1DContext(num_chan,
				options.smoothing.kernel_type, options.smoothing.kernel_width,
				true, &convolve1d_context) != sakura_Status_kOK) {
			throw "";
		}
		assert(convolve1d_context != nullptr);
		std::unique_ptr<sakura_Convolve1DContext, Convolve1DContextDeleter> convolve1d_context_holder(
				convolve1d_context);

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
				ParallelJob(i, num_chan, v, f, options, calibration_context,
						baseline_context, convolve1d_context);
			} else {
				group.async([=] {
					ParallelJob(i, num_chan, v, f, options, calibration_context,
							baseline_context, convolve1d_context);
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
