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
#include <stdexcept>
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
	sakura_StatisticsResult result;
	sakura_Status status = sakura_ComputeStatistics(num_data, input_data, input_mask, &result);
	if (status != sakura_Status_kOK) {
		throw std::runtime_error("");
	}
}

void JobFinished(int i) {
	std::cout << "Job ";
	std::cout << i;
	std::cout << " have done.\n";
}

struct SharedWorkingSet {
	size_t num_chan;
	size_t num_data;
	CalibrationContext calibration_context;
	sakura_BaselineContext *baseline_context;
	sakura_Convolve1DContext *convolve1d_context;

	casa::Table table;
	casa::ROArrayColumn<float> spectra_column;
	casa::ROArrayColumn<unsigned char> flagtra_column;

	void Release() {
		LOG4CXX_INFO(logger, "Releasing SharedWorkingSet");
		if (baseline_context != nullptr) {
			sakura_DestroyBaselineContext(baseline_context);
			baseline_context = nullptr;
		}
		if (convolve1d_context != nullptr) {
			sakura_DestroyConvolve1DContext(convolve1d_context);
			convolve1d_context = nullptr;
		}
	}
};

void ParallelJob(int job_id, size_t num_v, float const v[], uint8_t const f[],
		E2EOptions const &option_list, SharedWorkingSet const *shared) {
	// generate temporary storage
	AlignedArrayGenerator generator;
	float *out_data = generator.GetAlignedArray<float>(num_v);
	bool *mask = generator.GetAlignedArray<bool>(num_v);

	// Convert bit flag to boolean mask
	ExecuteBitFlagToMask(num_v, f, mask);

	// Execute Calibration
	ExecuteCalibration(num_v, v, shared->calibration_context, out_data);

	// Execute Flag NaN/Inf
	ExecuteNanOrInfFlag(num_v, out_data, mask);

	// Execute Flag edge
	ExecuteFlagEdge(option_list.flagging.edge_channels, num_v, mask);

	// Execute Baseline
	ExecuteBaseline(shared->baseline_context, num_v, out_data, out_data, mask);

	// Execute Clipping
	ExecuteClipping(option_list.flagging.clipping_threshold, num_v, out_data,
			mask);

	// Execute Smoothing
	ExecuteSmoothing(shared->convolve1d_context, num_v, out_data, out_data,
			mask);

	// Calculate Statistics
	CalculateStatistics(num_v, out_data, mask);

//	float sum = 0.0;
//	for (unsigned int i = 0; i < num_v; ++i)
//		sum += v[i];
//	std::cout << "Job " << job_id << ": v[0]=" << v[0] << ", sum(v)=" << sum
//			<< std::endl;
	// sleep(job_id % 3); // my job is sleeping.
	xdispatch::main_queue().sync([=]() {
		// run casa operations in main_queue.
		JobFinished(job_id);
	});
}

SharedWorkingSet *InitializeSharedWorkingSet(E2EOptions const &options,
		AlignedArrayGenerator *array_generator) {
	LOG4CXX_INFO(logger, "Enter: InitializeSharedWorkingSet");
	assert(array_generator != nullptr);

	SharedWorkingSet *shared = new SharedWorkingSet { 0, 0, { }, nullptr, nullptr,
			GetSelectedTable(options.input_file, options.ifno) };
	shared->spectra_column.attach(shared->table, "SPECTRA");
	shared->flagtra_column.attach(shared->table, "FLAGTRA");
	shared->num_chan = shared->spectra_column(0).nelements();
	shared->num_data = shared->table.nrow();
	try {
		// Create Context and struct for calibration
		// calibration context
		FillCalibrationContext(options.calibration.sky_table,
				options.calibration.tsys_table, options.ifno,
				options.calibration.tsys_ifno, array_generator,
				&shared->calibration_context);
		if (shared->calibration_context.num_channel_sky != shared->num_chan) {
			throw "";
		}
		// baseline context
		if (sakura_CreateBaselineContext(options.baseline.baseline_type,
				options.baseline.order, shared->num_chan,
				&shared->baseline_context) != sakura_Status_kOK) {
			throw "";
		}
		assert(shared->baseline_context != nullptr);
		// convolve1d context
		if (sakura_CreateConvolve1DContext(shared->num_chan,
				options.smoothing.kernel_type, options.smoothing.kernel_width,
				options.smoothing.use_fft, &shared->convolve1d_context)
				!= sakura_Status_kOK) {
			throw "";
		}
		assert(shared->convolve1d_context != nullptr);
	} catch (...) {
		ScopeGuard guard_for_shared([shared] {
			delete shared;
		});
		shared->Release();
		throw;
	}
	LOG4CXX_INFO(logger, "Leave: InitializeSharedWorkingSet");
	return shared;
}

struct WorkingSet {
	SharedWorkingSet *shared;
	float *spectrum;
	unsigned char *flag;
};

WorkingSet *InitializeWorkingSet(SharedWorkingSet *shared, size_t elements,
		AlignedArrayGenerator *array_generator) {
	std::unique_ptr<WorkingSet[]> result(new WorkingSet[elements]);
	for (size_t i = 0; i < elements; ++i) {
		result[i].shared = shared;
		result[i].spectrum = array_generator->GetAlignedArray<float>(
				shared->num_chan);
		result[i].flag = array_generator->GetAlignedArray<unsigned char>(
				shared->num_chan);
	}
	return result.release();
}

void Reduce(E2EOptions const &options) {
	LOG4CXX_INFO(logger, "Enter: Reduce");
	AlignedArrayGenerator array_generator;
	SharedWorkingSet *shared = nullptr;
	xdispatch::main_queue().sync([&shared, &options, &array_generator] {
		// run casa operations in main_queue.
			shared = InitializeSharedWorkingSet(options, &array_generator);
		});
	ScopeGuard guard_for_shared([shared] {
		ScopeGuard guard_for_shared([shared] {
					delete shared;
				});
		shared->Release();
	});

	size_t num_threads = std::min(shared->num_data, size_t(20));
	std::unique_ptr<WorkingSet[]> working_sets(
			InitializeWorkingSet(shared, num_threads, &array_generator));

	std::vector<size_t> available_workers(num_threads);
	for (size_t i = 0; i < num_threads; ++i) {
		available_workers[i] = i;
	}

	LOG4CXX_INFO(logger, "Initialized");

	bool const serialize = false;
	xdispatch::queue serial_queue("serial");
	xdispatch::semaphore semaphore(num_threads);
	size_t rows = shared->spectra_column.nrow();

	auto group = xdispatch::group();
	for (size_t i = 0; i < rows; ++i) {
		semaphore.acquire();
		size_t working_set_id = 0;
		serial_queue.sync([&working_set_id, &available_workers] {
			working_set_id = available_workers.back();
			available_workers.pop_back();
		});
		LOG4CXX_INFO(logger, "working_set_id: " << working_set_id);
		// get data from the table
		xdispatch::main_queue().sync([=, &working_sets] {
			// run casa operations in main_queue.
			LOG4CXX_INFO(logger, "Reading record " << i);
				GetArrayCell(working_sets[working_set_id].spectrum, i, shared->spectra_column,
						casa::IPosition(1, shared->num_chan));
				GetArrayCell(working_sets[working_set_id].flag, i, shared->flagtra_column,
						casa::IPosition(1, shared->num_chan));
			});

		const float *v = working_sets[working_set_id].spectrum;
		const uint8_t *f =
				reinterpret_cast<const uint8_t*>(working_sets[working_set_id].flag);
		if (serialize) {
			ParallelJob(i, shared->num_chan, v, f, options, shared);
			serial_queue.sync([&working_set_id, &available_workers] {
				available_workers.push_back(working_set_id);
			});
			semaphore.release();
		} else {
			group.async(
					[=, &available_workers, &serial_queue, &semaphore, &options] {
						ScopeGuard cleanup([=,&available_workers, &serial_queue, &semaphore] {
									serial_queue.sync([=,&working_set_id, &available_workers] {
												available_workers.push_back(working_set_id);
											});
									semaphore.release();
								});
						ParallelJob(i, shared->num_chan, v, f, options, shared);
					});
		}
	}
	group.wait(xdispatch::time_forever);
	LOG4CXX_INFO(logger, "Leave: Reduce");
}

void E2eReduce(int argc, char const* const argv[]) {
	LOG4CXX_INFO(logger, "Enter: E2eReduce");

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
		Reduce(options);
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
