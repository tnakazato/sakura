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
struct StaticInitializer {
	StaticInitializer() {
		::log4cxx::PropertyConfigurator::configure("config.log4j");
	}
};
static StaticInitializer initializer;

auto logger = log4cxx::Logger::getLogger("app");

inline void ExecuteBitFlagToMask(size_t num_data, uint8_t const input_flag[],
bool output_mask[]) {
	sakura_Status status = sakura_Uint8ToBool(num_data, input_flag,
			output_mask);
	if (status != sakura_Status_kOK) {
		throw std::runtime_error("bitflag_to_bool");
	}
	status = sakura_InvertBool(num_data, output_mask, output_mask);
	if (status != sakura_Status_kOK) {
		throw std::runtime_error("bitflag_to_bool");
	}
}

inline void ExecuteCalibration(double timestamp, size_t num_data,
		float const input_data[], CalibrationContext const &calibration_context,
		float out_data[]) {
	// interpolate Tsys to timestamp
	size_t num_channel = calibration_context.num_channel_sky;
	float *tsys = calibration_context.tsys_work;
	double *t = calibration_context.timestamp_work;
	t[0] = timestamp;
	sakura_Status status = sakura_InterpolateYAxisFloat(
			sakura_InterpolationMethod_kLinear, 0,
			calibration_context.num_data_tsys,
			calibration_context.timestamp_tsys, num_channel,
			calibration_context.tsys, 1, t, tsys);

	if (status != sakura_Status_kOK) {
		throw std::runtime_error("calibration");
	}

	// interpolated sky spectrum to timestamp
	auto *sky = calibration_context.sky_work;
	status = sakura_InterpolateYAxisFloat(sakura_InterpolationMethod_kLinear, 0,
			calibration_context.num_data_sky, calibration_context.timestamp_sky,
			num_channel, calibration_context.sky_spectra, 1, t, sky);

	if (status != sakura_Status_kOK) {
		throw std::runtime_error("calibration");
	}

	// apply calibration
	// TODO: allow reference == result?
	status = sakura_ApplyPositionSwitchCalibration(num_channel, tsys,
			num_channel, input_data, sky, out_data);

	if (status != sakura_Status_kOK) {
		throw std::runtime_error("calibration");
	}
}

inline void ExecuteBaseline(float clipping_threshold_sigma,
		uint16_t num_fitting_max, bool mask_after_clipping[],
		sakura_BaselineContext const *context, size_t num_data,
		float const input_data[], float output_data[], bool mask[]) {
	sakura_Status status = sakura_SubtractBaseline(num_data, input_data, mask,
			context, clipping_threshold_sigma, num_fitting_max, true,
			mask_after_clipping, output_data);
	if (status != sakura_Status_kOK) {
		throw std::runtime_error("baseline");
	}
}

inline void ExecuteSmoothing(sakura_Convolve1DContext const *context,
		size_t num_data, float const input_data_arg[], float output_data_arg[],
		bool mask_arg[]) {
	auto input_data = AssumeAligned(input_data_arg);
	auto output_data = AssumeAligned(output_data_arg);
	auto mask = AssumeAligned(mask_arg);

	// apply mask by filling zero like CASA
	for (size_t i = 0; i < num_data; ++i) {
		if (mask[i]) {
			output_data[i] = input_data[i];
		} else {
			output_data[i] = 0.0f;
		}
	}
	sakura_Status status = sakura_Convolve1D(context, num_data, output_data,
			mask, output_data);
	if (status != sakura_Status_kOK) {
		throw std::runtime_error("convolution");
	}
}

inline void AlignedBoolAnd(size_t elements, bool const src[], bool dst[]) {
	auto dst_u8 = AssumeAligned(reinterpret_cast<uint8_t *>(dst));
	auto src_u8 = AssumeAligned(reinterpret_cast<uint8_t const*>(src));
	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);
	STATIC_ASSERT(sizeof(*dst_u8) == sizeof(*dst));
	for (size_t i = 0; i < elements; ++i) {
		dst_u8[i] &= src_u8[i];
	}
}
inline void ExecuteNanOrInfFlag(size_t num_data, float const input_data[],
bool mask_arg[]) {
	auto mask = AssumeAligned(mask_arg);

//	std::cout << "ExecuteFlagNanOrInf" << std::endl;
	AlignedArrayGenerator generator;
	bool *mask_local = generator.GetAlignedArray<bool>(num_data); // TODO move mask_local to struct WorkingSet
	sakura_Status status = sakura_SetFalseFloatIfNanOrInf(num_data, input_data,
			mask_local);
	if (status != sakura_Status_kOK) {
		throw std::runtime_error("NaN and Inf flag");
	}
	// Merging mask
	AlignedBoolAnd(num_data, mask_local, mask);
}

inline void ExecuteFlagEdge(size_t num_edge, size_t num_data, bool mask[]) {
//	std::cout << "ExecuteFlagEdge" << std::endl;
	AlignedArrayGenerator generator;
	alignas(32)
	int lower[] = { static_cast<int>(num_edge - 1) }; // it's ok to be larger than num_data
	alignas(32)
	int upper[] = { static_cast<int>(num_data - num_edge) }; // it's ok to be negative
	STATIC_ASSERT(ELEMENTSOF(lower) == ELEMENTSOF(upper));
	size_t num_condition = ELEMENTSOF(lower);
	if (lower[0] > upper[0]) {
		// All channels should be False (False)
		num_condition = 0;
		lower[0] = -1;
		upper[0] = -1;
	}
	bool *mask_local = generator.GetAlignedArray<bool>(num_data);
	int *channel_id = generator.GetAlignedArray<int>(num_data); // TODO make channel_id static and shared among threads, and initialize once
	for (size_t i = 0; i < num_data; ++i) {
		channel_id[i] = static_cast<int>(i);
	}
	sakura_Status status = sakura_SetTrueIntInRangesExclusive(num_data,
			channel_id, num_condition, lower, upper, mask_local);
	if (status != sakura_Status_kOK) {
		throw std::runtime_error("Flag Edge");
	}
	// Merging mask
	AlignedBoolAnd(num_data, mask_local, mask);
}

inline void ExecuteClipping(float threshold, size_t num_data,
		float const input_data[], bool mask[]) {
//	std::cout << "ExecuteClipping" << std::endl;
	AlignedArrayGenerator generator;
	bool *mask_local = generator.GetAlignedArray<bool>(num_data);
	alignas(32)
	float lower[] = { -std::abs(threshold) };
	alignas(32)
	float upper[] = { std::abs(threshold) };
	STATIC_ASSERT(ELEMENTSOF(lower) == ELEMENTSOF(upper));
	// simple [-threshold, threshold] clipping
	sakura_Status status = sakura_SetTrueFloatInRangesExclusive(num_data,
			input_data, ELEMENTSOF(lower), lower, upper, mask_local);
	if (status != sakura_Status_kOK) {
		throw std::runtime_error("Clip Flag");
	}
	// Merging mask
	AlignedBoolAnd(num_data, mask_local, mask);
}

inline void CalculateStatistics(size_t num_data, float const input_data[],
bool const input_mask[]) {
	sakura_StatisticsResult result;
	sakura_Status status = sakura_ComputeStatistics(num_data, input_data,
			input_mask, &result);
	if (status != sakura_Status_kOK) {
		throw std::runtime_error("statistics");
	}
}

void JobFinished(size_t i) {
	LOG4CXX_DEBUG(logger, "Row " << i << " have done.");
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
	casa::ROScalarColumn<double> time_column;

	~SharedWorkingSet() {
		LOG4CXX_DEBUG(logger, "Destructing SharedWorkingSet");
		if (baseline_context != nullptr) {
			sakura_DestroyBaselineContext(baseline_context);
		}
		if (convolve1d_context != nullptr) {
			sakura_DestroyConvolve1DContext(convolve1d_context);
		}
	}
};

size_t const rows_per_processing = 8;

void ParallelJob(size_t job_id, size_t jobs, size_t num_v,
		double const tarray[/*jobs*/],
		float const * const varray[rows_per_processing],
		uint8_t const * const farray[rows_per_processing],
		E2EOptions const &option_list, SharedWorkingSet const *shared) {
// generate temporary storage
	AlignedArrayGenerator generator;
	float *out_data[rows_per_processing];
	for (size_t i = 0; i < jobs; ++i) {
		out_data[i] = generator.GetAlignedArray<float>(num_v);
	}
	bool *mask = generator.GetAlignedArray<bool>(num_v);
	bool *baseline_after_mask = generator.GetAlignedArray<bool>(num_v);

	for (size_t i = 0; i < jobs; ++i) {
		auto const v = varray[i];
		auto const f = farray[i];
		double const t = tarray[i];
		// Convert bit flag to boolean mask
		ExecuteBitFlagToMask(num_v, f, mask);

		// Execute Calibration
		ExecuteCalibration(t, num_v, v, shared->calibration_context,
				out_data[i]);

		// Execute Flag NaN/Inf
		ExecuteNanOrInfFlag(num_v, out_data[i], mask);

		// Execute Flag edge
		ExecuteFlagEdge(option_list.flagging.edge_channels, num_v, mask);

		// Execute Baseline
		ExecuteBaseline(option_list.baseline.clipping_threshold,
				option_list.baseline.num_fitting_max, baseline_after_mask,
				shared->baseline_context, num_v, out_data[i], out_data[i],
				mask);

		// Execute Clipping
		ExecuteClipping(option_list.flagging.clipping_threshold, num_v,
				out_data[i], mask);

		// Execute Smoothing
		ExecuteSmoothing(shared->convolve1d_context, num_v, out_data[i],
				out_data[i], mask);

		// Calculate Statistics
		CalculateStatistics(num_v, out_data[i], mask);

	}
// sleep(job_id % 3); // my job is sleeping.
	xdispatch::main_queue().sync([=]() {
		// run casa operations in main_queue.
		// write out_data[job_id + i] where i = 0 .. jobs-1
			JobFinished(job_id);
		});
}

SharedWorkingSet *InitializeSharedWorkingSet(E2EOptions const &options,
		unsigned int polno, AlignedArrayGenerator *array_generator) {
	LOG4CXX_DEBUG(logger, "Enter: InitializeSharedWorkingSet");
	assert(array_generator != nullptr);

// TODO: need loop on POLNO
	SharedWorkingSet *shared = new SharedWorkingSet { 0, 0, { }, nullptr,
			nullptr, GetSelectedTable(options.input_file, options.ifno, polno) };
	try {
		shared->num_data = shared->table.nrow();
		if (shared->num_data > 0) {
			shared->spectra_column.attach(shared->table, "SPECTRA");
			shared->flagtra_column.attach(shared->table, "FLAGTRA");
			shared->time_column.attach(shared->table, "TIME");
			shared->num_chan = shared->spectra_column(0).nelements();
			// Create Context and struct for calibration
			// calibration context
			FillCalibrationContext(options.calibration.sky_table,
					options.calibration.tsys_table, options.ifno,
					options.calibration.tsys_ifno, polno, array_generator,
					&shared->calibration_context);
			if (shared->calibration_context.num_channel_sky
					!= shared->num_chan) {
				throw std::runtime_error("calibration");
			}
			// baseline context
			if (sakura_CreateBaselineContext(options.baseline.baseline_type,
					options.baseline.order, shared->num_chan,
					&shared->baseline_context) != sakura_Status_kOK) {
				throw std::runtime_error("baseline");
			}
			assert(shared->baseline_context != nullptr);
			// convolve1d context
			if (sakura_CreateConvolve1DContext(shared->num_chan,
					options.smoothing.kernel_type,
					options.smoothing.kernel_width, options.smoothing.use_fft,
					&shared->convolve1d_context) != sakura_Status_kOK) {
				throw std::runtime_error("convolution");
			}
			assert(shared->convolve1d_context != nullptr);
		}
	} catch (...) {
		delete shared;
		LOG4CXX_DEBUG(logger, "Leave: InitializeSharedWorkingSet");
		throw;
	}
	LOG4CXX_DEBUG(logger, "Leave: InitializeSharedWorkingSet");
	return shared;
}

struct WorkingSet {
	SharedWorkingSet *shared;
	float *spectrum[rows_per_processing];
	unsigned char *flag[rows_per_processing];
	double timestamp[rows_per_processing];
};

WorkingSet *InitializeWorkingSet(SharedWorkingSet *shared, size_t elements,
		AlignedArrayGenerator *array_generator) {
	std::unique_ptr<WorkingSet[]> result(new WorkingSet[elements]);
	for (size_t i = 0; i < elements; ++i) {
		result[i].shared = shared;
		for (size_t j = 0; j < rows_per_processing; ++j) {
			result[i].spectrum[j] = array_generator->GetAlignedArray<float>(
					shared->num_chan);
			result[i].flag[j] = array_generator->GetAlignedArray<unsigned char>(
					shared->num_chan);
		}
	}
	return result.release();
}

void Reduce(E2EOptions const &options) {
	LOG4CXX_DEBUG(logger, "Enter: Reduce");
	unsigned int const num_polarizations = 4;
	for (unsigned int polno = 0; polno < num_polarizations; polno++) {
		LOG4CXX_INFO(logger, "Processing POLNO " << polno);
		AlignedArrayGenerator array_generator;
		SharedWorkingSet *shared = nullptr;
		xdispatch::main_queue().sync(
				[&shared, &options, &polno, &array_generator] {
					// run casa operations in main_queue.
					shared = InitializeSharedWorkingSet(options, polno, &array_generator);
				});
		ScopeGuard guard_for_shared([shared] {
			delete shared;
		});

		if (shared->num_data == 0) {
			continue;
		}

		size_t num_threads = std::min(shared->num_data,
				size_t(options.max_threads));
		std::unique_ptr<WorkingSet[]> working_sets(
				InitializeWorkingSet(shared, num_threads, &array_generator));

		std::vector<size_t> available_workers(num_threads);
		for (size_t i = 0; i < num_threads; ++i) {
			available_workers[i] = i;
		}

		LOG4CXX_DEBUG(logger, "Initialized");

		bool const serialize = options.serialize;
		xdispatch::queue serial_queue("my serial queue");
		xdispatch::semaphore semaphore(num_threads);

		auto group = xdispatch::group();
		for (size_t i = 0; i < shared->num_data; i += rows_per_processing) {
			size_t rows = std::min(rows_per_processing, shared->num_data - i);
			semaphore.acquire();
			size_t working_set_id = 0;
			serial_queue.sync([&working_set_id, &available_workers] {
				assert(available_workers.size() > 0);
				working_set_id = available_workers.back();
				available_workers.pop_back();
			});
			LOG4CXX_DEBUG(logger, "working_set_id: " << working_set_id);
			// get data from the table
			xdispatch::main_queue().sync([=, &working_sets] {
				// run casa operations in main_queue.
					LOG4CXX_DEBUG(logger, "Reading record " << i);
					for (size_t j = 0; j < rows; ++j) {
						GetArrayCell(working_sets[working_set_id].spectrum[j], i+j, shared->spectra_column,
								casa::IPosition(1, shared->num_chan));
					}
					for (size_t j = 0; j < rows; ++j) {
						GetArrayCell(working_sets[working_set_id].flag[j], i+j, shared->flagtra_column,
								casa::IPosition(1, shared->num_chan));
					}
					GetScalarCells(working_sets[working_set_id].timestamp, i, rows, shared->time_column);
				});

			auto v =
					reinterpret_cast<float const * const *>(working_sets[working_set_id].spectrum);
			auto f =
					reinterpret_cast<uint8_t const * const *>(&working_sets[working_set_id].flag);
			auto t =
					reinterpret_cast<double const *>(working_sets[working_set_id].timestamp);
			if (serialize) {
				ScopeGuard cleanup(
						[=,&available_workers, &serial_queue, &semaphore] {
							serial_queue.sync([=,&working_set_id, &available_workers] {
										available_workers.push_back(working_set_id);
									});
							semaphore.release();
						});
				ParallelJob(i, rows, shared->num_chan, t, v, f, options,
						shared);
			} else {
				group.async(
						[=, &available_workers, &serial_queue, &semaphore, &options] {
							ScopeGuard cleanup([=,&available_workers, &serial_queue, &semaphore] {
										serial_queue.sync([=,&working_set_id, &available_workers] {
													available_workers.push_back(working_set_id);
												});
										semaphore.release();
									});
							ParallelJob(i, rows, shared->num_chan, t, v, f, options, shared);
						});
			}
		}
		group.wait(xdispatch::time_forever);
	}
	LOG4CXX_DEBUG(logger, "Leave: Reduce");
}

void E2eReduce(int argc, char const* const argv[]) {
	LOG4CXX_DEBUG(logger, "Enter: E2eReduce");

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
	LOG4CXX_DEBUG(logger, "Leave: E2eReduce");
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
	xdispatch::global_queue().async([=] {
		main_(argc, argv);
	});
	xdispatch::exec(); // never returns
	return 0;
}
