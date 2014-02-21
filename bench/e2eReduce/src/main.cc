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

namespace {
struct StaticInitializer {
	StaticInitializer() {
		::log4cxx::PropertyConfigurator::configure("config.log4j");
	}
};
static StaticInitializer initializer;

auto logger = log4cxx::Logger::getLogger("app");
}

#include "option_parser.h"
#include "context_handler.h"
#include "utils.h"

namespace {
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

inline void ExecuteMaskToBitFlag(size_t num_data, bool input_mask[],
		uint8_t output_flag_arg[]) {
	auto output_flag = AssumeAligned(output_flag_arg);
	sakura_Status status = sakura_InvertBool(num_data, input_mask, input_mask);
	// Now input_mask[i]=T if i-th channel should be flagged.
	if (status != sakura_Status_kOK) {
		throw std::runtime_error("sakura_InvertBool");
	}
	for (size_t i = 0; i < num_data; ++i) {
		output_flag[i] = 0;
	}
	// Now all elements in output_flag is false.
	uint8_t bit_mask = 1 << 7;
	status = sakura_OperateBitsUint8Or(bit_mask, num_data, output_flag,
			input_mask, output_flag);
	if (status != sakura_Status_kOK) {
		throw std::runtime_error("sakura_OperateBitsUint8Or");
	}
}

inline void ExecuteCalibration(size_t num_row, double const timestamp[], size_t num_data,
		float const input_data[], CalibrationContext const &calibration_context,
		float sky[], float out_data[], float out_tsys[]) {
	assert(num_data == calibration_context.num_channel);

	// interpolate Tsys to timestamp
	sakura_Status status = sakura_InterpolateYAxisFloat(
			sakura_InterpolationMethod_kLinear, 0,
			calibration_context.num_data_tsys,
			calibration_context.timestamp_tsys, num_data,
			calibration_context.tsys, num_row, timestamp, out_tsys);

	if (status != sakura_Status_kOK) {
		throw std::runtime_error("calibration");
	}

	// interpolated sky spectrum to timestamp
	status = sakura_InterpolateYAxisFloat(sakura_InterpolationMethod_kLinear, 0,
			calibration_context.num_data_sky, calibration_context.timestamp_sky,
			num_data, calibration_context.sky_spectra, num_row, timestamp, sky);

	if (status != sakura_Status_kOK) {
		throw std::runtime_error("calibration");
	}

	// apply calibration
	// TODO: allow reference == result?
	status = sakura_ApplyPositionSwitchCalibration(num_data * num_row, out_tsys,
			num_data * num_row, input_data, sky, out_data);

	if (status != sakura_Status_kOK) {
		throw std::runtime_error("calibration");
	}
}

inline void ExecuteBaseline(float clipping_threshold_sigma,
		uint16_t num_fitting_max, bool const mask[],
		sakura_BaselineContext const *context, size_t num_data,
		float const input_data[], float output_data[],
		bool mask_after_clipping[]) {
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

struct SharedWorkingSet {
	size_t num_chan;
	size_t num_data;
	CalibrationContext calibration_context;
	sakura_BaselineContext *baseline_context;
	sakura_Convolve1DContext *convolve1d_context;
	bool *line_mask;

	// tables
	casa::Table input_table;
	casa::Table output_table;

	// input columns
	casa::ROArrayColumn<float> input_spectra_column;
	casa::ROArrayColumn<unsigned char> input_flagtra_column;
	casa::ROScalarColumn<double> input_time_column;

	// output columns
	casa::ArrayColumn<float> output_spectra_column;
	casa::ArrayColumn<unsigned char> output_flagtra_column;
	casa::ArrayColumn<float> output_tsys_column;

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

void PutResult(SharedWorkingSet *shared, size_t row, size_t num_channels,
		float out_data[], unsigned char out_flag[], float out_tsys[]) {
	LOG4CXX_DEBUG(logger, "Writing result to row " << row);
	size_t row_number = shared->input_table.rowNumbers()[row];
	LOG4CXX_DEBUG(logger, "Row number of original table is " << row_number);
	shared->output_spectra_column.put(row_number,
			casa::Array<float>(casa::IPosition(1, num_channels), out_data,
					casa::SHARE));
	shared->output_flagtra_column.put(row_number,
			casa::Array<unsigned char>(casa::IPosition(1, num_channels),
					out_flag, casa::SHARE));
	shared->output_tsys_column.put(row_number,
			casa::Array<float>(casa::IPosition(1, num_channels), out_tsys,
					casa::SHARE));
}

size_t const rows_per_processing = 64;
struct WorkArea {
	bool *mask;
	bool *baseline_mask;
	bool *baseline_after_mask;
	double *timestamp;
	float *sky;
};

void ParallelJob(size_t row_id, size_t rows, size_t num_v,
		double const tarray[/*jobs*/],
		float const * const varray[rows_per_processing],
		uint8_t const * const farray[rows_per_processing],
		E2EOptions const &option_list, SharedWorkingSet const *shared,
		WorkArea *work_area, float **out_data, uint8_t **out_flag,
		float **out_tsys) {

	float *sky = AssumeAligned(work_area->sky);
	double *timestamp = AssumeAligned(work_area->timestamp);

	bool *mask = AssumeAligned(work_area->mask);
	bool *baseline_mask = AssumeAligned(work_area->baseline_mask);
	bool *baseline_after_mask = AssumeAligned(work_area->baseline_after_mask);

	for (size_t i = 0; i < rows; ++i) {
		auto const v = AssumeAligned(varray[i]);
		auto const f = AssumeAligned(farray[i]);
		double const t = tarray[i];
		// Convert bit flag to boolean mask
		ExecuteBitFlagToMask(num_v, f, mask);

		// Execute Calibration
		timestamp[0] = t;
		ExecuteCalibration(1, timestamp, num_v, v, shared->calibration_context,
				sky, out_data[i], out_tsys[i]);

		// Execute Flag NaN/Inf
		ExecuteNanOrInfFlag(num_v, out_data[i], mask);

		// Execute Flag edge
		ExecuteFlagEdge(option_list.flagging.edge_channels, num_v, mask);

		// Execute Baseline
		for (size_t ichan = 0; ichan < num_v; ++ichan) {
			baseline_mask[ichan] = mask[ichan];
		}
		AlignedBoolAnd(num_v, shared->line_mask, baseline_mask);
		ExecuteBaseline(option_list.baseline.clipping_threshold,
				option_list.baseline.num_fitting_max, baseline_mask,
				shared->baseline_context, num_v, out_data[i], out_data[i],
				baseline_after_mask);

		// Execute Clipping
		ExecuteClipping(option_list.flagging.clipping_threshold, num_v,
				out_data[i], mask);

		// Execute Smoothing
		ExecuteSmoothing(shared->convolve1d_context, num_v, out_data[i],
				out_data[i], mask);

		// Calculate Statistics
		CalculateStatistics(num_v, out_data[i], mask);

		// Convert boolean mask back to bit flag
		uint8_t *uint8_flag = reinterpret_cast<uint8_t *>(out_flag[i]);
		ExecuteMaskToBitFlag(num_v, mask, uint8_flag);

	}
}

SharedWorkingSet *InitializeSharedWorkingSet(E2EOptions const &options,
		unsigned int polno, AlignedArrayGenerator *array_generator) {
	LOG4CXX_DEBUG(logger, "Enter: InitializeSharedWorkingSet");
	assert(array_generator != nullptr);

	SharedWorkingSet *shared =
			new SharedWorkingSet { 0, 0, { }, nullptr, nullptr, nullptr,
					GetSelectedTable(options.input_file, options.ifno, polno,
					true), casa::Table(options.output_file, casa::Table::Update) };
	try {
		shared->num_data = shared->input_table.nrow();
		if (shared->num_data > 0) {
			// attach input columns
			shared->input_spectra_column.attach(shared->input_table, "SPECTRA");
			shared->input_flagtra_column.attach(shared->input_table, "FLAGTRA");
			shared->input_time_column.attach(shared->input_table, "TIME");
			shared->num_chan = shared->input_spectra_column(0).nelements();
			// attach output columns
			shared->output_spectra_column.attach(shared->output_table,
					"SPECTRA");
			shared->output_flagtra_column.attach(shared->output_table,
					"FLAGTRA");
			shared->output_tsys_column.attach(shared->output_table, "TSYS");
			// line mask for baseline fitting
			bool *mask = array_generator->GetAlignedArray<bool>(
					shared->num_chan);
			for (size_t i = 0; i < shared->num_chan; ++i) {
				mask[i] = false;
			}
			for (size_t i = 0; i < options.baseline.line_mask.size(); i += 2) {
				size_t start = options.baseline.line_mask[i];
				size_t end = std::min(options.baseline.line_mask[i + 1],
						shared->num_chan);
				for (size_t j = start; j <= end; ++j) {
					mask[j] = true;
				}
			}
			shared->line_mask = mask;
			// Create Context and struct for calibration
			// calibration context
			FillCalibrationContext(options.calibration.sky_table,
					options.calibration.tsys_table, options.ifno,
					options.calibration.tsys_ifno, polno, array_generator,
					&shared->calibration_context);
			if (shared->calibration_context.num_channel != shared->num_chan) {
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

	float *out_data[rows_per_processing];
	uint8_t *out_flag[rows_per_processing];
	float *out_tsys[rows_per_processing];

	WorkArea work_area;
};

WorkingSet *InitializeWorkingSet(SharedWorkingSet *shared, size_t elements,
		size_t num_chan, AlignedArrayGenerator *array_generator) {
	std::unique_ptr<WorkingSet[]> result(new WorkingSet[elements]);
	for (size_t i = 0; i < elements; ++i) {
		result[i].shared = shared;
		for (size_t j = 0; j < rows_per_processing; ++j) {
			result[i].spectrum[j] = array_generator->GetAlignedArray<float>(
					num_chan);
			result[i].flag[j] = array_generator->GetAlignedArray<unsigned char>(
					num_chan);
			result[i].out_data[j] = array_generator->GetAlignedArray<float>(
					num_chan);
			result[i].out_flag[j] = array_generator->GetAlignedArray<
					unsigned char>(num_chan);
			result[i].out_tsys[j] = array_generator->GetAlignedArray<float>(
					num_chan);
		}
		result[i].work_area.mask = array_generator->GetAlignedArray<bool>(
				num_chan);
		result[i].work_area.baseline_mask = array_generator->GetAlignedArray<
		bool>(num_chan);
		result[i].work_area.baseline_after_mask =
				array_generator->GetAlignedArray<bool>(num_chan);
		result[i].work_area.timestamp = array_generator->GetAlignedArray<double>(1);
		result[i].work_area.sky = array_generator->GetAlignedArray<float>(num_chan);
	}
	return result.release();
}

void Reduce(E2EOptions const &options) {
	LOG4CXX_DEBUG(logger, "Enter: Reduce");
	unsigned int const num_polarizations = 4;
	xdispatch::main_queue().sync([&options] {
		// create output table by copying input table
			CreateOutputTable(options.input_file, options.output_file);
		});
	for (unsigned int polno = 0; polno < num_polarizations; polno++) {
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
			LOG4CXX_DEBUG(logger, "Skip POLNO " << polno);
			continue;
		}

		LOG4CXX_INFO(logger, "Processing POLNO " << polno);

		size_t num_threads = std::min(shared->num_data,
				size_t(options.max_threads));
		std::unique_ptr<WorkingSet[]> working_sets(
				InitializeWorkingSet(shared, num_threads, shared->num_chan,
						&array_generator));

		std::vector<size_t> available_workers(num_threads);
		for (size_t i = 0; i < num_threads; ++i) {
			available_workers[i] = i;
		}

		LOG4CXX_DEBUG(logger, "Initialized");

		bool const serialize = options.serialize;
		xdispatch::queue serial_queue("my serial queue");
		xdispatch::semaphore semaphore(num_threads);
		volatile bool ok = true;

		auto group = xdispatch::group();
		for (size_t i = 0; ok && i < shared->num_data; i += rows_per_processing) {
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
						GetArrayCell(working_sets[working_set_id].spectrum[j], i+j, shared->input_spectra_column,
								casa::IPosition(1, shared->num_chan));
					}
					for (size_t j = 0; j < rows; ++j) {
						GetArrayCell(working_sets[working_set_id].flag[j], i+j, shared->input_flagtra_column,
								casa::IPosition(1, shared->num_chan));
					}
					GetScalarCells(working_sets[working_set_id].timestamp, i, rows, shared->input_time_column);
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
							serial_queue.sync([working_set_id, &semaphore, &available_workers] {
										available_workers.push_back(working_set_id);
									});
							semaphore.release();
						});
				ParallelJob(i, rows, shared->num_chan, t, v, f, options, shared,
						&working_sets[working_set_id].work_area,
						working_sets[working_set_id].out_data,
						working_sets[working_set_id].out_flag,
						working_sets[working_set_id].out_tsys);
				for (size_t j = 0; j < rows; ++j) {
					PutResult(const_cast<SharedWorkingSet *>(shared), i + j,
							shared->num_chan,
							working_sets[working_set_id].out_data[j],
							working_sets[working_set_id].out_flag[j],
							working_sets[working_set_id].out_tsys[j]);
					LOG4CXX_DEBUG(logger, "Row " << i + j << " have done.");
				}
			} else {
				group.async(
						[=, &ok, &working_sets, &available_workers, &serial_queue, &semaphore, &options] {
							ScopeGuard cleanup([=,&available_workers, &serial_queue, &semaphore] {
										serial_queue.sync([=, &available_workers] {
													available_workers.push_back(working_set_id);
												});
										semaphore.release();
									});
							do {
								try {
									ParallelJob(i, rows, shared->num_chan, t, v, f, options, shared,
											&working_sets[working_set_id].work_area,
											working_sets[working_set_id].out_data,
											working_sets[working_set_id].out_flag,
											working_sets[working_set_id].out_tsys);
								} catch (std::runtime_error &e) {
									LOG4CXX_ERROR(logger, "Exception was thrown:" << e.what());
									ok = false;
									break;
								}

								xdispatch::main_queue().sync([=, &shared, &working_sets]() {
											// run casa operations in main_queue.
											for (size_t j = 0; j < rows; ++j) {
												PutResult(const_cast<SharedWorkingSet *>(shared), i + j, shared->num_chan,
														working_sets[working_set_id].out_data[j],
														working_sets[working_set_id].out_flag[j],
														working_sets[working_set_id].out_tsys[j]);
												LOG4CXX_DEBUG(logger, "Row " << i + j << " have done.");
											}
										});
							}while (false);
						});
			}
		} // i
		if (!serialize) {
			group.wait(xdispatch::time_forever);
			xdispatch::main_queue().sync([] {}); // to ensure other jobs submitted to main_queue finished
		}
		if (!ok) {
			throw std::runtime_error("Failed");
		}
	} // polno
	LOG4CXX_DEBUG(logger, "Leave: Reduce");
}

void BatchReduce(E2EOptions const &options) {
	LOG4CXX_DEBUG(logger, "Enter: Reduce");
	xdispatch::main_queue().sync([&options] {
		// create output table by copying input table
			CreateOutputTable(options.input_file, options.output_file);
		});
	unsigned int const num_polarizations = 4;
	std::unique_ptr<SharedWorkingSet> shared[num_polarizations];
	struct {
		size_t num_channels;
		size_t num_rows;
		float **spectrum;
		uint8_t **flag;
		double *timestamp;
	} input_data[num_polarizations];
	struct {
		float **spectrum;
		uint8_t **flag;
		float **tsys;
	} output_data[num_polarizations];
	size_t max_num_chan = 0;
	AlignedArrayGenerator array_generator;
	{
		LOG4CXX_INFO(logger, "Start reading data...");
		auto start = sakura_GetCurrentTime();
		for (unsigned int polno = 0; polno < num_polarizations; polno++) {
			xdispatch::main_queue().sync(
					[&shared, &options, &polno, &array_generator] {
						shared[polno].reset(InitializeSharedWorkingSet(options, polno, &array_generator));
					});
			input_data[polno].spectrum = nullptr;
			input_data[polno].flag = nullptr;
			input_data[polno].timestamp = nullptr;
			input_data[polno].num_channels = 0;
			input_data[polno].num_rows = shared[polno]->num_data;
			output_data[polno].spectrum = nullptr;
			output_data[polno].flag = nullptr;
			output_data[polno].tsys = nullptr;
			if (input_data[polno].num_rows > 0) {
				input_data[polno].num_channels = shared[polno]->num_chan;
				input_data[polno].spectrum = array_generator.GetAlignedArray<
						float*>(input_data[polno].num_rows);
				input_data[polno].flag = array_generator.GetAlignedArray<
						uint8_t*>(input_data[polno].num_rows);
				input_data[polno].timestamp = array_generator.GetAlignedArray<
						double>(input_data[polno].num_rows);

				output_data[polno].spectrum = array_generator.GetAlignedArray<
						float*>(input_data[polno].num_rows);
				output_data[polno].flag = array_generator.GetAlignedArray<
						uint8_t*>(input_data[polno].num_rows);
				output_data[polno].tsys =
						array_generator.GetAlignedArray<float*>(
								input_data[polno].num_rows);

				max_num_chan = std::max(max_num_chan,
						input_data[polno].num_channels);
				xdispatch::main_queue().sync(
						[&input_data, &output_data, &shared, &options, polno, &array_generator] {
							for (size_t row = 0; row < input_data[polno].num_rows; ++row) {
								LOG4CXX_DEBUG(logger, "Reading record " << row);
								input_data[polno].spectrum[row] = array_generator.GetAlignedArray<float>(input_data[polno].num_channels);
								GetArrayCell(input_data[polno].spectrum[row], row, shared[polno]->input_spectra_column,
										casa::IPosition(1, input_data[polno].num_channels));
								input_data[polno].flag[row] = array_generator.GetAlignedArray<uint8_t>(input_data[polno].num_channels);
								GetArrayCell(input_data[polno].flag[row], row, shared[polno]->input_flagtra_column,
										casa::IPosition(1, input_data[polno].num_channels));

								output_data[polno].spectrum[row] = array_generator.GetAlignedArray<float>(input_data[polno].num_channels);
								output_data[polno].flag[row] = array_generator.GetAlignedArray<uint8_t>(input_data[polno].num_channels);
								output_data[polno].tsys[row] = array_generator.GetAlignedArray<float>(input_data[polno].num_channels);
							}
							GetScalarCells(input_data[polno].timestamp, 0, input_data[polno].num_rows, shared[polno]->input_time_column);
						});
			} else {
				LOG4CXX_DEBUG(logger, "Skip POLNO " << polno);
			}
		}
		LOG4CXX_INFO(logger,
				"All data have read within " << sakura_GetCurrentTime() - start << " sec");
	}

	LOG4CXX_INFO(logger, "Start analyzing...");
	double start = sakura_GetCurrentTime();
	size_t num_threads = options.max_threads;
	std::unique_ptr<WorkingSet[]> working_sets(
			InitializeWorkingSet(nullptr, num_threads, max_num_chan,
					&array_generator));

	std::vector<size_t> available_workers(num_threads);
	for (size_t i = 0; i < num_threads; ++i) {
		available_workers[i] = i;
	}
	bool const serialize = options.serialize;
	xdispatch::queue serial_queue("my serial queue");
	xdispatch::semaphore semaphore(num_threads);
	auto group = xdispatch::group();

	volatile bool ok = true;
	for (unsigned int polno = 0; ok && polno < num_polarizations; polno++) {

		if (input_data[polno].num_rows > 0) {
			LOG4CXX_INFO(logger, "Processing POLNO " << polno);
		}

		for (size_t i = 0; i < num_threads; ++i) {
			working_sets[i].shared = shared[polno].get();
		}

		for (size_t i = 0; ok && i < input_data[polno].num_rows; i +=
				rows_per_processing) {
			size_t rows = std::min(rows_per_processing,
					input_data[polno].num_rows - i);
			semaphore.acquire();
			size_t working_set_id = 0;
			serial_queue.sync([&working_set_id, &available_workers] {
				assert(available_workers.size() > 0);
				working_set_id = available_workers.back();
				available_workers.pop_back();
			});
			LOG4CXX_DEBUG(logger, "working_set_id: " << working_set_id);
			// get data from the table
			xdispatch::main_queue().sync([=, &shared, &working_sets] {
				// run casa operations in main_queue.
					LOG4CXX_DEBUG(logger, "Reading record " << i);
					for (size_t j = 0; j < rows; ++j) {
						GetArrayCell(working_sets[working_set_id].spectrum[j], i+j, shared[polno]->input_spectra_column,
								casa::IPosition(1, input_data[polno].num_channels));
					}
					for (size_t j = 0; j < rows; ++j) {
						GetArrayCell(working_sets[working_set_id].flag[j], i+j, shared[polno]->input_flagtra_column,
								casa::IPosition(1, input_data[polno].num_channels));
					}
					GetScalarCells(working_sets[working_set_id].timestamp, i, rows, shared[polno]->input_time_column);
				});

			auto v =
					reinterpret_cast<float const * const *>(&input_data[polno].spectrum[i]);
			auto f =
					reinterpret_cast<uint8_t const * const *>(&input_data[polno].flag[i]);
			auto t =
					reinterpret_cast<double const *>(&input_data[polno].timestamp[i]);
			if (serialize) {
				ScopeGuard cleanup(
						[=,&available_workers, &serial_queue, &semaphore] {
							serial_queue.sync([=,&working_set_id, &available_workers] {
										available_workers.push_back(working_set_id);
									});
							semaphore.release();
						});
				ParallelJob(i, rows, shared[polno]->num_chan, t, v, f, options,
						shared[polno].get(),
						&working_sets[working_set_id].work_area,
						&output_data[polno].spectrum[i],
						&output_data[polno].flag[i],
						&output_data[polno].tsys[i]);
			} else {
				group.async(
						[=, &ok, &shared, &input_data, &output_data, &working_sets, &available_workers, &serial_queue, &semaphore, &options] {
							ScopeGuard cleanup([=,&available_workers, &serial_queue, &semaphore] {
										serial_queue.sync([working_set_id, &available_workers, &semaphore] {
													available_workers.push_back(working_set_id);
												});
										semaphore.release();
									});
							try {
							ParallelJob(i, rows, input_data[polno].num_channels, t, v, f, options, shared[polno].get(),
									&working_sets[working_set_id].work_area,
									&output_data[polno].spectrum[i],
									&output_data[polno].flag[i],
									&output_data[polno].tsys[i]);
							} catch (std::runtime_error &e) {
								LOG4CXX_ERROR(logger, "Exception was thrown:" << e.what());
								ok = false;
							}
						});
			}
		} // i
	} // polno
	if (!serialize) {
		group.wait(xdispatch::time_forever);
	}
	if (!ok) {
		throw std::runtime_error("Failed");
	}
	LOG4CXX_INFO(logger,
			"Analyzed within " << sakura_GetCurrentTime() - start << " sec");
	{
		LOG4CXX_INFO(logger, "Start writing data...");
		double start = sakura_GetCurrentTime();
		xdispatch::main_queue().sync(
				[=, &shared, &input_data, &output_data, &working_sets]() {
					for (unsigned int polno = 0; polno < num_polarizations; polno++) {
						for (size_t i = 0; i < input_data[polno].num_rows; ++i) {
							PutResult(const_cast<SharedWorkingSet *>(shared[polno].get()), i, input_data[polno].num_channels,
									output_data[polno].spectrum[i],
									output_data[polno].flag[i],
									output_data[polno].tsys[i]);
							LOG4CXX_DEBUG(logger, "Row " << i << " have done.");
						}
					}
				});
		LOG4CXX_INFO(logger,
				"All data have written within " << sakura_GetCurrentTime() - start << " sec");
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

	if (options.input_file.size() > 0) {
		if (options.batch) {
			BatchReduce(options);
		} else {
			Reduce(options);
		}
	}
	LOG4CXX_DEBUG(logger, "Leave: E2eReduce");
}

void main_(int argc, char const* const argv[]) {
	LOG4CXX_INFO(logger, "start");
	double start_time = sakura_GetCurrentTime();

	sakura_Status result;
	xdispatch::main_queue().sync([&] {
		LOG4CXX_DEBUG(logger, "Initializing libsakura");
		result = sakura_Initialize(nullptr, nullptr);
	});
	if (result == sakura_Status_kOK) {
		try {
			E2eReduce(argc, argv);
		} catch (casa::AipsError &e) {
			LOG4CXX_ERROR(logger, "Exception raised: " << e.getMesg());
		} catch (std::runtime_error &e) {
			LOG4CXX_ERROR(logger, "Exception raised: " << e.what());
		} catch (...) {
			LOG4CXX_ERROR(logger, "Exception raised");
		}
		xdispatch::main_queue().sync([start_time] {
			LOG4CXX_DEBUG(logger, "Cleaning up libsakura");
			sakura_CleanUp();
			LOG4CXX_INFO(logger, "Sync-ing...");
			sync();
			double end_time = sakura_GetCurrentTime();
			std::cout << "finished " << end_time - start_time << "secs\n";
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
