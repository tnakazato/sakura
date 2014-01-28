#ifndef _SAKURA_E2E_CONTEXT_HANDLER_H_
#define _SAKURA_E2E_CONTEXT_HANDLER_H_

#include <iostream>

#include <libsakura/sakura.h>

#include "utils.h"

namespace {
// Context for calibration
struct CalibrationContext {
	size_t num_data_sky;
	size_t num_data_tsys;
	size_t num_channel;
	float *sky_spectra;
	float *tsys;
	double *timestamp_sky;
	double *timestamp_tsys;
};

void FillCalibrationContext(std::string const sky_table,
		std::string const tsys_table, unsigned int sky_ifno,
		unsigned int tsys_ifno, unsigned int polno,
		AlignedArrayGenerator *array_generator,
		CalibrationContext *calibration_context) {
	// local AlignedArrayGenerator
	AlignedArrayGenerator local_array_generator;

	// sky table
	float *sky_spectra = nullptr;
	double *sky_time = nullptr;
	size_t num_chan_sky, num_row_sky;
	GetFromCalTable(sky_table, sky_ifno, polno, "SPECTRA", array_generator,
			&sky_spectra, &sky_time, &num_chan_sky, &num_row_sky);
	assert(sky_spectra != nullptr && sky_time != nullptr);

	// tsys table
	float *tsys = nullptr;
	double *tsys_time = nullptr;
	size_t num_chan_tsys, num_row_tsys;
	GetFromCalTable(tsys_table, tsys_ifno, polno, "TSYS",
			&local_array_generator, &tsys, &tsys_time, &num_chan_tsys,
			&num_row_tsys);
	assert(tsys != nullptr && tsys_time != nullptr);
	array_generator->Transfer(local_array_generator.index() - 2,
			&local_array_generator);
	array_generator->Transfer(local_array_generator.index() - 1,
			&local_array_generator);

	// get frequency label from the table
	// for spectral data
	double *frequency_label_target = local_array_generator.GetAlignedArray<
			double>(num_chan_sky);
	GetFrequencyLabelFromScantable(sky_table, sky_ifno, num_chan_sky,
			frequency_label_target);

	// for Tsys
	double *frequency_label_tsys =
			local_array_generator.GetAlignedArray<double>(num_chan_tsys);
	GetFrequencyLabelFromScantable(tsys_table, tsys_ifno, num_chan_tsys,
			frequency_label_tsys);

	float *interpolated_tsys = array_generator->GetAlignedArray<float>(
			num_chan_sky * num_row_tsys);
	sakura_InterpolateXAxisFloat(sakura_InterpolationMethod_kSpline, 0,
			num_chan_tsys, frequency_label_tsys, num_row_tsys, tsys,
			num_chan_sky, frequency_label_target, interpolated_tsys);

	// Create Context and struct for calibration
	// calibration context
	calibration_context->num_channel = num_chan_sky;
	calibration_context->num_data_sky = num_row_sky;
	calibration_context->num_data_tsys = num_row_tsys;
	calibration_context->timestamp_sky = sky_time;
	calibration_context->timestamp_tsys = tsys_time;
	calibration_context->sky_spectra = sky_spectra;
	calibration_context->tsys = interpolated_tsys;
}

} /* anonymous namespace */

#endif /* _SAKURA_E2E_CONTEXT_HANDLER_H_ */
