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
	size_t num_channel_sky;
	size_t num_channel_tsys;
	float *sky_spectra;
	float *tsys;
	double *timestamp_sky;
	double *timestamp_tsys;
	double *frequency_label_sky;
	double *frequency_label_tsys;
};

void FillCalibrationContext(
		std::string const sky_table, std::string const tsys_table,
		unsigned int sky_ifno, unsigned int tsys_ifno,
		AlignedArrayGenerator *array_generator, CalibrationContext *calibration_context) {
	// sky table
	float *sky_spectra = nullptr;
	double *sky_time = nullptr;
	size_t num_chan_sky, num_row_sky;
	GetFromCalTable(sky_table, sky_ifno, "SPECTRA", array_generator,
			&sky_spectra, &sky_time, &num_chan_sky, &num_row_sky);
	assert(sky_spectra != nullptr && sky_time != nullptr);

	// tsys table
	float *tsys = nullptr;
	double *tsys_time = nullptr;
	size_t num_chan_tsys, num_row_tsys;
	GetFromCalTable(tsys_table, tsys_ifno, "TSYS", array_generator, &tsys,
			&tsys_time, &num_chan_tsys, &num_row_tsys);
	assert(tsys != nullptr && tsys_time != nullptr);

	// get frequency label from the table
	// for spectral data
	double *frequency_label_target = array_generator->GetAlignedArray<double>(
			num_chan_sky);
	GetFrequencyLabelFromScantable(sky_table, sky_ifno, num_chan_sky,
			frequency_label_target);

	// for Tsys
	double *frequency_label_tsys = array_generator->GetAlignedArray<double>(
			num_chan_tsys);
	GetFrequencyLabelFromScantable(tsys_table, tsys_ifno, num_chan_tsys,
			frequency_label_tsys);

	// Create Context and struct for calibration
	// calibration context
	calibration_context->num_channel_sky = num_chan_sky;
	calibration_context->num_channel_tsys = num_chan_tsys;
	calibration_context->num_data_sky = num_row_sky;
	calibration_context->num_data_tsys = num_row_tsys;
	calibration_context->timestamp_sky = sky_time;
	calibration_context->timestamp_tsys = tsys_time;
	calibration_context->sky_spectra = sky_spectra;
	calibration_context->tsys = tsys;
	calibration_context->frequency_label_sky = frequency_label_target;
	calibration_context->frequency_label_tsys = frequency_label_tsys;
}

} /* anonymous namespace */

#endif /* _SAKURA_E2E_CONTEXT_HANDLER_H_ */
