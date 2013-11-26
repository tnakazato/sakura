#ifndef _SAKURA_E2E_CONTEXT_HANDLER_H_
#define _SAKURA_E2E_CONTEXT_HANDLER_H_

#include <iostream>

#include <libsakura/sakura.h>

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

// Deleter for baseline context
struct BaselineContextDeleter {
	inline void operator()(sakura_BaselineContext *context) const noexcept {
//		std::cout << "destroy baseline context" << std::endl;
		sakura_DestroyBaselineContext(context);
	}
};

// Deleter for convolve1d context
struct Convolve1DContextDeleter {
	inline void operator()(sakura_Convolve1DContext *context) const noexcept {
//		std::cout << "destroy convolve1d context" << std::endl;
		sakura_DestroyConvolve1DContext(context);
	}
};


} /* anonymous namespace */

#endif /* _SAKURA_E2E_CONTEXT_HANDLER_H_ */
