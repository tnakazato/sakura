#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "libsakura/sakura.h"
#include "libsakura/memory_manager.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/localdef.h"

#define EIGEN_DENSEBASE_PLUGIN "eigen_binary_visitor_plugin.h"
#include <Eigen/Core>

using ::Eigen::Map;
using ::Eigen::Array;
using ::Eigen::Dynamic;
using ::Eigen::Aligned;

namespace {

inline void GetBaselineModel(
		size_t num_chan, int order, double *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	size_t num_model = order + 1;

	for (size_t i = 0; i < num_model; ++i) {
		for (size_t j = 0; j < num_chan; ++j) {
			out[num_chan*i+j] = pow(static_cast<double>(j), static_cast<double>(i));
		}
	}
}

inline void DoSubtractBaseline(
		size_t num_chan, float const *in_data, bool const *in_mask, size_t num_model,
		double *model_data, float clipping_threshold_sigma, unsigned int num_clipping_max,
		bool get_residual, float *out) {
	//std::cout << "DoSubtractBaselineEigen function is called" << std::endl;

	assert(LIBSAKURA_SYMBOL(IsAligned)(in_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	bool *clip_mask = reinterpret_cast<bool *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(bool)*num_chan));
	for (size_t i = 0; i < num_chan; ++i) {
		clip_mask[i] = true;
	}

	bool *composite_mask = reinterpret_cast<bool *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(bool)*num_chan));
	float *best_fit_model = reinterpret_cast<float *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(float)*num_chan));
	float *residual_data = reinterpret_cast<float *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(float)*num_chan));
	//bool *new_clip_mask = reinterpret_cast<bool *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(bool)*num_chan));

	for (size_t i = 0; i <= num_clipping_max; ++i) {
		LIBSAKURA_SYMBOL(OperateLogicalAnd)(num_chan, in_mask, clip_mask, composite_mask);
		LIBSAKURA_SYMBOL(GetBestFitModel)(num_chan, in_data, composite_mask, num_model, model_data, best_fit_model);
		LIBSAKURA_SYMBOL(OperateFloatSubtraction)(num_chan, in_data, best_fit_model, residual_data);

		//GetRms(); // calculate rms
		//LIBSAKURA_SYMBOL(SetTrueFloatInRangesInclusive)(); //new_clip_mask generated based on rms and residual_data

		/*
		bool no_mask_change_after_clipping = true;
		for (size_t j = 0; j < num_chan; ++j) {
			if (new_clip_mask[j] != clip_mask[j]) {
				no_mask_change_after_clipping = false;
				break;
			}
		}

		if (no_mask_change_after_clipping) {
			break;
		} else {
			for (size_t j = 0; j < num_chan; ++j) {
				clip_mask[j] = new_clip_mask[j];
			}
		}
		*/
	}

	if (get_residual) {
		for (size_t i = 0; i < num_chan; ++i) {
			out[i] = residual_data[i];
		}
	} else {
		for (size_t i = 0; i < num_chan; ++i) {
			out[i] = best_fit_model[i];
		}
	}

	LIBSAKURA_PREFIX::Memory::Free(clip_mask);
	LIBSAKURA_PREFIX::Memory::Free(composite_mask);
	LIBSAKURA_PREFIX::Memory::Free(best_fit_model);
	LIBSAKURA_PREFIX::Memory::Free(residual_data);
	//LIBSAKURA_PREFIX::Memory::Free(new_clip_mask);
}

inline void SubtractBaselinePolynomial(
		size_t num_chan, float const *in_data, bool const *in_mask,
		int order, float clipping_threshold_sigma,
		unsigned int num_clipping_max, bool get_residual,
		float *out) {
	size_t num_model = order + 1;
	double *model_data = reinterpret_cast<double *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(double)*num_chan*num_model));

	GetBaselineModel(num_chan, order, model_data);
	DoSubtractBaseline(num_chan, in_data, in_mask, num_model, model_data,
			clipping_threshold_sigma, num_clipping_max, get_residual, out);

	LIBSAKURA_PREFIX::Memory::Free(model_data);
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(Baseline, ARCH_SUFFIX)::SubtractBaselinePolynomial(
		size_t num_chan, float const in_data[/*num_chan*/],
		bool const in_mask[/*num_chan*/], int order,
		float clipping_threshold_sigma, unsigned int num_clipping_max,
		bool get_residual, float out[/*num_chan*/]) const {
	SubtractBaselinePolynomial(num_chan, in_data, in_mask, order,
			clipping_threshold_sigma, num_clipping_max, get_residual, out);
}

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::GetBaselineModel(
		size_t num_chan, int order, double out[/*(order+1)*num_chan*/]) const {
	GetBaselineModel(num_chan, order, out);
}

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::DoSubtractBaseline(
		size_t num_chan, float const in_data[/*num_chan*/],
		bool const in_mask[/*num_chan*/], size_t num_model,
		double model_data[/*num_model * num_chan*/],
		float clipping_threshold_sigma, unsigned int num_clipping_max,
		bool get_residual, float out[/*num_chan*/]) const {
	DoSubtractBaseline(num_chan, in_data, in_mask, num_model, model_data,
		clipping_threshold_sigma, num_clipping_max, get_residual, out);
}

}
